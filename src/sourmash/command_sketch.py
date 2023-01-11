"""
Functions implementing the 'sketch' subcommands and related functions.
"""
import sys
import os
from collections import defaultdict, Counter
import csv
import shlex

import screed

import sourmash
from .signature import SourmashSignature
from .logging import notify, error, set_quiet, print_results
from .command_compute import (_compute_individual, _compute_merged,
                              ComputeParameters, add_seq, set_sig_name,
                              DEFAULT_MMHASH_SEED)
from sourmash import sourmash_args
from sourmash.sourmash_args import check_scaled_bounds, check_num_bounds
from sourmash.sig.__main__ import _summarize_manifest, _SketchInfo
from sourmash.manifest import CollectionManifest

DEFAULTS = dict(
    dna='k=31,scaled=1000,noabund',
    protein='k=10,scaled=200,noabund',
    dayhoff='k=16,scaled=200,noabund',
    hp='k=42,scaled=200,noabund'
)


def _parse_params_str(params_str):
    "Parse a parameter string of the form 'k=ks,num=num,scaled=scaled,abund'."
    moltype = None
    params = {}
    params['ksize'] = []
    items = params_str.split(',')
    for item in items:
        if item == 'abund':
            params['track_abundance'] = True
        elif item == 'noabund':
            params['track_abundance'] = False
        elif item.startswith('k'):
            if len(item) < 3 or item[1] != '=':
                raise ValueError("k takes a parameter, e.g. 'k=31'")
            params['ksize'].append(int(item[2:]))
        elif item.startswith('num'):
            if len(item) < 5 or item[3] != '=':
                raise ValueError("num takes a parameter, e.g. 'num=500'")
            if params.get('scaled'):
                raise ValueError("cannot set both num and scaled in a single minhash")
            try:
                num = item[4:]
                num = int(num)
            except ValueError:
                raise ValueError(f"cannot parse num='{num}' as a number")

            num = check_num_bounds(num)

            params['num'] = int(item[4:])
            params['scaled'] = 0
        elif item.startswith('scaled'):
            if len(item) < 8 or item[6] != '=':
                raise ValueError("scaled takes a parameter, e.g. 'scaled=1000'")
            if params.get('num'):
                raise ValueError("cannot set both num and scaled in a single minhash")
            try:
                scaled = item[7:]
                scaled = int(scaled)
            except ValueError:
                raise ValueError(f"cannot parse scaled='{scaled}' as an integer")

            scaled = check_scaled_bounds(scaled)

            params['scaled'] = scaled
            params['num'] = 0
        elif item.startswith('seed'):
            if len(item) < 6 or item[4] != '=':
                raise ValueError("seed takes a parameter, e.g. 'seed=42'")
            params['seed'] = int(item[5:])
        elif item in ('protein', 'dayhoff', 'hp', 'dna'):
            moltype = item
        else:
            raise ValueError(f"unknown component '{item}' in params string")

    return moltype, params


class _signatures_for_sketch_factory:
    "Build sigs on demand, based on args input to 'sketch'."
    def __init__(self, params_str_list, default_moltype):
        # first, set up defaults per-moltype
        defaults = {}
        for moltype, pstr in DEFAULTS.items():
            mt, d = _parse_params_str(pstr)
            assert mt is None             # defaults cannot have moltype set!
            defaults[moltype] = d
        self.defaults = defaults

        # next, fill out params_list
        self.params_list = []
        self.mult_ksize_by_3 = True

        if params_str_list:
            # parse each params_str passed in, using default_moltype if none
            # provided.
            for params_str in params_str_list:
                moltype, params = _parse_params_str(params_str)
                if moltype and moltype != 'dna' and default_moltype == 'dna':
                    raise ValueError(f"Incompatible sketch type ({default_moltype}) and parameter override ({moltype}) in '{params_str}'; maybe use 'sketch translate'?")
                elif moltype == 'dna' and default_moltype and default_moltype != 'dna':
                    raise ValueError(f"Incompatible sketch type ({default_moltype}) and parameter override ({moltype}) in '{params_str}'")
                elif moltype is None:
                    if default_moltype is None:
                        raise ValueError(f"No default moltype and none specified in param string")
                    moltype = default_moltype

                self.params_list.append((moltype, params))
        else:
            if default_moltype is None:
                raise ValueError(f"No default moltype and none specified in param string")
            # no params str? default to a single sig, using default_moltype.
            self.params_list.append((default_moltype, {}))

    def get_compute_params(self, *, split_ksizes=False):
        for moltype, params_d in self.params_list:
            # get defaults for this moltype from self.defaults:
            default_params = self.defaults[moltype]
            def_seed = default_params.get('seed', DEFAULT_MMHASH_SEED)
            def_num = default_params.get('num', 0)
            def_abund = default_params['track_abundance']
            def_scaled = default_params.get('scaled', 0)
            def_dna = default_params.get('is_dna', moltype == 'dna')
            def_protein = default_params.get('is_protein',
                                             moltype == 'protein')
            def_dayhoff = default_params.get('is_dayhoff',
                                             moltype == 'dayhoff')
            def_hp = default_params.get('is_hp', moltype == 'hp')

            # handle ksize specially, for now - multiply by three?
            def_ksizes = default_params['ksize']
            ksizes = params_d.get('ksize')
            if not ksizes:
                ksizes = def_ksizes

            # 'command sketch' adjusts k-mer sizes by 3 if non-DNA sketch.
            if self.mult_ksize_by_3 and not def_dna:
                ksizes = [ k*3 for k in ksizes ]

            make_param = lambda ksizes: ComputeParameters(
                                            ksizes=ksizes,
                                            seed=params_d.get('seed', def_seed),
                                            protein=def_protein,
                                            dayhoff=def_dayhoff,
                                            hp=def_hp,
                                            dna=def_dna,
                                            num_hashes=params_d.get('num', def_num),
                                            track_abundance=params_d.get('track_abundance',
                                                                         def_abund),
                                            scaled=params_d.get('scaled', def_scaled))

            if split_ksizes:
                for ksize in ksizes:
                    params_obj = make_param([ksize])
                    yield params_obj
            else:
                params_obj = make_param(ksizes)
                yield params_obj

    def __call__(self, *, split_ksizes=False):
        "Produce a new set of signatures built to match the param strings."
        sigs = []
        for params in self.get_compute_params(split_ksizes=split_ksizes):
            sig = SourmashSignature.from_params(params)
            sigs.append(sig)

        return sigs


def _add_from_file_to_filenames(args):
    "Add filenames from --from-file to args.filenames"
    from .sourmash_args import load_pathlist_from_file
    if args.from_file:
        file_list = load_pathlist_from_file(args.from_file)
        args.filenames.extend(file_list)


def _execute_sketch(args, signatures_factory):
    "Once configured, run 'sketch' the same way underneath."
    set_quiet(args.quiet)

    if not args.filenames:
        error('error: no input filenames provided! nothing to do - exiting.')
        sys.exit(-1)

    if args.license != 'CC0':
        error('error: sourmash only supports CC0-licensed signatures. sorry!')
        sys.exit(-1)

    notify(f'computing signatures for files: {", ".join(args.filenames)}')

    if args.merge and not args.output:
        error("ERROR: must specify -o with --merge")
        sys.exit(-1)

    if args.output and args.output_dir:
        error("ERROR: --output-dir doesn't make sense with -o/--output")
        sys.exit(-1)

    # get number of output sigs:
    num_sigs = len(signatures_factory.params_list)
    notify(f'Computing a total of {num_sigs} signature(s) for each input.')

    if num_sigs == 0:
        error('...nothing to calculate!? Exiting!')
        sys.exit(-1)

    if args.merge:               # single name specified - combine all
        _compute_merged(args, signatures_factory)
    else:                        # compute individual signatures
        _compute_individual(args, signatures_factory)


def dna(args):
    """Compute a DNA signature for one or more files.

    CTB: make usable via Python?
    """
    # for dna:
    args.input_is_protein = False

    try:
        signatures_factory = _signatures_for_sketch_factory(args.param_string,
                                                            'dna')
    except ValueError as e:
        error(f"Error creating signatures: {str(e)}")
        sys.exit(-1)

    _add_from_file_to_filenames(args)
    _execute_sketch(args, signatures_factory)


def protein(args):
    """Compute a protein signature for one or more files.

    CTB: make usable via Python?
    """
    # for protein:
    args.input_is_protein = True
    args.check_sequence = False

    # provide good defaults for dayhoff/hp/protein!
    if args.dayhoff and args.hp:
        raise ValueError("cannot set both --dayhoff and --hp")
    if args.dayhoff:
        moltype = 'dayhoff'
    elif args.hp:
        moltype = 'hp'
    else:
        moltype = 'protein'

    try:
        signatures_factory = _signatures_for_sketch_factory(args.param_string,
                                                            moltype)
    except ValueError as e:
        error(f"Error creating signatures: {str(e)}")
        sys.exit(-1)

    _add_from_file_to_filenames(args)
    _execute_sketch(args, signatures_factory)


def translate(args):
    """Compute protein signatures from DNA/RNA, for one or more files.

    CTB: make usable via Python?
    """
    # for translate:
    args.input_is_protein = False

    # provide good defaults for dayhoff/hp/protein!
    if args.dayhoff and args.hp:
        raise ValueError("cannot set both --dayhoff and --hp")
    if args.dayhoff:
        moltype = 'dayhoff'
    elif args.hp:
        moltype = 'hp'
    else:
        moltype = 'protein'

    try:
        signatures_factory = _signatures_for_sketch_factory(args.param_string,
                                                            moltype)
    except ValueError as e:
        error(f"Error creating signatures: {str(e)}")
        sys.exit(-1)

    _add_from_file_to_filenames(args)
    _execute_sketch(args, signatures_factory)


def _compute_sigs(to_build, output, *, check_sequence=False):
    "actually build the signatures in 'to_build' and output them to 'output'"
    save_sigs = sourmash_args.SaveSignaturesToLocation(output)
    save_sigs.open()

    for (name, filename), param_objs in to_build.items():
        assert param_objs

        # now, set up to iterate over sequences.
        with screed.open(filename) as screed_iter:
            if not screed_iter:
                error(f"ERROR: no sequences found in '{filename}'?!")
                sys.exit(-1)

            # build the set of empty sigs
            sigs = []

            is_dna = param_objs[0].dna
            for p in param_objs:
                if p.dna: assert is_dna
                sig = SourmashSignature.from_params(p)
                sigs.append(sig)

            input_is_protein = not is_dna

            # read sequence records & sketch
            notify(f'... reading sequences from {filename}')
            for n, record in enumerate(screed_iter):
                if n % 10000 == 0:
                    if n:
                        notify('\r...{} {}', filename, n, end='')

                try:
                    add_seq(sigs, record.sequence, input_is_protein,
                            check_sequence)
                except ValueError as exc:
                    error(f"ERROR when reading from '{filename}' - ")
                    error(str(exc))
                    sys.exit(-1)

            notify('...{} {} sequences', filename, n, end='')

            set_sig_name(sigs, filename, name)
            for sig in sigs:
                save_sigs.add(sig)

            notify(f'calculated {len(sigs)} signatures for {n+1} sequences in {filename}')


    save_sigs.close()
    notify(f"saved {len(save_sigs)} signature(s) to '{save_sigs.location}'. Note: signature license is CC0.")


def _output_csv_info(filename, sigs_to_build):
    "output information about what signatures to build, in CSV format"
    output_n = 0
    with sourmash_args.FileOutputCSV(filename) as csv_fp:
        w = csv.DictWriter(csv_fp, fieldnames=['filename', 'sketchtype',
                                               'output_index', 'name',
                                               'param_strs'])
        w.writeheader()

        output_n = 0
        for (name, filename), param_objs in sigs_to_build.items():
            param_strs = []

            # should all be the same!
            if param_objs[0].dna:
                assert all( ( p.dna for p in param_objs ) )
                sketchtype = "dna"
            else:
                assert not any( ( p.dna for p in param_objs ) )
                sketchtype = "protein"

            for p in param_objs:
                param_strs.append(p.to_param_str())

            row = dict(filename=filename, sketchtype=sketchtype,
                       param_strs="-p " + " -p ".join(param_strs),
                       name=name, output_index=output_n)

            w.writerow(row)

            output_n += 1


def fromfile(args):
    if args.license != 'CC0':
        error('error: sourmash only supports CC0-licensed signatures. sorry!')
        sys.exit(-1)

    if args.output_signatures and os.path.exists(args.output_signatures):
        if not args.force_output_already_exists:
            error(f"** ERROR: output location '{args.output_signatures}' already exists!")
            error(f"** Not overwriting/appending.")
            error(f"** Use --force-output-already-exists if you want to overwrite/append.")
            sys.exit(-1)

    # now, create the set of desired sketch specs.
    try:
        # omit a default moltype - must be provided in param string.
        sig_factory = _signatures_for_sketch_factory(args.param_string, None)
    except ValueError as e:
        error(f"Error creating signatures: {str(e)}")
        sys.exit(-1)

    # take the signatures factory => convert into a bunch of ComputeParameters
    # objects.
    build_params = list(sig_factory.get_compute_params(split_ksizes=True))

    # confirm that they do not adjust seed, which is not supported in
    # 'fromfile' b/c we don't store that info in manifests. (see #1849)
    for p in build_params:
        if p.seed != DEFAULT_MMHASH_SEED:
            error("** ERROR: cannot set 'seed' in 'sketch fromfile'")
            sys.exit(-1)

    # cross-product all of the names in the input CSV file
    # with the sketch spec(s) provided on the command line.

    to_build = defaultdict(list)
    all_names = {}
    total_rows = 0
    skipped_sigs = 0
    n_missing_name = 0
    n_duplicate_name = 0

    for csvfile in args.csvs:
        with sourmash_args.FileInputCSV(csvfile) as r:
            for row in r:
                name = row['name']
                if not name:
                    n_missing_name += 1
                    continue

                genome = row['genome_filename']
                proteome = row['protein_filename']
                total_rows += 1

                if name in all_names:
                    n_duplicate_name += 1
                else:
                    all_names[name] = (genome, proteome)

    fail_exit = False
    if n_duplicate_name:
        error(f"** ERROR: {n_duplicate_name} entries have duplicate 'name' records. Exiting!")
        fail_exit = True

    if n_missing_name:
        error(f"** ERROR: {n_missing_name} entries have blank 'name's? Exiting!")
        fail_exit = True

    if fail_exit:
        sys.exit(-1)

    # load manifests from '--already-done' databases => turn into
    # ComputeParameters objects, indexed by name.

    already_done = defaultdict(list)
    already_done_rows = []
    for filename in args.already_done:
        idx = sourmash.load_file_as_index(filename)
        manifest = idx.manifest
        assert manifest

        # for each manifest row,
        for row in manifest.rows:
            name = row['name']
            if name:
                # build a ComputeParameters object for later comparison
                p = ComputeParameters.from_manifest_row(row)

                # add to list for this name
                already_done[name].append(p)

                # matching name? check if we already have sig. if so, store!
                if name in all_names:
                    if p in build_params:
                        already_done_rows.append(row)

    already_done_manifest = CollectionManifest(already_done_rows)
    if args.already_done:
        notify(f"Loaded {len(already_done)} pre-existing names from manifest(s)")
        notify(f"collected {len(already_done_rows)} rows for already-done signatures.")

    ## now check which are already done and track only those that are
    ## need to be done.

    total_sigs = 0
    missing = defaultdict(list)
    missing_count = 0
    for name, (genome, proteome) in all_names.items():
        plist = already_done.get(name, [])

        # check list of already done against build parameters
        for p in build_params:
            total_sigs += 1

            # does this signature already exist?
            if p not in plist:
                # nope - figure out genome/proteome needed
                filename = genome if p.dna else proteome
                filetype = 'genome' if p.dna else 'proteome'

                if filename:
                    # add to build list
                    to_build[(name, filename)].append(p)
                else:
                    notify(f"WARNING: fromfile entry '{name}' is missing a {filetype}")
                    missing[name].append(p)
                    missing_count += 1
            else:
                skipped_sigs += 1

    ## we now have 'to_build' which contains the things we can build,
    ## and 'missing', which contains anything we cannot build. Report!

    notify(f"Read {total_rows} rows, requesting that {total_sigs} signatures be built.")

    if already_done_manifest:
        info_d = _summarize_manifest(already_done_manifest)
        print_results('---')
        print_results("summary of already-done sketches:")

        for ski in info_d['sketch_info']:
            mh_type = f"num={ski['num']}" if ski['num'] else f"scaled={ski['scaled']}"
            mh_abund = ", abund" if ski['abund'] else ""

            sketch_str = f"{ski['count']} sketches with {ski['moltype']}, k={ski['ksize']}, {mh_type}{mh_abund}"

            print_results(f"   {sketch_str: <50} {ski['n_hashes']} total hashes")

        print_results('---')

    if args.output_manifest_matching:
        already_done_manifest.write_to_filename(args.output_manifest_matching)
        notify(f"output {len(already_done_manifest)} already-done signatures to '{args.output_manifest_matching}' in manifest format.")

    if missing:
        error("** ERROR: we cannot build some of the requested signatures.")
        error(f"** {missing_count} total signatures (for {len(missing)} names) cannot be built.")
        if args.ignore_missing:
            error("** (continuing past this error because --ignore-missing was set)")
        else:
            sys.exit(-1)

    notify(f"** {total_sigs - skipped_sigs} new signatures to build from {len(to_build)} files;")
    if not to_build:
        notify(f"** Nothing to build. Exiting!")
        sys.exit(0)

    if skipped_sigs:
        notify(f"** {skipped_sigs} already exist, so skipping those.")
    else:
        notify(f"** we found no pre-existing signatures that match.")

    ## first, print out a summary of to_build:

    print_results('---')
    print_results("summary of sketches to build:")

    counter = Counter()
    build_info_d = {}
    for filename, param_objs in to_build.items():
        for p in param_objs:
            moltype = p.moltype
            assert len(p.ksizes) == 1
            ksize = p.ksizes[0]
            if not p.dna: ksize //= 3

            ski = _SketchInfo(ksize=ksize, moltype=p.moltype,
                              scaled=p.scaled, num=p.num_hashes,
                              abund=p.track_abundance)
            counter[ski] += 1

    for ski, count in counter.items():
        mh_type = f"num={ski.num}" if ski.num else f"scaled={ski.scaled}"
        mh_abund = ", abund" if ski.abund else ""

        sketch_str = f"{count} sketches with {ski.moltype}, k={ski.ksize}, {mh_type}{mh_abund}"

        print_results(f"   {sketch_str: <50}")

    print_results('---')

    ## now, onward ho - do we build anything, or output stuff, or just exit?

    if args.output_signatures:                   # actually compute
        _compute_sigs(to_build, args.output_signatures,
                      check_sequence=args.check_sequence)

    if args.output_csv_info: # output info necessary to construct
        _output_csv_info(args.output_csv_info, to_build)

    notify(f"** {total_sigs} total requested; output {total_sigs - skipped_sigs}, skipped {skipped_sigs}")
