"""
Functions implementing the 'sketch' subcommands and related functions.
"""
import sys
import os
from collections import defaultdict
import csv
import shlex

import screed

import sourmash
from .signature import SourmashSignature
from .logging import notify, error, set_quiet
from .command_compute import (_compute_individual, _compute_merged,
                              ComputeParameters, add_seq, set_sig_name)
from sourmash import sourmash_args
from sourmash.sourmash_args import check_scaled_bounds, check_num_bounds

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


class _signatures_for_sketch_factory(object):
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
            def_seed = default_params.get('seed', 42)
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

            make_param = lambda ksizes: ComputeParameters(ksizes,
                                            params_d.get('seed', def_seed),
                                            def_protein,
                                            def_dayhoff,
                                            def_hp,
                                            def_dna,
                                            params_d.get('num', def_num),
                                            params_d.get('track_abundance',
                                                         def_abund),
                                            params_d.get('scaled', def_scaled))

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
                notify(f"no sequences found in '{filename}'?!")
                continue

            # build the set of empty sigs
            sigs = []

            is_dna = param_objs[0].dna
            for p in param_objs:
                if p.dna: assert is_dna
                sig = SourmashSignature.from_params(p)
                sigs.append(sig)

            input_is_protein = not is_dna

            # read sequence records & sketch
            notify('... reading sequences from {}', filename)
            for n, record in enumerate(screed_iter):
                if n % 10000 == 0:
                    if n:
                        notify('\r...{} {}', filename, n, end='')

                add_seq(sigs, record.sequence, input_is_protein,
                        check_sequence)

            notify('...{} {} sequences', filename, n, end='')

            set_sig_name(sigs, filename, name)
            for sig in sigs:
                save_sigs.add(sig)

            notify(f'calculated {len(sigs)} signatures for {n+1} sequences in {filename}')


    save_sigs.close()
    notify(f"saved {len(save_sigs)} signature(s) to '{save_sigs.location}'. Note: signature license is CC0.")


def fromfile(args):
    # TODO:
    # check license
    # check-sequence
    if args.output_signatures and args.output_commands:
        error(f"** ERROR: --output-signatures and --output-commands cannot both be specified")
        sys.exit(-1)

    if args.output_signatures and os.path.exists(args.output_signatures):
        if not args.force_output_already_exists:
            error(f"** ERROR: output location '{args.output_signatures}' already exists!")
            error(f"** Not overwriting/appending.")
            error(f"** Use --force-output-already-exists if you want to overwrite/append.")
            sys.exit(-1)

    # load manifests from '--already-done' databases => turn into
    # ComputeParameters objects, indexed by name.
    #
    # CTB: note: 'seed' is not tracked by manifests currently. Oops.
    # so we'll have to block 'seed' from being passed in by '-p'.

    already_done = defaultdict(list)
    for filename in args.already_done:
        idx = sourmash.load_file_as_index(filename)
        manifest = idx.manifest
        assert manifest

        # for each manifest row,
        for row in manifest.rows:
            name = row['name']
            if not name:
                continue

            # build a ComputeParameters object for later comparison
            p = ComputeParameters.from_manifest_row(row)

            # add to list for this name
            already_done[name].append(p)

    if args.already_done:
        notify(f"Loaded {len(already_done)} pre-existing names from manifest(s)")

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

    #
    # the big loop - cross-product all of the names in the input CSV file
    # with the sketch spec(s) provided on the command line, figure out
    # which ones do not yet exist, and record them to be calculated.
    #

    to_build = defaultdict(list)
    with open(args.csv, newline="") as fp:
        r = csv.DictReader(fp)

        n_skipped = 0
        total = 0
        for row in r:
            name = row['name']
            genome = row['genome_filename']
            proteome = row['protein_filename']

            plist = already_done[name]
            for p in build_params:
                total += 1

                # has this been done?
                if p not in plist:
                    # nope - figure out genome/proteome needed
                    filename = None
                    if p.dna:
                        filename = genome
                    else:
                        assert p.dayhoff or p.protein or p.hp
                        filename = proteome
                    assert filename

                    # add to build list
                    to_build[(name, filename)].append(p)
                else:
                    n_skipped += 1

    if to_build:
        if args.output_signatures:                   # actually compute
            notify(f"** Building {total - n_skipped} sketches for {len(to_build)} files")
            _compute_sigs(to_build, args.output_signatures)

        elif args.output_commands: # output sourmash commands
            out_obj = sourmash_args.FileOutput(args.output_commands)
            out_fp = out_obj.open()

            output_n = 0
            for (name, filename), param_objs in to_build.items():
                param_strs = []
                is_dna = None

                for p in param_objs:
                    # support straight up sourmash command output,
                    # and also CSV output
                    if p.dna:
                        assert is_dna != False
                        is_dna = True
                    else:
                        is_dna = False
                        assert is_dna != True
                    param_strs.append(p.to_param_str())

                assert is_dna is not None
                sketchtype = "dna" if is_dna else "protein"

                param_strs = "-p " + "-p ".join(param_strs)
                name = shlex.quote(name)
                output_loc = f"XXX_{output_n}.zip"

                print(f"sourmash sketch {sketchtype} {filename} --name {name} -o {output_loc} {param_strs}", file=out_fp)
                output_n += 1

            out_obj.close()


    notify(f"** {total} total requested; built {total - n_skipped}, skipped {n_skipped}")
