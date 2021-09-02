"""
Functions implementing the 'sketch' subcommands and related functions.
"""
import sys

from .signature import SourmashSignature
from .logging import notify, error, set_quiet
from .command_compute import (_compute_individual, _compute_merged,
                              ComputeParameters)
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
        elif item.startswith('k='):
            params['ksize'].append(int(item[2:]))
        elif item.startswith('num='):
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
        elif item.startswith('scaled='):
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
        elif item.startswith('seed='):
            params['seed'] = int(item[5:])
        elif item == 'protein':
            moltype = 'protein'
        elif item == 'dayhoff':
            moltype = 'dayhoff'
        elif item == 'hp':
            moltype = 'hp'
        elif item == 'dna':
            moltype = 'dna'
        else:
            raise ValueError(f"unknown component '{item}' in params string")

    return moltype, params


class _signatures_for_sketch_factory(object):
    "Build sigs on demand, based on args input to 'sketch'."
    def __init__(self, params_str_list, default_moltype, mult_ksize_by_3):

        # first, set up defaults per-moltype
        defaults = {}
        for moltype, pstr in DEFAULTS.items():
            mt, d = _parse_params_str(pstr)
            assert mt is None             # defaults cannot have moltype set!
            defaults[moltype] = d
        self.defaults = defaults

        # next, fill out params_list
        self.params_list = []
        self.mult_ksize_by_3 = mult_ksize_by_3

        if params_str_list:
            # parse each params_str passed in, using default_moltype if none
            # provided.
            for params_str in params_str_list:
                moltype, params = _parse_params_str(params_str)
                if moltype and moltype != 'dna' and default_moltype == 'dna':
                    raise ValueError(f"Incompatible sketch type ({default_moltype}) and parameter override ({moltype}) in '{params_str}'; maybe use 'sketch translate'?")
                elif moltype == 'dna' and default_moltype != 'dna':
                    raise ValueError(f"Incompatible sketch type ({default_moltype}) and parameter override ({moltype}) in '{params_str}'")
                elif moltype is None:
                    moltype = default_moltype

                self.params_list.append((moltype, params))
        else:
            # no params str? default to a single sig, using default_moltype.
            self.params_list.append((default_moltype, {}))

    def get_compute_params(self):
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

            if self.mult_ksize_by_3:
                ksizes = [ k*3 for k in ksizes ]

            params_obj = ComputeParameters(ksizes,
                                           params_d.get('seed', def_seed),
                                           def_protein,
                                           def_dayhoff,
                                           def_hp,
                                           def_dna,
                                           params_d.get('num', def_num),
                                           params_d.get('track_abundance',
                                                        def_abund),
                                           params_d.get('scaled', def_scaled))

            yield params_obj

    def __call__(self):
        "Produce a new set of signatures built to match the param strings."
        sigs = []
        for params in self.get_compute_params():
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

    if args.output and args.outdir:
        error("ERROR: --outdir doesn't make sense with -o/--output")
        sys.exit(-1)

    # get number of output sigs:
    num_sigs = len(signatures_factory.params_list)
    notify(f'Computing a total of {num_sigs} signature(s).')

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
                                                            'dna',
                                                            mult_ksize_by_3=False)
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
                                                            moltype,
                                                            mult_ksize_by_3=True)
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
                                                            moltype,
                                                            mult_ksize_by_3=True)
    except ValueError as e:
        error(f"Error creating signatures: {str(e)}")
        sys.exit(-1)

    _add_from_file_to_filenames(args)
    _execute_sketch(args, signatures_factory)
