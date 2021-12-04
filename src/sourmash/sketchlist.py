"Sketchlist code for sketching many files."
import csv
#from enum import Enum
from collections import namedtuple

from command_sketch import DEFAULTS, _parse_param_str, _signatures_for_sketch_factory


# set up preprocessing functions for column stuff
#preprocess = {}

# exact matches
#preprocess['name'] = lambda x: x
# does param string here override cli `-p`?
#preprocess['param_string'] = lambda x: x

# input and output dir/filenames (no spaces)
#preprocess['input'] = lambda x: x.split(' ')[0]
#preprocess['output'] = lambda x: x.split(' ')[0]

#preprocess['license'] = lambda x: x


### Experimenting with a way to sketch from sketchlist
# current approach:
# 1. parse sketchlist csv info --> sketchinfo
#   using SourmashSignature would NOT work bc want to allow multiple compatible sketch param strings.
# 2. parse cli info --> sketchinfo
# 3. since each sketch has a separate sketchinfo associated, build sig factories independently (is this bad, memory-wise etc?)
# 4. execute sketches

# maybe this would be better as class? default values of name, output, etc???
SketchInfo = namedtuple("SketchInfo", "name, input, output, param_string, exists, valid_params, license")


class SignatureSketchInfo:
    """Sketchlist class for sketching collections of signatures.

    Initialize using ``SignatureSketchlist.from_sketchlist_args(argstr)``.

    Here, 'sketchfile' is the path to a CSV file.

    'coltype's that are currently supported:
    * 'name' - name for the signature
    * 'input' - input filename or directory
    * 'output' - output filename or directory
    * 'license' - signature license. Currently only CC0 is supported
    * 'param_string' - parameter string for sketching

    """
    supported_coltypes = ('name', 'input', 'output',
                          'license', 'param_string')

    def __init__(self, require_all=False):
        "create a sketchlist"

        self.sketchlist_file = None
        self.sketchlist = []
        self.all_param_str = set()
        #self.all_moltype = set() # pot make work with mult moltype

    def add(self, sketchinfo):
        "Add a sketchinfo namedtuple to this sketchlist."
        if isinstance(sketchinfo, SketchInfo):
            self.sketchlist.append(sketchinfo)
            #self.all_moltype.add(sketchinfo.moltype)

def check_inputs(self, inp):
        "does the input file exist?"
        #if any([os.path.exists(x) for x in inp]): # does this do more work?
        for inpF in inp:
            if os.path.exists(inpF):
                # if any inputfile exists, we can build sketch.
                return True
        return False

    # implement so we can initialize from py api?
    #def init(self, sketchinfo=[]):
    #    "initialize a Sketchlist object with given sketch info/params."
    #    if self.sketchlist is not None:
    #        raise ValueError("already initialized?")
    #    self.sketchinfo = sketchinfo
    #    return self.sketchinfo

    # use preprocess here instead??
    def _parse_sketchinfo(row):
        name,inp,outp,param_str,exists,valid_params,license= [None]*7
        inp = row['input'].split(' ')[0].split(',') # allow comma-separated for merging files --> sigs?
        if isinstance(inp, str):
            inp = [inp]
        exists = check_inputs(inp)
        if row['name']:
            name= row['name']
        if row['license']:
            license= row['license']
        if row['output']:
            outp = row['output'].split(' ')[0]
        if row['param_string']:
            param_str = row['param_string']
            # run this just to check that param str is ok? or wait till sketch factory?
            try:
                mt, par = _parse_params_str(param_str)
                valid_params=True
            except ValueError:
                valid_params=False
        # do we want to organize by param string here, to create sigs for sketch factory just once?
        return SketchInfo(name,inp,outp,param_str,exists,valid_params,license)


    def parse_cli_to_sketchlist(self, args):
        name,inp,outp,param_str,exists,valid_params,license= [None]*7

        if args.merge:
            if not args.output:
                error("ERROR: must specify -o with --merge")
                sys.exit(-1)
            else:
                name=args.merge

        if args.output and args.outdir:
            error("ERROR: --outdir doesn't make sense with -o/--output")
            sys.exit(-1)
        elif args.output:
            outp = args.output
        else:
            # can we just handle file vs dir later, or do we need to keep track of outdir,output separately?
            outp = args.outdir

        # ntp note: inp needs to be list. I think it is, but check.
        inp =  _add_from_file_to_filenames(args)
        exists = check_inputs(inp)

        if args.license != 'CC0':
            error('error: sourmash only supports CC0-licensed signatures. sorry!')
            sys.exit(-1)
        license=args.license

        if args.param_string:
            param_str = args.param_str
            # run this just to check that param str is ok?
            try:
                mt, par = _parse_params_str(param_str)
                valid_params=True
            except ValueError:
                valid_params=False

         # do i want to do any param checking here?
         if name:
             # if these are merged sigs, add single sketchinfo
             this_sketchinfo = SketchInfo(name,inp,outp,param_str,exists,valid_params,license)
             self.add(this_sketchinfo)
         else:
             # add sketchinfo for each file
             for inputfile in inp:
                 this_sketchinfo = SketchInfo(name,[inputfile],outp,param_str,exists,valid_params,license)
                 self.add(this_sketchinfo)


    def load_sketchlist_file(self, args):
        "load sketchfile, return num empty vals, and set of duplicate vals."
        self.sketchlist_file = args.sketchlist

        #dup_vals = set() #check for duplicates?
        with open(sketchlist_file, newline='') as csvfile:
            x = csvfile.readline()

            # skip leading comment line in case there's a header
            if x[0] == '#':
                pass
            else:
                csvfile.seek(0)

            r = csv.DictReader(csvfile)

            # make sure input column exists (or process as filelist, but eh, seems better to force colnames to exist)
            if 'input' not in r.fieldnames:
                raise ValueError("please check that the sketchlist includes an 'input' column of filenames or directories")

            for row in r:
                # parse sketch info
                this_sketchinfo = _parse_sketchinfo(row)
               # look for duplicate values or empty values
                #if this_sketchinfo in sketchlist:
                #    dup_vals.add(this_sketchinfo)
                #else:
                    #self.add(this_sketchinfo)
                self.add(this_sketchinfo)

        # do we need to return the duplicate vals?
        #return n_input_missing, dup_vals


        # could probably handle mult moltypes. ignore for now and only sketch cli moltype
        #def _sketchlist_signatures_for_sketch_factory(self, moltype, mult_ksize_by_3):
        #    _signatures_for_sketch_factory(self.all_param_str, moltype, mult_ksize_by_3)

        def _execute_sketches(self, args, mult_ksize_by_3=False):
            "Run 'sketch' for all params and inputs in sketchlist."
            set_quiet(args.quiet)

            if not self.sketchlist:
                    error('error: no sketch info provided! nothing to do - exiting.')
                    sys.exit(-1)

            # get number of output sigs:
            num_sigs = len(signatures_factory.params_list)
            notify(f'Computing a total of {num_sigs} signature(s).')

            if num_sigs == 0:
                error('...nothing to calculate!? Exiting!')
                sys.exit(-1)
            # keep track of some info --> notify of these after sketching!
            n_input_missing = 0
            n_invalid_params = 0
            n_invalid_license = 0

            for n, sinfo in enumerate(self.sketchlist):
                if not sinfo.exists:
                    n_input_missing += 1
                if not sinfo.valid_params:
                    n_invalid_params += 1
                if sinfo.license != 'CC0'
                    notify('sourmash only supports CC0-licensed signatures. Skipping {n_invalid_license + 1}th invalid signature. sorry!')
                    n_invalid_license += 1
                    continue
                # how often to notify?
                if n % 1000 = 0:
                    notify(f'computing signatures for files: {", ".join(sinfo.inp)}')

                # use args.param string if not provided individually. Do we need to parse this again here, or ...?
                if not sinfo.param_string:
                    sinfo.param_string = args.param_str
                # parse the param string
                sinfo.moltype, sinfo.param_string = _parse_params_str(sinfo.param_string)
                # what to do about mult_ksize_by_3?

                # we've already separated everything out - any files in the same sketchinfo should be merged into same sketch

                try:
                    sig_factory = _signatures_for_sketch_factory(sinfo.param_str, moltype, mult_ksize_by_3)
                except ValueError as e:
                    error(f"Error creating signatures: {str(e)}")
                    continue # ignore error and keep going with remaining sketches?
                    #sys.exit(-1)
                #_compute_merged(args, sig_factory)
                # NEED TO PASS sketchinfo, not args to get proper sketch.HMM
                _compute_merged(sinfo, sig_factory)

                # can we build /add to a manifest as we go?



