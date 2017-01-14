import sys
from sourmash_lib import signature


def add_moltype_args(parser, default_dna=None):
    parser.add_argument('--protein', dest='protein', action='store_true')
    parser.add_argument('--no-protein', dest='protein',
                        action='store_false')
    parser.set_defaults(protein=False)

    parser.add_argument('--dna', dest='dna', default=None,
                        action='store_true')
    parser.add_argument('--no-dna', dest='dna', action='store_false')
    parser.set_defaults(dna=default_dna)


def get_moltype(sig, require=False):
    if sig.estimator.is_molecule_type('dna'):
        moltype = 'DNA'
    elif sig.estimator.is_molecule_type('protein'):
        moltype = 'protein'
    else:
        raise ValueError('unknown molecule type for sig {}'.format(sig.name()))

    return moltype


class LoadSingleSignatures(object):
    def __init__(self, filelist,  select_ksize=None, select_moltype=None,
                 ignore_files=set()):
        self.filelist = filelist
        self.select_ksize = select_ksize
        self.select_moltype = select_moltype
        self.ignore_files = ignore_files

        self.skipped_ignore = 0
        self.skipped_nosig = 0
        self.ksizes = set()
        self.moltypes = set()

    def __iter__(self):
        for filename in self.filelist:
            if filename in self.ignore_files:
                self.skipped_ignore += 1
                continue

            sl = signature.load_signatures(filename,
                                           select_ksize=self.select_ksize,
                                           select_moltype=self.select_moltype)
            sl = list(sl)
            if len(sl) == 0:
                self.skipped_nosig += 1
                continue

            for query in sl:
                query_moltype = get_moltype(query)
                query_ksize = query.estimator.ksize

                self.ksizes.add(query_ksize)
                self.moltypes.add(query_moltype)

                yield filename, query, query_moltype, query_ksize

            if len(self.ksizes) > 1 or len(self.moltypes) > 1:
                raise ValueError('multiple k-mer sizes/molecule types present')

