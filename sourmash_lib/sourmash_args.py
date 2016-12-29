import sys
from sourmash_lib import signature as sig


def add_moltype_args(parser, default_dna=None):
    parser.add_argument('--protein', dest='protein', action='store_true')
    parser.add_argument('--no-protein', dest='protein',
                        action='store_false')
    parser.set_defaults(protein=False)

    parser.add_argument('--dna', dest='dna', default=None,
                        action='store_true')
    parser.add_argument('--no-dna', dest='dna', action='store_false')
    parser.set_defaults(dna=default_dna)


class LoadSingleSignatures(object):
    def __init__(self, filelist,  select_ksize=None, select_moltype=None,
                 ignore_files=set()):
        self.filelist = filelist
        self.select_ksize = select_ksize
        self.select_moltype = select_moltype
        self.ignore_files = ignore_files

        self.skipped_ignore = 0
        self.skipped_nosig = 0

    def __iter__(self):
        for filename in self.filelist:
            if filename in self.ignore_files:
                self.skipped_iignore += 1
                continue

            sl = sig.load_signatures(filename,
                                     select_ksize=self.select_ksize,
                                     select_moltype=self.select_moltype)
            sl = list(sl)
            if len(sl) != 1:
                self.skipped_nosig += 1
                continue

            query = sl[0]
            query_moltype = 'UNKNOWN'
            if query.estimator.is_molecule_type('dna'):
                query_moltype = 'DNA'
            elif query.estimator.is_molecule_type('protein'):
                query_moltype = 'protein'
            query_ksize = query.estimator.ksize

            yield filename, query, query_moltype, query_ksize
