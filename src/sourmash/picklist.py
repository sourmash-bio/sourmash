"Picklist code for extracting subsets of signatures."
import csv
from enum import Enum

# set up preprocessing functions for column stuff
preprocess = {}

# exact matches
preprocess['name'] = lambda x: x
preprocess['md5'] = lambda x: x

# identifier matches/prefix foo - space delimited identifiers
preprocess['identprefix'] = lambda x: x.split(' ')[0].split('.')[0]
preprocess['ident'] = lambda x: x.split(' ')[0]

# match 8 characters
preprocess['md5prefix8'] = lambda x: x[:8]
preprocess['md5short'] = lambda x: x[:8]


class PickStyle(Enum):
    INCLUDE = 1
    EXCLUDE = 2


class SignaturePicklist:
    """Picklist class for subsetting collections of signatures.

    Initialize using ``SignaturePicklist.from_picklist_args(argstr)``,
    which takes an argument str like so: 'pickfile:column:coltype'.

    Here, 'pickfile' is the path to a CSV file; 'column' is the name of
    the column to select from the CSV file; and 'coltype' is the type of
    matching to do on that column.

    'coltype's that are currently supported:
    * 'name' - exact match to signature's name
    * 'md5' - exact match to signature's md5sum
    * 'md5prefix8' - match to 8-character prefix of signature's md5sum
    * 'md5short' - same as md5prefix8
    * 'ident' - exact match to signature's identifier
    * 'identprefix' - match to signature's identifier, before '.'

    Identifiers are constructed by using the first space delimited word in
    the signature name.

    You can also use 'gather', 'prefetch', 'search' and 'manifest' as
    column types; these take the CSV output of 'gather', 'prefetch',
    'search', and 'sig manifest' as picklists. 'column' must be left
    blank in this case: e.g. use 'pickfile.csv::gather'.
    """
    meta_coltypes = ('manifest', 'gather', 'prefetch', 'search')
    supported_coltypes = ('md5', 'md5prefix8', 'md5short',
                          'name', 'ident', 'identprefix')

    def __init__(self, coltype, *, pickfile=None, column_name=None,
                 pickstyle=PickStyle.INCLUDE):
        "create a picklist of column type 'coltype'."

        # first, check coltype...
        valid_coltypes = set(self.meta_coltypes)
        valid_coltypes.update(self.supported_coltypes)
        if coltype not in valid_coltypes:
            raise ValueError(f"invalid picklist column type '{coltype}'")

        # if we're using gather or prefetch or manifest, set column_name
        # automatically (after checks).
        if coltype in self.meta_coltypes:
            if column_name:
                raise ValueError(f"no column name allowed for coltype '{coltype}'")
            if coltype == 'gather':
                # for now, override => md5short in column md5
                coltype = 'md5prefix8'
                column_name = 'md5'
            elif coltype == 'prefetch':
                # for now, override => md5short in column match_md5
                coltype = 'md5prefix8'
                column_name = 'match_md5'
            elif coltype == 'manifest' or coltype == 'search':
                # for now, override => md5
                coltype = 'md5'
                column_name = 'md5'
            else:               # should never be reached!
                assert 0

        self.coltype = coltype
        self.pickfile = pickfile
        self.column_name = column_name
        self.pickstyle = pickstyle

        self.preprocess_fn = preprocess[coltype]
        self.pickset = None
        self.found = set()
        self.n_queries = 0

    @classmethod
    def from_picklist_args(cls, argstr):
        "load a picklist from an argument string 'pickfile:col:coltype:style'"
        picklist = argstr.split(':')
        pickstyle = PickStyle.INCLUDE

        # pickstyle specified?
        if len(picklist) == 4:
            pickstyle_str = picklist.pop()
            if pickstyle_str == 'include':
                pickstyle = PickStyle.INCLUDE
            elif pickstyle_str == 'exclude':
                pickstyle = PickStyle.EXCLUDE
            else:
                raise ValueError(f"invalid picklist 'pickstyle' argument, '{pickstyle_str}': must be 'include' or 'exclude'")

        if len(picklist) != 3:
            raise ValueError(f"invalid picklist argument '{argstr}'")

        assert len(picklist) == 3
        pickfile, column, coltype = picklist

        return cls(coltype, pickfile=pickfile, column_name=column,
                   pickstyle=pickstyle)

    def _get_sig_attribute(self, ss):
        "for a given SourmashSignature, return attribute for this picklist."
        coltype = self.coltype
        if coltype in ('md5', 'md5prefix8', 'md5short'):
            q = ss.md5sum()
        elif coltype in ('name', 'ident', 'identprefix'):
            q = ss.name
        else:
            assert 0

        return q

    def init(self, values=[]):
        "initialize a Picklist object with given values."
        if self.pickset is not None:
            raise ValueError("already initialized?")
        self.pickset = set(values)
        return self.pickset

    def load(self, pickfile, column_name):
        "load pickset, return num empty vals, and set of duplicate vals."
        pickset = self.init()

        n_empty_val = 0
        dup_vals = set()
        with open(pickfile, newline='') as csvfile:
            x = csvfile.readline()

            # skip leading comment line in case there's a manifest header
            if x[0] == '#':
                pass
            else:
                csvfile.seek(0)

            r = csv.DictReader(csvfile)

            if column_name not in r.fieldnames:
                raise ValueError(f"column '{column_name}' not in pickfile '{pickfile}'")

            for row in r:
                # pick out values from column
                col = row[column_name]
                if not col:
                    n_empty_val += 1
                    continue

                col = self.preprocess_fn(col)

                # look for duplicate values or empty values
                if col in pickset:
                    dup_vals.add(col)
                else:
                    self.add(col)

        return n_empty_val, dup_vals

    def add(self, value):
        "Add a value to this picklist."
        self.pickset.add(value)

    def __contains__(self, ss):
        "does this signature match anything in the picklist?"
        # pull out the relevant signature attribute
        q = self._get_sig_attribute(ss)

        # mangle into the kinds of values we support here
        q = self.preprocess_fn(q)

        # add to the number of queries performed,
        self.n_queries += 1

        # determine if ok or not.
        if self.pickstyle == PickStyle.INCLUDE:
            if q in self.pickset:
                self.found.add(q)
                return True
        elif self.pickstyle == PickStyle.EXCLUDE:
            if q not in self.pickset:
                self.found.add(q)
                return True
        return False

    def matches_manifest_row(self, row):
        "does the given manifest row match this picklist?"
        if self.coltype == 'md5':
            colkey = 'md5'
        elif self.coltype in ('md5prefix8', 'md5short'):
            colkey = 'md5short'
        elif self.coltype in ('name', 'ident', 'identprefix'):
            colkey = 'name'
        else:
            assert 0

        q = row[colkey]
        q = self.preprocess_fn(q)
        self.n_queries += 1

        if self.pickstyle == PickStyle.INCLUDE:
            if q in self.pickset:
                self.found.add(q)
                return True
        elif self.pickstyle == PickStyle.EXCLUDE:
            if q not in self.pickset:
                self.found.add(q)
                return True
        return False

    def filter(self, it):
        "yield all signatures in the given iterator that are in the picklist"
        for ss in it:
            if self.__contains__(ss):
                yield ss


def passes_all_picklists(ss, picklists):
    "does the signature 'ss' pass all of the picklists?"
    for picklist in picklists:
        if ss not in picklist:
            return False
    return True
