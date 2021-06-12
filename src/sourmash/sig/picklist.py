"Picklist code for extracting subsets of signatures."
import csv


# set up preprocessing functions for column stuff
preprocess = {}

# exact matches
preprocess['name'] = lambda x: x
preprocess['md5'] = lambda x: x

# identifier matches/prefix foo - space delimited identifiers
preprocess['ident.'] = lambda x: x.split(' ')[0].split('.')[0]
preprocess['ident'] = lambda x: x.split(' ')[0]

# match 8 characters
preprocess['md5prefix8'] = lambda x: x[:8]


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
    * 'ident' - exact match to signature's identifier
    * 'ident.' - match to signature's identifier, before '.'

    Identifiers are constructed by using the first space delimited word in
    the signature name.
    """
    def __init__(self, pickfile, column_name, coltype):
        self.pickfile = pickfile
        self.column_name = column_name
        self.coltype = coltype

        if coltype not in ('md5', 'md5prefix8', 'name', 'ident', 'ident.'):
            raise ValueError(f"invalid picklist column type '{coltype}'")

        self.preprocess_fn = preprocess[coltype]
        self.pickset = None

    @classmethod
    def from_picklist_args(cls, argstr):
        "load a picklist from an argument string 'pickfile:column:coltype'"
        picklist = argstr.split(':')
        if len(picklist) != 3:
            raise ValueError(f"invalid picklist argument '{argstr}'")

        assert len(picklist) == 3
        pickfile, column, coltype = picklist

        return cls(pickfile, column, coltype)

    def _get_sig_attribute(self, ss):
        "for a given SourmashSignature, return attribute for this picklist."
        coltype = self.coltype
        if coltype == 'md5':
            q = ss.md5sum()
        elif coltype == 'md5prefix8':
            q = ss.md5sum()
        elif coltype == 'name':
            q = ss.name
        elif coltype == 'ident':
            q = ss.name
        elif coltype == 'ident.':
            q = ss.name

        return q

    def load(self):
        "load pickset, return num empty vals, and set of duplicate vals."
        pickset = set()
        n_empty_val = 0
        dup_vals = set()
        with open(self.pickfile, newline='') as csvfile:
            r = csv.DictReader(csvfile)

            if self.column_name not in r.fieldnames:
                raise ValueError("column '{self.column_name}' not in pickfile '{self.pickfile}'")

            for row in r:
                # pick out values from column
                col = row[self.column_name]
                if not col:
                    n_empty_val += 1
                    continue

                col = self.preprocess_fn(col)

                # look for duplicate values or empty values
                if col in pickset:
                    dup_vals.add(col)
                else:
                    pickset.add(col)

        self.pickset = pickset
        return n_empty_val, dup_vals

    def __contains__(self, ss):
        "does this signature match anything in the picklist?"
        q = self._get_sig_attribute(ss)
        q = self.preprocess_fn(q)

        if q in self.pickset:
            return True
        return False

    def filter(self, it):
        for ss in it:
            if self.__contains__(ss):
                yield ss
