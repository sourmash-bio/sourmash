import json
import gzip
from collections import OrderedDict

taxlist = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus',
           'species']
null_names = set(['[Blank]', 'na', 'null'])


_print_debug = False
def set_debug(state):
    global _print_debug
    _print_debug = True

def debug(*args):
    if _print_debug:
        print(*args)


def build_tree(assignments, initial=None):
    """
    Builds a tree of dictionaries from lists of (rank, name) tuples
    in 'assignments'.  This tree can then be used to find least common
    ancestor agreements/confusion.
    """
    if initial is None:
        tree = {}
    else:
        tree = initial

    for assignment in assignments:
        node = tree

        for rank, name in assignment:
            if name:
                child = node.get((rank, name), {})
                node[(rank, name)] = child

                # shift -> down in tree
                node = child

    return tree


def find_lca(tree):
    """
    Given a tree produced by 'find_tree', find the first node with multiple
    children, OR the only leaf in the tree.  Return ((rank, name), reason),
    where 'reason' is the number of children of the returned node, i.e.e
    0 if it's a leaf and > 1 if it's an internal node.
    """

    node = tree
    cur = ('root', 'root')
    while 1:
        if len(node) == 1:                # descend to only child
            cur = next(iter(node.keys()))
            node = node[cur]
        elif len(node) == 0:              # at leaf; end
            return cur, 0
        else:                             # len(node) > 1 => confusion!!
            return cur, len(node)


def build_reverse_tree(assignments, initial=None):
    """
    Builds a child -> parent dictionary (a reverse DAG) from lists of
    (rank, name) tuples in 'assignments'.
    """
    if initial is None:
        parents = {}
    else:
        parents = initial

    for assignment in assignments:
        last_node = ('root', 'root')
        for rank, name in assignment:
            if name:
                parents[(rank, name)] = last_node
                last_node = (rank, name)

    return parents


class LCA_Database(object):
    """
    Wrapper class for taxonomic database.

    obj.lineage_dict: key 'lineage_id' => lineage tuple [(name, rank), ...]
    obj.hashval_to_lineage_id: key 'hashval' => 'lineage_id'
    obj.ksize: k-mer size
    obj.scaled: scaled value
    obj.signatures_to_lineage: key 'md5sum' => 'lineage_id'
    """
    def __init__(self):
        self.lineage_dict = None
        self.hashval_to_lineage_id = None
        self.ksize = None
        self.scaled = None
        self.signatures_to_lineage = None

    def load(self, db_name):
        "Load from a JSON file."
        xopen = open
        if db_name.endswith('.gz'):
            xopen = gzip.open

        with xopen(db_name, 'rt') as fp:
            load_d = json.load(fp)
            version = load_d['version']
            assert version == '1.0'

            type = load_d['type']
            assert type == 'sourmash_lca'

            ksize = load_d['ksize']
            scaled = load_d['scaled']

            lineage_dict_2 = load_d['lineages']
            lineage_dict = {}
            for k, v in lineage_dict_2.items():
                vv = []
                for rank in taxlist:
                    name = v.get(rank, '')
                    vv.append((rank, name))

                lineage_dict[int(k)] = tuple(vv)

            hashval_to_lineage_id_2 = load_d['hashval_assignments']
            hashval_to_lineage_id = {}
            for k, v in hashval_to_lineage_id_2.items():
                hashval_to_lineage_id[int(k)] = v

            signatures_to_lineage = load_d['signatures_to_lineage']

        self.lineage_dict = lineage_dict
        self.hashval_to_lineage_id = hashval_to_lineage_id
        self.ksize = ksize
        self.scaled = scaled
        self.signatures_to_lineage = signatures_to_lineage

    def save(self, db_name):
        "Save to a JSON file."
        xopen = open
        if db_name.endswith('.gz'):
            xopen = gzip.open

        with xopen(db_name, 'wt') as fp:
            # use an OrderedDict to preserve output order
            save_d = OrderedDict()
            save_d['version'] = '1.0'
            save_d['type'] = 'sourmash_lca'
            save_d['license'] = 'CC0'
            save_d['ksize'] = self.ksize
            save_d['scaled'] = self.scaled

            # convert lineage internals from tuples to dictionaries
            save_d['lineages'] = OrderedDict([ (k, OrderedDict(v)) \
                                     for k, v in self.lineage_dict.items() ])

            # convert values from sets to lists, so that JSON knows how to save
            save_d['hashval_assignments'] = \
               dict((k, list(v)) for (k, v) in self.hashval_to_lineage_id.items())
            save_d['signatures_to_lineage'] = self.signatures_to_lineage
            json.dump(save_d, fp)
