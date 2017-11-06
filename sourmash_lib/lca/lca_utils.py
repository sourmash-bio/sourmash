import json

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


def test_build_tree():
    tree = build_tree([[('rank1', 'name1'), ('rank2', 'name2')]])
    assert tree == { ('rank1', 'name1'): { ('rank2', 'name2') : {}} }


def test_build_tree_2():
    tree = build_tree([[('rank1', 'name1'), ('rank2', 'name2a')],
                       [('rank1', 'name1'), ('rank2', 'name2b')],
                      ])

    assert tree == { ('rank1', 'name1'): { ('rank2', 'name2a') : {},
                                           ('rank2', 'name2b') : {}} }


def test_build_tree_3():                  # empty 'rank2' name
    tree = build_tree([[('rank1', 'name1'), ('rank2', '')]])
    assert tree == { ('rank1', 'name1'): {} }


def test_build_tree_4():
    tree = build_tree([[('rank1', 'name1'), ('rank2', 'name2a')],
                      ])

    tree = build_tree([[('rank1', 'name1'), ('rank2', 'name2b')],
                      ], tree)

    assert tree == { ('rank1', 'name1'): { ('rank2', 'name2a') : {},
                                           ('rank2', 'name2b') : {}} }


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


def test_find_lca():
    tree = build_tree([[('rank1', 'name1'), ('rank2', 'name2')]])
    lca = find_lca(tree)

    assert lca == (('rank2', 'name2'), 0)


def test_find_lca_2():
    tree = build_tree([[('rank1', 'name1'), ('rank2', 'name2a')],
                       [('rank1', 'name1'), ('rank2', 'name2b')],
                      ])
    lca = find_lca(tree)

    assert lca == (('rank1', 'name1'), 2)


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


def test_build_reverse_tree():
    parents = build_reverse_tree([[('rank1', 'name1'), ('rank2', 'name2')]])

    print(parents)
    assert parents == { ('rank2', 'name2'): ('rank1', 'name1'),
                        ('rank1', 'name1'): ('root', 'root') }


def test_build_reverse_tree_2():
    parents = build_reverse_tree([[('rank1', 'name1'), ('rank2', 'name2a')],
                                 [('rank1', 'name1'), ('rank2', 'name2b')],
                                 ])

    assert parents == { ('rank2', 'name2a'): ('rank1', 'name1'),
                        ('rank2', 'name2b'): ('rank1', 'name1'),
                        ('rank1', 'name1'): ('root', 'root') }


def test_build_reverse_tree_3():
    parents = build_reverse_tree([[('rank1', 'name1'), ('rank2', 'name2a')],
                                 ])
    parents = build_reverse_tree([[('rank1', 'name1'), ('rank2', 'name2b')],
                                 ], parents)

    assert parents == { ('rank2', 'name2a'): ('rank1', 'name1'),
                        ('rank2', 'name2b'): ('rank1', 'name1'),
                        ('rank1', 'name1'): ('root', 'root') }


def test_build_reverse_tree_4():          # empty 'rank2' name
    parents = build_reverse_tree([[('rank1', 'name1'), ('rank2', '')]])

    print(parents)
    assert parents == { ('rank1', 'name1'): ('root', 'root') }


class LCA_Database(object):
    def __init__(self):
        self.lineage_dict = None
        self.hashval_to_lineage_id = None
        self.ksize = None
        self.scaled = None
        self.signatures_to_lineage = None

    def load(self, db_name):
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


