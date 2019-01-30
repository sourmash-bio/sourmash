"""
Utility functions for lowest-common-ancestor analysis tools.
"""
from __future__ import print_function, division
import sys
import json
import gzip
from os.path import exists
from collections import OrderedDict, namedtuple, defaultdict, Counter

__all__ = ['taxlist', 'zip_lineage', 'build_tree', 'find_lca',
           'load_single_database', 'load_databases', 'gather_assignments',
           'count_lca_for_assignments', 'LineagePair']

try:                                      # py2/py3 compat
    from itertools import zip_longest
except ImportError:
    from itertools import izip_longest as zip_longest
import pprint

from .._minhash import get_max_hash_for_scaled
from ..logging import notify, error, debug

# type to store an element in a taxonomic lineage
LineagePair = namedtuple('LineagePair', ['rank', 'name'])


def check_files_exist(*files):
    ret = True
    not_found = []
    for f in files:
        if not exists(f):
            not_found.append(f)
            ret = False

    if len(not_found):
        error('Error! Could not find the following files.'
              ' Make sure the file paths are specified correctly.\n{}'.format('\n'.join(not_found)))

    return ret


# ordered list of taxonomic ranks
def taxlist(include_strain=True):
    """
    Provide an ordered list of taxonomic ranks.
    """
    for k in ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus',
              'species']:
        yield k
    if include_strain:
        yield 'strain'


# produce an ordered list of tax names from lineage
def zip_lineage(lineage, include_strain=True, truncate_empty=False):
    """
    Given an iterable of LineagePair objects, return list of lineage names.

    This utility function handles species/strain and empty lineage entries
    gracefully.

    >>> x = [ LineagePair('superkingdom', 'a'), LineagePair('phylum', 'b') ]
    >>> zip_lineage(x)
    ['a', 'b', '', '', '', '', '', '']

    >>> x = [ LineagePair('superkingdom', 'a'), LineagePair(None, ''), LineagePair('class', 'c') ]
    >>> zip_lineage(x)
    ['a', '', 'c', '', '', '', '', '']
    """

    row = []
    empty = LineagePair(None, '')
    for taxrank, lineage_tup in zip_longest(taxlist(include_strain=include_strain), lineage, fillvalue=empty):
        if lineage_tup == empty:
            if truncate_empty:
                break
        else:
            # validate non-empty tax, e.g. superkingdom/phylum/class in order.
            if lineage_tup.rank != taxrank:
                raise ValueError('incomplete lineage at {} - is {} instead'.format(taxrank, lineage_tup.rank))

        row.append(lineage_tup.name)
    return row


# filter function toreplace blank/na/null with 'unassigned'
filter_null = lambda x: 'unassigned' if x.strip() in \
  ('[Blank]', 'na', 'null', '') else x
null_names = set(['[Blank]', 'na', 'null'])


def build_tree(assignments, initial=None):
    """
    Builds a tree of dictionaries from lists of LineagePair objects
    in 'assignments'.  This tree can then be used to find lowest common
    ancestor agreements/confusion.
    """
    if initial is None:
        tree = {}
    else:
        tree = initial

    if not assignments:
        raise ValueError("empty assignment passed to build_tree")

    for assignment in assignments:
        node = tree

        for lineage_tup in assignment:
            if lineage_tup.name:
                child = node.get(lineage_tup, {})
                node[lineage_tup] = child

                # shift -> down in tree
                node = child

    return tree


def find_lca(tree):
    """
    Given a tree produced by 'find_tree', find the first node with multiple
    children, OR the only leaf in the tree.  Return (lineage_tup, reason),
    where 'reason' is the number of children of the returned node, i.e.e
    0 if it's a leaf and > 1 if it's an internal node.
    """

    node = tree
    lineage = []
    while 1:
        if len(node) == 1:                # descend to only child; track path
            lineage_tup = next(iter(node.keys()))
            lineage.append(lineage_tup)
            node = node[lineage_tup]
        elif len(node) == 0:              # at leaf; end
            return tuple(lineage), 0
        else:                             # len(node) > 1 => confusion!!
            return tuple(lineage), len(node)


class LCA_Database(object):
    """
    Wrapper class for taxonomic database.

    obj.ident_to_idx: key 'identifier' to 'idx'
    obj.idx_to_lid: key 'idx' to 'lid'
    obj.lid_to_lineage: key 'lid' to tuple of LineagePair objects
    obj.hashval_to_idx: key 'hashval' => set('idx')
    obj.lineage_to_lid: key (tuple of LineagePair objects) to 'lid'
    """
    def __init__(self):
        self.ksize = None
        self.scaled = None
        
        self.ident_to_idx = None
        self.idx_to_lid = None
        self.lineage_to_lid = None
        self.lid_to_lineage = None
        self.hashval_to_idx = None
    
        self.filename = None

    def __repr__(self):
        return "LCA_Database('{}')".format(self.filename)

    def load(self, db_name):
        "Load from a JSON file."
        xopen = open
        if db_name.endswith('.gz'):
            xopen = gzip.open

        with xopen(db_name, 'rt') as fp:
            load_d = {}
            try:
                load_d = json.load(fp)
            except json.decoder.JSONDecodeError:
                pass

            if not load_d:
                raise ValueError("cannot parse database file '{}' as JSON; invalid format.")

            version = None
            db_type = None
            try:
                version = load_d.get('version')
                db_type = load_d.get('type')
            except AttributeError:
                pass

            if db_type != 'sourmash_lca':
                raise ValueError("database file '{}' is not an LCA db.".format(db_name))

            if version != '2.0' or 'lid_to_lineage' not in load_d:
                raise ValueError("Error! This is an old-style LCA DB. You'll need to build or download a newer one.")

            ksize = int(load_d['ksize'])
            scaled = int(load_d['scaled'])
            self.ksize = ksize
            self.scaled = scaled

            # convert lineage_dict to proper lineages (tuples of LineagePairs)
            lid_to_lineage_2 = load_d['lid_to_lineage']
            lid_to_lineage = {}
            for k, v in lid_to_lineage_2.items():
                v = dict(v)
                vv = []
                for rank in taxlist():
                    name = v.get(rank, '')
                    vv.append(LineagePair(rank, name))

                lid_to_lineage[int(k)] = tuple(vv)
            self.lid_to_lineage = lid_to_lineage

            # convert hashval -> lineage index keys to integers (looks like
            # JSON doesn't have a 64 bit type so stores them as strings)
            hashval_to_idx_2 = load_d['hashval_to_idx']
            hashval_to_idx = {}

            for k, v in hashval_to_idx_2.items():
                hashval_to_idx[int(k)] = v
            self.hashval_to_idx = hashval_to_idx

            self.ident_to_name = load_d['ident_to_name']
            self.ident_to_idx = load_d['ident_to_idx']

            self.idx_to_lid = {}
            for k, v in load_d['idx_to_lid'].items():
                self.idx_to_lid[int(k)] = v

        self.filename = db_name

    def save(self, db_name):
        "Save to a JSON file."
        xopen = open
        if db_name.endswith('.gz'):
            xopen = gzip.open

        with xopen(db_name, 'wt') as fp:
            # use an OrderedDict to preserve output order
            save_d = OrderedDict()
            save_d['version'] = '2.0'
            save_d['type'] = 'sourmash_lca'
            save_d['license'] = 'CC0'
            save_d['ksize'] = self.ksize
            save_d['scaled'] = self.scaled

            # convert lineage internals from tuples to dictionaries
            d = OrderedDict()
            for k, v in self.lid_to_lineage.items():
                d[k] = dict([ (vv.rank, vv.name) for vv in v ])
            save_d['lid_to_lineage'] = d

            # convert values from sets to lists, so that JSON knows how to save
            save_d['hashval_to_idx'] = \
               dict((k, list(v)) for (k, v) in self.hashval_to_idx.items())

            save_d['ident_to_name'] = self.ident_to_name
            save_d['ident_to_idx'] = self.ident_to_idx
            save_d['idx_to_lid'] = self.idx_to_lid
            save_d['lid_to_lineage'] = self.lid_to_lineage
            
            json.dump(save_d, fp)

    def downsample_scaled(self, scaled):
        """
        Downsample to the provided scaled value, i.e. eliminate all hashes
        that don't fall in the required range.
        """
        if scaled == self.scaled:
            return
        elif scaled < self.scaled:
            raise ValueError("cannot decrease scaled from {} to {}".format(self.scaled, scaled))

        max_hash = get_max_hash_for_scaled(scaled)
        new_hashvals = {}
        for k, v in self.hashval_to_idx.items():
            if k < max_hash:
                new_hashvals[k] = v
        self.hashval_to_idx = new_hashvals
        self.scaled = scaled

    def get_lineage_assignments(self, hashval):
        """
        Get a list of lineages for this hashval.
        """
        x = []

        idx_list = self.hashval_to_idx.get(hashval, [])
        for idx in idx_list:
            lid = self.idx_to_lid.get(idx, None)
            if lid is not None:
                lineage = self.lid_to_lineage[lid]
                x.append(lineage)

        return x

    def find(self, minhash, threshold, containment=False, ignore_scaled=False):
        """
        Do a Jaccard similarity or containment search.
        """
        # make sure we're looking at the same scaled value as database
        if self.scaled > minhash.scaled:
            minhash = minhash.downsample_scaled(self.scaled)
        elif self.scaled < minhash.scaled and not ignore_scaled:
            raise ValueError("lca db scaled is {} vs query {}; must downsample".format(self.scaled, minhash.scaled))

        if not hasattr(self, 'signatures'):
            debug('creating signatures for LCA DB...')
            sigd = defaultdict(minhash.copy_and_clear)

            for (k, v) in self.hashval_to_idx.items():
                for vv in v:
                    sigd[vv].add_hash(k)

            self.signatures = sigd

        debug('=> {} signatures!', len(self.signatures))

        # build idx_to_ident from ident_to_idx
        if not hasattr(self, 'idx_to_ident'):
            idx_to_ident = {}
            for k, v in self.ident_to_idx.items():
                idx_to_ident[v] = k

            self.idx_to_ident = idx_to_ident

        query_mins = set(minhash.get_mins())

        # collect matching hashes:
        c = Counter()
        for hashval in query_mins:
            idx_list = self.hashval_to_idx.get(hashval, [])
            for idx in idx_list:
                c[idx] += 1

        debug('number of matching signatures for hashes: {}', len(c))

        for idx, count in c.items():
            ident = self.idx_to_ident[idx]
            name = self.ident_to_name[ident]
            debug('looking at {} ({})', ident, name)

            match_mh = self.signatures[idx]
            match_size = len(match_mh)

            debug('count: {}; query_mins: {}; match size: {}',
                  count, len(query_mins), match_size)

            if containment:
                score = count / len(query_mins)
            else:
                score = count / (len(query_mins) + match_size - count)

            debug('score: {} (containment? {})', score, containment)

            if score >= threshold:
                # reconstruct signature... ugh.
                from .. import SourmashSignature
                match_sig = SourmashSignature(match_mh, name=name)

                yield score, match_sig, match_sig.md5sum(), self.filename, name


def load_single_database(filename, verbose=False):
    "Load a single LCA database; return (db, ksize, scaled)"
    dblist, ksize, scaled = load_databases([filename], verbose=verbose)
    return dblist[0], ksize, scaled


def load_databases(filenames, scaled=None, verbose=True):
    "Load multiple LCA databases; return (dblist, ksize, scaled)"
    ksize_vals = set()
    scaled_vals = set()
    dblist = []

    # load all the databases
    for db_name in filenames:
        if verbose:
            notify(u'\r\033[K', end=u'', file=sys.stderr)
            notify('... loading database {}'.format(db_name), end='\r',
                  file=sys.stderr)

        lca_db = LCA_Database()
        lca_db.load(db_name)

        ksize_vals.add(lca_db.ksize)
        if len(ksize_vals) > 1:
            raise Exception('multiple ksizes, quitting')

        if scaled and scaled > lca_db.scaled:
            lca_db.downsample_scaled(scaled)
        scaled_vals.add(lca_db.scaled)

        dblist.append(lca_db)

    ksize = ksize_vals.pop()
    scaled = scaled_vals.pop()

    if verbose:
        notify(u'\r\033[K', end=u'')
        notify('loaded {} LCA databases. ksize={}, scaled={}', len(dblist),
               ksize, scaled)

    return dblist, ksize, scaled


def gather_assignments(hashvals, dblist):
    """
    Gather assignments from across all the databases for all the hashvals.
    """
    assignments = defaultdict(set)
    for hashval in hashvals:
        for lca_db in dblist:
            lineages = lca_db.get_lineage_assignments(hashval)
            if lineages:
                assignments[hashval].update(lineages)

    return assignments


def count_lca_for_assignments(assignments):
    """
    For each hashval, count the LCA across its assignments.
    """
    counts = Counter()
    for hashval in assignments:

        # for each list of tuple_info [(rank, name), ...] build
        # a tree that lets us discover lowest-common-ancestor.
        lineages = assignments[hashval]
        tree = build_tree(lineages)

        # now find either a leaf or the first node with multiple
        # children; that's our lowest-common-ancestor node.
        lca, reason = find_lca(tree)
        counts[lca] += 1

    return counts
