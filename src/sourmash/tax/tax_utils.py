"""
Utility functions for taxonomy analysis tools.
"""
import os
import csv
from collections import namedtuple, defaultdict
from collections import abc
from itertools import zip_longest
from typing import NamedTuple
from dataclasses import dataclass, field, replace, asdict
import gzip

from sourmash import sqlite_utils, sourmash_args
from sourmash.exceptions import IndexNotSupported
from sourmash.distance_utils import containment_to_distance

import sqlite3


__all__ = ['get_ident', 'ascending_taxlist', 'collect_gather_csvs',
           'load_gather_results', 'check_and_load_gather_csvs',
           'find_match_lineage', 'summarize_gather_at',
           'find_missing_identities', 'make_krona_header',
           'aggregate_by_lineage_at_rank', 'format_for_krona',
           'write_krona', 'write_summary', 'write_classifications',
           'combine_sumgather_csvs_by_lineage', 'write_lineage_sample_frac',
           'MultiLineageDB', 'RankLineageInfo']

from sourmash.logging import notify
from sourmash.sourmash_args import load_pathlist_from_file

# CTB: these could probably usefully be converted into dataclasses.
QInfo = namedtuple("QInfo", "query_md5, query_filename, query_bp, query_hashes, total_weighted_hashes")
SumGathInf = namedtuple("SumGathInf", "query_name, rank, fraction, lineage, query_md5, query_filename, f_weighted_at_rank, bp_match_at_rank, query_ani_at_rank, total_weighted_hashes")
ClassInf = namedtuple("ClassInf", "query_name, status, rank, fraction, lineage, query_md5, query_filename, f_weighted_at_rank, bp_match_at_rank, query_ani_at_rank")

# Essential Gather column names that must be in gather_csv to allow `tax` summarization
EssentialGatherColnames = ('query_name', 'name', 'f_unique_weighted', 'f_unique_to_query', 'unique_intersect_bp', 'remaining_bp', 'query_md5', 'query_filename')

# import lca utils as needed for now
from sourmash.lca import lca_utils
from sourmash.lca.lca_utils import (taxlist, display_lineage, pop_to_rank)

class LineagePair(NamedTuple):
    rank: str
    name: str = None
    taxid: int = None

@dataclass(frozen=True, order=True)
class BaseLineageInfo:
    """BaseLineageInfo Class, defining a set of functions that can be used to handle
       hierarchical set of LineagePairs. """
    # need to set compare=False for any mutable type to keep this class hashable
    ranks: tuple() # require ranks
    lineage: tuple = field(default=()) #tuple of LineagePairs
    lineage_str: str = field(default=None, compare=False) # ';'- or ','-separated str of lineage names
    lineage_dict: dict = field(default=None, compare=False) # dict of rank: name

    def __post_init__(self):
        "Initialize according to passed values"
        # ranks must be tuple for hashability
        if isinstance(self.ranks, list):
            object.__setattr__(self, "ranks", tuple(self.ranks))
        if self.lineage:
            self._init_from_lineage_tuples()
        elif self.lineage_str is not None:
            self._init_from_lineage_str()
        elif self.lineage_dict is not None:
            self._init_from_lineage_dict()
        else:
            self._init_empty()

    def __eq__(self, other):
        if other == (): # just handy: if comparing to a null tuple, don't try to find it's lineage before returning False
            return False
        return all([self.ranks == other.ranks and self.lineage==other.lineage])

    @property
    def taxlist(self):
        return self.ranks

    @property
    def ascending_taxlist(self):
        return self.ranks[::-1]

    @property
    def lowest_rank(self):
        if not self.filled_ranks:
            return None
        return self.filled_ranks[-1]

    def rank_index(self, rank):
        return self.ranks.index(rank)

    @property
    def filled_lineage(self):
        """Return lineage down to lowest non-empty rank. Preserves missing ranks above."""
        # Would we prefer this to be the default returned by lineage??
        if not self.filled_ranks:
            return ()
        lowest_filled_rank_idx = self.rank_index(self.filled_ranks[-1])
        return self.lineage[:lowest_filled_rank_idx+1]

    @property
    def lowest_lineage_name(self):
        """Return the name of the lowest filled lineage"""
        if not self.filled_ranks:
            return ""
        return self.filled_lineage[-1].name

    @property
    def lowest_lineage_taxid(self):
        """Return the taxid of the lowest filled lineage"""
        if not self.filled_ranks:
            return ""
        return self.filled_lineage[-1].taxid

    def _init_empty(self):
        'initialize empty genome lineage'
        new_lineage = []
        for rank in self.ranks:
            new_lineage.append(LineagePair(rank=rank))
        # set lineage and filled_ranks (because frozen, need to do it this way)
        object.__setattr__(self, "lineage", new_lineage)
        object.__setattr__(self, "filled_ranks", ())

    def _init_from_lineage_tuples(self):
        'initialize from tuple/list of LineagePairs, allowing empty ranks and reordering if necessary'
        new_lineage = []
        # check this is a list or tuple of lineage tuples:
        for rank in self.ranks:
            new_lineage.append(LineagePair(rank=rank))
        for lin_tup in self.lineage:
            # now add input tuples in correct spots. This corrects for order and allows empty values.
            if not isinstance(lin_tup, (LineagePair, lca_utils.LineagePair)):
                raise ValueError(f"{lin_tup} is not LineagePair.")
                # find index for this rank
            if lin_tup.rank: # skip this tuple if rank is None or "" (empty lineage tuple. is this needed?)
                try:
                    rank_idx = self.rank_index(lin_tup.rank)
                except ValueError as e:
                    raise ValueError(f"Rank '{lin_tup.rank}' not present in {', '.join(self.ranks)}") from e
                # make sure we're adding tax_utils.LineagePairs, not lca_utils.LineagePairs for consistency
                if isinstance(lin_tup, lca_utils.LineagePair):
                    new_lineage[rank_idx] =  LineagePair(rank=lin_tup.rank, name=lin_tup.name)
                else:
                    new_lineage[rank_idx] =  lin_tup
        # build list of filled ranks
        filled_ranks = [a.rank for a in new_lineage if a.name]
        # set lineage and filled_ranks
        object.__setattr__(self, "lineage", tuple(new_lineage))
        object.__setattr__(self, "filled_ranks", filled_ranks)

    def _init_from_lineage_dict(self):
        'initialize from lineage dict, e.g. from gather csv, allowing empty ranks and reordering if necessary'
        if not isinstance(self.lineage_dict, (dict)):
            raise ValueError(f"{self.lineage_dict} is not dictionary")
        # first, initialize_empty
        new_lineage = []
        # build empty lineage
        for rank in self.ranks:
            new_lineage.append(LineagePair(rank=rank))
        # now add input information in correct spots. This corrects for order and allows empty values.
        for rank, info in self.lineage_dict.items():
            try:
                rank_idx = self.rank_index(rank)
            except ValueError as e:
                raise ValueError(f"Rank '{rank}' not present in {', '.join(self.ranks)}") from e

            name, taxid = None, None
            if isinstance(info, dict):
                if 'name' in info.keys():
                    name = info['name']
                if 'taxid' in info.keys():
                    taxid = info['taxid']
            elif isinstance(info, str):
                name = info
            new_lineage[rank_idx] =  LineagePair(rank=rank, name=name, taxid=taxid)
        # build list of filled ranks
        filled_ranks = [a.rank for a in new_lineage if a.name]
        # set lineage and filled_ranks
        object.__setattr__(self, "lineage", tuple(new_lineage))
        object.__setattr__(self, "filled_ranks", filled_ranks)

    def _init_from_lineage_str(self):
        """
        Turn a ; or ,-separated set of lineages into a list of LineagePair objs.
        """
        new_lineage = self.lineage_str.split(';')
        if len(new_lineage) == 1:
            new_lineage = self.lineage_str.split(',')
        new_lineage = [ LineagePair(rank=rank, name=n) for (rank, n) in zip_longest(self.ranks, new_lineage) ]
        # build list of filled ranks
        filled_ranks = [a.rank for a in new_lineage if a.name]
        object.__setattr__(self, "lineage", tuple(new_lineage))
        object.__setattr__(self, "filled_ranks", filled_ranks)

    def zip_lineage(self, truncate_empty=False):
        """
        Return lineage names as a list
        """
        if truncate_empty:
            zipped = [a.name for a in self.filled_lineage]
        else:
            zipped = [a.name for a in self.lineage]
        # replace None with empty string ("")
        if None in zipped:
            zipped = ['' if x is None else x for x in zipped]

        return zipped

    def zip_taxid(self, truncate_empty=False):
        """
        Return taxids as a list
        """
        if truncate_empty:
            zipped = [a.taxid for a in self.filled_lineage]
        else:
            zipped = [a.taxid for a in self.lineage]
        # replace None with empty string (""); cast taxids to str
        zipped = ['' if x is None else str(x) for x in zipped]

        return zipped

    def display_lineage(self, truncate_empty=True, null_as_unclassified=False):
        "Return lineage names as ';'-separated list"
        lin = ";".join(self.zip_lineage(truncate_empty=truncate_empty))
        if null_as_unclassified:
            if lin == "":
                return "unclassified"
        else:
            return lin

    def display_taxid(self, truncate_empty=True):
        "Return lineage taxids as ';'-separated list"
        return ";".join(self.zip_taxid(truncate_empty=truncate_empty))

    def check_rank_availability(self, rank):
        if rank in self.ranks: # rank is available
            return True
        raise ValueError(f"Desired Rank '{rank}' not available for this lineage.")
 
    def rank_is_filled(self, rank, other=None):
        self.check_rank_availability(rank)
        if other is not None:
            if rank in self.filled_ranks and rank in other.filled_ranks:
                return True
        elif rank in self.filled_ranks:
            return True
        return False

    def is_lineage_match(self, other, rank):
        """
        check to see if two lineages are a match down to given rank.
        """
        if not other.ranks == self.ranks: # check same ranks
            raise ValueError("Cannot compare lineages from taxonomies with different ranks.")
        # always return false if rank is not filled in either of the two lineages
        if self.rank_is_filled(rank, other=other):
            rank_idx = self.rank_index(rank)
            a_lin = self.lineage[:rank_idx+1]
            b_lin = other.lineage[:rank_idx+1]
            if a_lin == b_lin:
                return 1
        return 0

    def pop_to_rank(self, rank):
        "Return new LineageInfo with ranks only filled to desired rank"
        # are we already above rank?
        if not self.rank_is_filled(rank):
            return replace(self)
        # if not, make filled_lineage at this rank + use to generate new LineageInfo
        new_lineage = self.lineage_at_rank(rank)
        new = replace(self, lineage = new_lineage)
        # replace doesn't run the __post_init__ properly. reinitialize.
        new._init_from_lineage_tuples()
        return new

    def lineage_at_rank(self, rank):
        # non-descructive pop_to_rank. Returns tuple of LineagePairs
        "Returns tuple of LineagePairs at given rank."
        # are we already above rank?
        if not self.rank_is_filled(rank):
            return self.filled_lineage
        # if not, return lineage tuples down to desired rank
        rank_idx = self.rank_index(rank)
        return self.filled_lineage[:rank_idx+1]


@dataclass(frozen=True, order=True)
class RankLineageInfo(BaseLineageInfo):
    """Class for storing multi-rank lineage information"""
    ranks: tuple = field(default_factory=lambda: ('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain'))

    def __post_init__(self):
        "Initialize according to passed values"
        # ranks must be tuple for hashability
        if isinstance(self.ranks, list):
            object.__setattr__(self, "ranks", tuple(self.ranks))
        if self.lineage:
            self._init_from_lineage_tuples()
        elif self.lineage_str is not None:
            self._init_from_lineage_str()
        elif self.lineage_dict is not None:
            self._init_from_lineage_dict()
        elif self.ranks:
            self._init_empty()


def get_ident(ident, *,
              keep_full_identifiers=False, keep_identifier_versions=False):
    # split identifiers = split on whitespace
    # keep identifiers = don't split .[12] from assembly accessions
    "Hack and slash identifiers."
    if not keep_full_identifiers:
        ident = ident.split(' ')[0]
        if not keep_identifier_versions:
            ident = ident.split('.')[0]
    return ident


def ascending_taxlist(include_strain=True):
    """
    Provide an ordered list of taxonomic ranks: strain --> superkingdom
    """
    ascending_taxlist = ['species', 'genus', 'family', 'order',
                         'class', 'phylum', 'superkingdom']
    if include_strain:
        ascending_taxlist = ['strain'] + ascending_taxlist
    for k in ascending_taxlist:
        yield k


def collect_gather_csvs(cmdline_gather_input, *, from_file=None):
    """
    collect gather files from cmdline; --from-file input
    """
    gather_csvs = []
    # ignore command line duplicates
    for gf in cmdline_gather_input:
        if gf not in gather_csvs:
            gather_csvs.append(gf)
        else:
            notify(f'ignoring duplicated reference to file: {gf}')
    # ignore pathlist duplicates
    if from_file:
        more_files = load_pathlist_from_file(from_file)
        for gf in more_files:
            if gf not in gather_csvs:
                gather_csvs.append(gf)
            else:
               notify(f'ignoring duplicated reference to file: {gf}')
    return gather_csvs


def load_gather_results(gather_csv, *, essential_colnames=EssentialGatherColnames,
                        seen_queries=None, force=False):
    "Load a single gather csv"
    if not seen_queries:
        seen_queries=set()
    header = []
    gather_results = []
    gather_queries = set()
    with sourmash_args.FileInputCSV(gather_csv) as r:
        header = r.fieldnames
        # check for empty file
        if not header:
            raise ValueError(f"Cannot read gather results from '{gather_csv}'. Is file empty?")

        # check for critical column names used by summarize_gather_at
        if not set(essential_colnames).issubset(header):
            raise ValueError(f"Not all required gather columns are present in '{gather_csv}'.")

        for n, row in enumerate(r):
            query_name = row['query_name']
            # check if we've seen this query already in a different gather CSV
            if query_name in seen_queries:
                # do not allow loading of same query from a second CSV.
                raise ValueError(f"Gather query {query_name} was found in more than one CSV. Cannot load from '{gather_csv}'.")
            else:
                gather_results.append(row)
            # add query name to the gather_queries from this CSV
            if query_name not in gather_queries:
                gather_queries.add(query_name)

    if not gather_results:
        raise ValueError(f'No gather results loaded from {gather_csv}.')
    else:
        notify(f"loaded {len(gather_results)} gather results from '{gather_csv}'.")
    return gather_results, header, gather_queries


def check_and_load_gather_csvs(gather_csvs, tax_assign, *, fail_on_missing_taxonomy=False, force=False):
    '''
    Load gather csvs, checking for empties and ids missing from taxonomic assignments.
    '''
    if not isinstance(gather_csvs, list):
        gather_csvs = [gather_csvs]
    gather_results = []
    total_missed = 0
    all_ident_missed = set()
    seen_queries = set()
    header = []
    n_ignored = 0
    for n, gather_csv in enumerate(gather_csvs):
        these_results = []
        try:
            these_results, header, seen_queries = load_gather_results(gather_csv, seen_queries=seen_queries, force=force)
        except ValueError as exc:
            if force:
                if "found in more than one CSV" in str(exc):
                    notify('Cannot force past duplicated gather query. Exiting.')
                    raise
                notify(str(exc))
                notify('--force is set. Attempting to continue to next set of gather results.')
                n_ignored+=1
                continue
            else:
                notify('Exiting.')
                raise

        # check for match identites in these gather_results not found in lineage spreadsheets
        ident_missed = find_missing_identities(these_results, tax_assign)
        if ident_missed:
            notify(f'The following are missing from the taxonomy information: {",".join(ident_missed)}')
            if fail_on_missing_taxonomy:
                raise ValueError('Failing on missing taxonomy, as requested via --fail-on-missing-taxonomy.')

            total_missed += len(ident_missed)
            all_ident_missed.update(ident_missed)
        # add these results to gather_results
        gather_results += these_results

    num_gather_csvs_loaded = n+1 - n_ignored
    notify(f'loaded {len(gather_results)} results total from {str(num_gather_csvs_loaded)} gather CSVs')

    return gather_results, all_ident_missed, total_missed, header


def find_match_lineage(match_ident, tax_assign, *, skip_idents = [],
                       keep_full_identifiers=False,
                       keep_identifier_versions=False):
    lineage=""
    match_ident = get_ident(match_ident, keep_full_identifiers=keep_full_identifiers, keep_identifier_versions=keep_identifier_versions)
    # if identity not in lineage database, and not --fail-on-missing-taxonomy, skip summarizing this match
    if match_ident in skip_idents:
        return lineage
    try:
        lineage = tax_assign[match_ident]
    except KeyError:
        raise ValueError(f"ident {match_ident} is not in the taxonomy database.")
    return lineage


def summarize_gather_at(rank, tax_assign, gather_results, *, skip_idents = [],
                        keep_full_identifiers=False,
                        keep_identifier_versions=False, best_only=False,
                        seen_perfect=set(),
                        estimate_query_ani=False):
    """
    Summarize gather results at specified taxonomic rank
    """
    # init dictionaries
    sum_uniq_weighted = defaultdict(lambda: defaultdict(float))
    # store together w/ ^ instead?
    sum_uniq_to_query = defaultdict(lambda: defaultdict(float))
    sum_uniq_bp = defaultdict(lambda: defaultdict(float))
    query_info = {}

    set_ksize = False
    ksize, scaled, query_nhashes = None, 0, None

    for row in gather_results:
        # get essential gather info
        if not set_ksize and "ksize" in row.keys():
            set_ksize = True
            ksize = int(row['ksize'])
            scaled = int(row['scaled'])
        
        query_name = row['query_name']
        f_unique_to_query = float(row['f_unique_to_query'])
        f_uniq_weighted = float(row['f_unique_weighted'])
        unique_intersect_bp = int(row['unique_intersect_bp'])
        total_weighted_hashes = int(row.get('total_weighted_hashes', 0))
        query_md5 = row['query_md5']
        query_filename = row['query_filename']
        # get query_bp
        if query_name not in query_info.keys(): #REMOVING THIS AFFECTS GATHER RESULTS!!! BUT query bp should always be same for same query? bug?
            if "query_nhashes" in row.keys():
                query_nhashes = int(row["query_nhashes"])
            if "query_bp" in row.keys():
                query_bp = int(row["query_bp"])
            else:
                query_bp = unique_intersect_bp + int(row['remaining_bp'])
        
        # store query info
        query_info[query_name] = QInfo(query_md5=query_md5, query_filename=query_filename, query_bp=query_bp, query_hashes=query_nhashes, total_weighted_hashes=total_weighted_hashes)
        
        if estimate_query_ani and (not ksize or not scaled):
            if not set_ksize:
                estimate_query_ani=False
                notify("WARNING: Please run gather with sourmash >= 4.4 to estimate query ANI at rank. Continuing without ANI...")
        
        match_ident = row['name']

        # 100% match? are we looking at something in the database?
        if f_unique_to_query >= 1.0 and query_name not in seen_perfect: # only want to notify once, not for each rank
            ident = get_ident(match_ident,
                              keep_full_identifiers=keep_full_identifiers,
                              keep_identifier_versions=keep_identifier_versions)
            seen_perfect.add(query_name)
            notify(f'WARNING: 100% match! Is query "{query_name}" identical to its database match, {ident}?')

        # get lineage for match
        lineage = find_match_lineage(match_ident, tax_assign,
                                    skip_idents=skip_idents,
                                    keep_full_identifiers=keep_full_identifiers,
                                    keep_identifier_versions=keep_identifier_versions)
        # ident was in skip_idents
        if not lineage:
            continue

        # summarize at rank!
        lineage = pop_to_rank(lineage, rank)
        assert lineage[-1].rank == rank, lineage[-1]
        # record info
        sum_uniq_to_query[query_name][lineage] += f_unique_to_query
        sum_uniq_weighted[query_name][lineage] += f_uniq_weighted
        sum_uniq_bp[query_name][lineage] += unique_intersect_bp

    # sort and store each as SumGathInf
    sum_uniq_to_query_sorted = []
    for query_name, lineage_weights in sum_uniq_to_query.items():
        qInfo = query_info[query_name]
        sumgather_items = list(lineage_weights.items())
        sumgather_items.sort(key = lambda x: -x[1])
        query_ani = None
        if best_only:
            lineage, fraction = sumgather_items[0]
            if fraction > 1:
                raise ValueError(f"The tax summary of query '{query_name}' is {fraction}, which is > 100% of the query!! This should not be possible. Please check that your input files come directly from a single gather run per query.")
            elif fraction == 0:
                continue
            f_weighted_at_rank = sum_uniq_weighted[query_name][lineage]
            bp_intersect_at_rank = sum_uniq_bp[query_name][lineage]
            if estimate_query_ani:
                query_ani = containment_to_distance(fraction, ksize, scaled,
                                                    n_unique_kmers= qInfo.query_hashes, sequence_len_bp= qInfo.query_bp).ani
            sres = SumGathInf(query_name, rank, fraction, lineage, qInfo.query_md5,
                                          qInfo.query_filename, f_weighted_at_rank, bp_intersect_at_rank, query_ani, qInfo.total_weighted_hashes * scaled)
            sum_uniq_to_query_sorted.append(sres)
        else:
            total_f_weighted= 0.0
            total_f_classified = 0.0
            total_bp_classified = 0
            for lineage, fraction in sumgather_items:
                query_ani = None
                if fraction > 1:
                    raise ValueError(f"The tax summary of query '{query_name}' is {fraction}, which is > 100% of the query!! This should not be possible. Please check that your input files come directly from a single gather run per query.")
                elif fraction == 0:
                    continue
                total_f_classified += fraction
                f_weighted_at_rank = sum_uniq_weighted[query_name][lineage]
                total_f_weighted += f_weighted_at_rank
                bp_intersect_at_rank = int(sum_uniq_bp[query_name][lineage])
                total_bp_classified += bp_intersect_at_rank
                if estimate_query_ani:
                    query_ani = containment_to_distance(fraction, ksize, scaled,
                                                        n_unique_kmers=qInfo.query_hashes, sequence_len_bp=qInfo.query_bp).ani
                sres = SumGathInf(query_name, rank, fraction, lineage, query_md5,
                                              query_filename, f_weighted_at_rank, bp_intersect_at_rank, query_ani, qInfo.total_weighted_hashes * scaled)
                sum_uniq_to_query_sorted.append(sres)

            # record unclassified
            lineage = ()
            query_ani = None
            fraction = 1.0 - total_f_classified
            if fraction > 0:
                f_weighted_at_rank = 1.0 - total_f_weighted
                bp_intersect_at_rank = qInfo.query_bp - total_bp_classified
                sres = SumGathInf(query_name, rank, fraction, lineage, query_md5,
                                              query_filename, f_weighted_at_rank, bp_intersect_at_rank, query_ani, qInfo.total_weighted_hashes*scaled)
                sum_uniq_to_query_sorted.append(sres)

    return sum_uniq_to_query_sorted, seen_perfect, estimate_query_ani


def find_missing_identities(gather_results, tax_assign):
    """
    Identify match ids/accessions from gather results
    that are not present in taxonomic assignments.
    """
    ident_missed= set()
    for row in gather_results:
        match_ident = row['name']
        match_ident = get_ident(match_ident)
        if match_ident not in tax_assign:
            ident_missed.add(match_ident)

    if ident_missed:
        notify(f'of {len(gather_results)} gather results, missed {len(ident_missed)} lineage assignments.')
    return ident_missed


# pass ranks; have ranks=[default_ranks]
def make_krona_header(min_rank, *, include_strain=False):
    "make header for krona output"
    header = ["fraction"]
    tl = list(taxlist(include_strain=include_strain))
    try:
        rank_index = tl.index(min_rank)
    except ValueError:
        raise ValueError(f"Rank {min_rank} not present in available ranks!")
    return tuple(header + tl[:rank_index+1])


def aggregate_by_lineage_at_rank(rank_results, *, by_query=False):
    '''
    Aggregate list of rank SumGathInfs,
    keeping query info or aggregating across queries.
    '''
    lineage_summary = defaultdict(float)
    if by_query:
        lineage_summary = defaultdict(dict)
    all_queries = []
    for res in rank_results:
        if res.query_name not in all_queries:
            all_queries.append(res.query_name)
        if by_query:
            lineage_summary[res.lineage][res.query_name] = res.fraction
        else:
            lineage_summary[res.lineage] += res.fraction
    return lineage_summary, all_queries, len(all_queries)


def format_for_krona(rank, summarized_gather):
    '''
    Aggregate list of SumGathInfs and format for krona output
    '''
    num_queries=0
    for res_rank, rank_results in summarized_gather.items():
        if res_rank == rank:
            lineage_summary, all_queries, num_queries = aggregate_by_lineage_at_rank(rank_results, by_query=False)
    # if aggregating across queries divide fraction by the total number of queries
    for lin, fraction in lineage_summary.items():
        # divide total fraction by total number of queries
        lineage_summary[lin] = fraction/num_queries

    # sort by fraction
    lin_items = list(lineage_summary.items())
    lin_items.sort(key = lambda x: -x[1])

    # reformat lineage for krona_results printing
    krona_results = []
    unclassified_fraction = 0
    for lin, fraction in lin_items:
        # save unclassified fraction for the end
        if lin == ():
            unclassified_fraction = fraction
            continue
        lin_list = display_lineage(lin).split(';')
        krona_results.append((fraction, *lin_list))

    # handle unclassified
    if unclassified_fraction:
        len_unclassified_lin = len(krona_results[-1]) -1
        unclassifed_lin = ["unclassified"]*len_unclassified_lin
        krona_results.append((unclassified_fraction, *unclassifed_lin))

    return krona_results


def write_krona(rank, krona_results, out_fp, *, sep='\t'):
    'write krona output'
    # CTB: do we want to optionally allow restriction to a specific rank
    # & above?
    header = make_krona_header(rank)
    tsv_output = csv.writer(out_fp, delimiter='\t')
    tsv_output.writerow(header)
    for res in krona_results:
        tsv_output.writerow(res)


def write_summary(summarized_gather, csv_fp, *, sep=',', limit_float_decimals=False):
    '''
    Write taxonomy-summarized gather results for each rank.
    '''
    header = SumGathInf._fields
    w = csv.DictWriter(csv_fp, header, delimiter=sep)
    w.writeheader()
    for rank, rank_results in summarized_gather.items():
        for res in rank_results:
            rD = res._asdict()
            if limit_float_decimals:
                rD['fraction'] = f'{res.fraction:.3f}'
                rD['f_weighted_at_rank'] = f'{res.f_weighted_at_rank:.3f}'
            rD['lineage'] = display_lineage(res.lineage)
            if rD['lineage'] == "":
                rD['lineage'] = "unclassified"
            w.writerow(rD)


def write_kreport(summarized_gather, csv_fp, *, sep='\t'):
    '''
    Write taxonomy-summarized gather results as kraken-style kreport.

    While this format typically records the percent of number of reads assigned to taxa,
    we can create comparable output by reporting the percent of k-mers (percent containment)
    and the total number of k-mers matched.

    standard reads-based `kreport` columns:
    - `Percent Reads Contained in Taxon`: The cumulative percentage of reads for this taxon and all descendants.
    - `Number of Reads Contained in Taxon`: The cumulative number of reads for this taxon and all descendants.
    - `Number of Reads Assigned to Taxon`: The number of reads assigned directly to this taxon (not a cumulative count of all descendants).
    - `Rank Code`: (U)nclassified, (R)oot, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. 
    - `NCBI Taxon ID`: Numerical ID from the NCBI taxonomy database.
    - `Scientific Name`: The scientific name of the taxon.

    Example reads-based `kreport` with all columns:
    ```
    88.41	2138742	193618	K	2	Bacteria
    0.16	3852	818	P	201174	  Actinobacteria
    0.13	3034	0	C	1760	    Actinomycetia
    0.13	3034	45	O	85009	      Propionibacteriales
    0.12	2989	1847	F	31957	        Propionibacteriaceae
    0.05	1142	352	G	1912216	          Cutibacterium
    0.03	790	790	S	1747	            Cutibacterium acnes
    ```

    sourmash `kreport` caveats:
    - `Percent k-mers Contained in Taxon`: weighted by k-mer abundance
    - `Estimated bp Contained in Taxon`: NOT WEIGHTED BY ABUNDANCE
    - `Number of Reads Assigned to Taxon` and `NCBI Taxon ID` will not be reported (blank entries).

    In the future, we may wish to report the NCBI taxid when we can (NCBI taxonomy only).
    '''
    columns = ["percent_containment", "num_bp_contained", "num_bp_assigned", "rank_code", "ncbi_taxid", "sci_name"]
    w = csv.DictWriter(csv_fp, columns, delimiter=sep)

    rankCode = { "superkingdom": "D", "kingdom": "K", "phylum": "P", "class": "C",
                 "order": "O", "family":"F", "genus": "G", "species": "S"} # , "": "U"

    # check - are we using v4.5.0 or later gather CSVs?
    for rank, rank_results in summarized_gather.items():
        for res in rank_results:
            if res.total_weighted_hashes == 0:
                raise ValueError("ERROR: cannot produce 'kreport' format from gather results before sourmash v4.5.0")

    unclassified_written=False
    for rank, rank_results in summarized_gather.items():
        rcode = rankCode[rank]
        for res in rank_results:
            # SummarizedGatherResults have an unclassified lineage at every rank, to facilitate reporting at a specific rank.
            # Here, we only need to report it once, since it will be the same fraction for all ranks
            if not res.lineage:
                rank_sciname = "unclassified"
                rcode = "U"
                # if we've already written the unclassified portion, skip and continue to next loop iteration
                if unclassified_written:
                    continue
                else:
                    unclassified_written=True
            else:
                rank_sciname = res.lineage[-1].name
            kresD = {"rank_code": rcode, "ncbi_taxid": "", "sci_name": rank_sciname,  "num_bp_assigned": 0}
            # total percent containment, weighted to include abundance info
            proportion = res.f_weighted_at_rank * 100
            kresD['percent_containment'] = f'{proportion:.2f}'
            # weighted bp
            kresD["num_bp_contained"] = int(res.f_weighted_at_rank * res.total_weighted_hashes)
            if rank == 'species' or rank_sciname == "unclassified":
                kresD["num_bp_assigned"] = kresD["num_bp_contained"]
            w.writerow(kresD)


def write_human_summary(summarized_gather, out_fp, display_rank):
    '''
    Write human-readable taxonomy-summarized gather results for a specific rank.
    '''
    header = SumGathInf._fields

    found_ANI = False
    results = [] 
    for rank, rank_results in summarized_gather.items():
        # only show results for a specified rank.
        if rank == display_rank:
            rank_results = list(rank_results)
            rank_results.sort(key=lambda res: -res.f_weighted_at_rank)

            for res in rank_results:
                rD = res._asdict()
                rD['fraction'] = f'{res.fraction:.3f}'
                rD['f_weighted_at_rank'] = f"{res.f_weighted_at_rank*100:>4.1f}%"
                if rD['query_ani_at_rank'] is not None:
                    found_ANI = True
                    rD['query_ani_at_rank'] = f"{res.query_ani_at_rank*100:>3.1f}%"
                else:
                    rD['query_ani_at_rank'] = '-    '
                rD['lineage'] = display_lineage(res.lineage)
                if rD['lineage'] == "":
                    rD['lineage'] = "unclassified"

                results.append(rD)


    if found_ANI:
        out_fp.write("sample name    proportion   cANI   lineage\n")
        out_fp.write("-----------    ----------   ----   -------\n")

        for rD in results:
            out_fp.write("{query_name:<15s}   {f_weighted_at_rank}     {query_ani_at_rank}  {lineage}\n".format(**rD))
    else:
        out_fp.write("sample name    proportion   lineage\n")
        out_fp.write("-----------    ----------   -------\n")

        for rD in results:
            out_fp.write("{query_name:<15s}   {f_weighted_at_rank}     {lineage}\n".format(**rD))


def write_lineage_csv(summarized_gather, csv_fp):
    '''
    Write a lineage-CSV format file suitable for use with sourmash tax ... -t.
    '''
    ranks = lca_utils.taxlist(include_strain=False)
    header = ['ident', *ranks]
    w = csv.DictWriter(csv_fp, header)
    w.writeheader()
    for rank, rank_results in summarized_gather.items():
        for res in rank_results:
            d = {}
            d[rank] = ""
            for rank, name in res.lineage:
                d[rank] = name

            d['ident'] = res.query_name
            w.writerow(d)


def write_classifications(classifications, csv_fp, *, sep=',', limit_float_decimals=False):
    '''
    Write taxonomy-classifed gather results.
    '''
    header = ClassInf._fields
    w = csv.DictWriter(csv_fp, header, delimiter=sep)
    w.writeheader()
    for rank, rank_results in classifications.items():
        for res in rank_results:
            rD = res._asdict()
            if limit_float_decimals:
                rD['fraction'] = f'{res.fraction:.3f}'
                rD['f_weighted_at_rank'] = f'{res.f_weighted_at_rank:.3f}'
            rD['lineage'] = display_lineage(res.lineage)
            # needed?
            if rD['lineage'] == "":
                rD['lineage'] = "unclassified"
            w.writerow(rD)


def combine_sumgather_csvs_by_lineage(gather_csvs, *, rank="species", accept_ranks = list(lca_utils.taxlist(include_strain=False)), force=False):
    '''
    Takes in one or more output csvs from `sourmash taxonomy summarize`
    and combines the results into a nested dictionary with lineages
    as the keys {lineage: {sample1: frac1, sample2: frac2}}.
    Uses the file basename (minus .csv extension) as sample identifier.

    usage:

        linD, all_samples = combine_sumgather_by_lineage(["sample1.csv", "sample2.csv"], rank="genus")

    output:

        linD = {lin_a: {'sample1': 0.4, 'sample2': 0.17, 'sample3': 0.6}
                lin_b: {'sample1': 0.0, 'sample2': 0.0,  'sample3': 0.1}
                lin_c: {'sample1': 0.3, 'sample2': 0.4,  'sample3': 0.2} }

        all_samples = ['sample1','sample2','sample3']

    '''
    if rank not in accept_ranks:
        raise ValueError(f"Rank {rank} not available.")

    sgD = defaultdict(dict)
    all_samples = []
    for g_csv in gather_csvs:
        # collect lineage info for this sample
        with open(g_csv, 'r') as fp:
            r = csv.DictReader(fp)
            for row in r:
                if row["rank"] == rank:
                    query_name = row["query_name"]
                    lin = row["lineage"]
                    frac = row["fraction"]
                    if query_name not in all_samples:
                        all_samples.append(query_name)
                    sgD[lin][query_name] = frac
            fp.close()
    return sgD, all_samples


def write_lineage_sample_frac(sample_names, lineage_dict, out_fp, *, format_lineage=False, sep='\t'):
    '''
    takes in a lineage dictionary with sample counts (output of combine_sumgather_by_lineage)
    and produces a tab-separated file with fractions for each sample.

    input: {lin_a: {sample1: 0.4, sample2: 0.17, sample3: 0.6}
            lin_b: {sample1: 0.0, sample2: 0.0, sample3: 0.1}
            lin_c: {sample1: 0.3, sample2: 0.4, sample3: 0.2}}

    output:

    lineage    sample1	sample2	sample3
    lin_a	  0.4    0.17     0.6
    lin_b	  0.0    0.0      0.1
    lin_c	  0.3    0.4      0.2
    '''

    header = ["lineage"] + sample_names
    w = csv.DictWriter(out_fp, header, delimiter=sep)
    w.writeheader()
    blank_row = {query_name: 0 for query_name in sample_names}
    unclassified_row = None
    for lin, sampleinfo in sorted(lineage_dict.items()):
        if format_lineage:
            lin = display_lineage(lin)

        #add lineage and 0 placeholders
        row = {'lineage': lin}
        row.update(blank_row)
        # add info for query_names that exist for this lineage
        row.update(sampleinfo)
        # if unclassified, save this row for the end
        if not lin:
            row.update({'lineage': 'unclassified'})
            unclassified_row = row
            continue
        # write row
        w.writerow(row)
    if unclassified_row:
        w.writerow(unclassified_row)


class LineageDB(abc.Mapping):
    "Base LineageDB class built around an assignments dictionary."
    def __init__(self, assign_d, avail_ranks):
        self.assignments = assign_d
        self.available_ranks = set(avail_ranks)

    def __getitem__(self, ident):
        "Retrieve the lineage tuple for identifer (or raise KeyError)"
        return self.assignments[ident]

    def __iter__(self):
        "Return all identifiers for this db."
        return iter(self.assignments)

    def __len__(self):
        "Return number of lineages"
        return len(self.assignments)

    def __bool__(self):
        "Are there any lineages at all in this database?"
        return bool(self.assignments)

    @classmethod
    def load(cls, filename, *, delimiter=',', force=False,
             keep_full_identifiers=False, keep_identifier_versions=True):
        """
        Load a taxonomy assignment CSV file into a LineageDB.

        'keep_full_identifiers=False' will split identifiers from strings
        using whitespace, e.g. 'IDENT other name stuff' => 'IDENT'

        'keep_identifier_versions=False' will remove trailing versions,
        e.g. 'IDENT.1' => 'IDENT'.
        """
        include_strain=False
        if not keep_identifier_versions and keep_full_identifiers:
            raise ValueError("keep_identifer_versions=False doesn't make sense with keep_full_identifiers=True")

        if not os.path.exists(filename):
            raise ValueError(f"'{filename}' does not exist")

        if os.path.isdir(filename):
            raise ValueError(f"'{filename}' is a directory")

        with sourmash_args.FileInputCSV(filename) as r:
            header = r.fieldnames
            if not header:
                raise ValueError(f'cannot read taxonomy assignments from {filename}')

            identifier = "ident"
            # check for ident/identifier, handle some common alternatives
            if "ident" not in header:
                # check for ident/identifier, handle some common alternatives
                if 'identifiers' in header:
                    identifier = 'identifiers'
                    header = ["ident" if "identifiers" == x else x for x in header]
                elif 'accession' in header:
                    identifier = 'accession'
                    header = ["ident" if "accession" == x else x for x in header]
                elif 'name' in header and 'lineage' in header:
                    return cls.load_from_gather_with_lineages(filename,
                                                              force=force)
                else:
                    header_str = ",".join([repr(x) for x in header])
                    raise ValueError(f'No taxonomic identifiers found; headers are {header_str}')
            # is "strain" an available rank?
            if "strain" in header:
                include_strain=True

            # check that all ranks are in header
            ranks = list(lca_utils.taxlist(include_strain=include_strain))
            if not set(ranks).issubset(header):
                # for now, just raise err if not all ranks are present.
                # in future, we can define `ranks` differently if desired
                # return them from this function so we can check the `available` ranks
                raise ValueError('Not all taxonomy ranks present')

            assignments = {}
            num_rows = 0
            n_species = 0
            n_strains = 0

            # now parse and load lineages
            for n, row in enumerate(r):
                num_rows += 1
                lineage = []
                # read row into a lineage pair
                for rank in lca_utils.taxlist(include_strain=include_strain):
                    lin = row[rank]
                    lineage.append(lca_utils.LineagePair(rank, lin))
                ident = row[identifier]

                # fold, spindle, and mutilate ident?
                ident = get_ident(ident,
                                  keep_full_identifiers=keep_full_identifiers,
                                  keep_identifier_versions=keep_identifier_versions)

                # clean lineage of null names, replace with 'unassigned'
                lineage = [ (a, lca_utils.filter_null(b)) for (a,b) in lineage ]
                lineage = [ lca_utils.LineagePair(a, b) for (a, b) in lineage ]

                # remove end nulls
                while lineage and lineage[-1].name == 'unassigned':
                    lineage = lineage[:-1]

                # store lineage tuple
                if lineage:
                    # check duplicates
                    if ident in assignments:
                        if assignments[ident] != tuple(lineage):
                            if not force:
                                raise ValueError(f"multiple lineages for identifier {ident}")
                    else:
                        assignments[ident] = tuple(lineage)

                        if lineage[-1].rank == 'species':
                            n_species += 1
                        elif lineage[-1].rank == 'strain':
                            n_species += 1
                            n_strains += 1

        return LineageDB(assignments, ranks)


    @classmethod
    def load_from_gather_with_lineages(cls, filename, *, force=False):
        """
        Load an annotated gather-with-lineages CSV file produced by
        'tax annotate' into a LineageDB.
        """
        include_strain = False

        if not os.path.exists(filename):
            raise ValueError(f"'{filename}' does not exist")

        if os.path.isdir(filename):
            raise ValueError(f"'{filename}' is a directory")

        with sourmash_args.FileInputCSV(filename) as r:
            header = r.fieldnames
            if not header:
                raise ValueError(f'cannot read taxonomy assignments from {filename}')

            if "name" not in header or "lineage" not in header:
                raise ValueError(f"Expected headers 'name' and 'lineage' not found. Is this a with-lineages file?")

            ranks = list(lca_utils.taxlist(include_strain=include_strain))
            assignments = {}
            num_rows = 0
            n_species = 0
            n_strains = 0

            # now parse and load lineages
            for n, row in enumerate(r):
                num_rows += 1

                name = row['name']
                ident = get_ident(name)
                lineage = row['lineage']
                lineage = lca_utils.make_lineage(lineage)

                # check duplicates
                if ident in assignments:
                    if assignments[ident] != tuple(lineage):
                        # this should not happen with valid
                        # sourmash tax annotate output, but check anyway.
                        if not force:
                            raise ValueError(f"multiple lineages for identifier {ident}")
                else:
                    assignments[ident] = tuple(lineage)

                    if lineage[-1].rank == 'species':
                        n_species += 1
                    elif lineage[-1].rank == 'strain':
                        n_species += 1
                        n_strains += 1

        return LineageDB(assignments, ranks)


class LineageDB_Sqlite(abc.Mapping):
    """
    A LineageDB based on a sqlite3 database with a 'sourmash_taxonomy' table.
    """
    # NOTE: 'order' is a reserved name in sql, so we have to use 'order_'.
    columns = ('superkingdom', 'phylum', 'order_', 'class', 'family',
               'genus', 'species', 'strain')
    table_name = 'sourmash_taxonomy'

    def __init__(self, conn, *, table_name=None):
        self.conn = conn

        # provide for legacy support for pre-sourmash_internal days...
        if table_name is not None:
            self.table_name = table_name

        # check that the right table is there.
        c = conn.cursor()
        try:
            c.execute(f'SELECT * FROM {self.table_name} LIMIT 1')
        except (sqlite3.DatabaseError, sqlite3.OperationalError):
            raise ValueError("not a taxonomy database")
            
        # check: can we do a 'select' on the right table?
        self.__len__()
        c = conn.cursor()

        # get available ranks...
        ranks = set()
        for column, rank in zip(self.columns, taxlist(include_strain=True)):
            query = f'SELECT COUNT({column}) FROM {self.table_name} WHERE {column} IS NOT NULL AND {column} != ""'
            c.execute(query)
            cnt, = c.fetchone()
            if cnt:
                ranks.add(rank)

        self.available_ranks = ranks
        self.cursor = c

    @classmethod
    def load(cls, location):
        "load taxonomy information from an existing sqlite3 database"
        conn = sqlite_utils.open_sqlite_db(location)
        if not conn:
            raise ValueError("not a sqlite taxonomy database")

        table_name = None
        c = conn.cursor()
        try:
            info = sqlite_utils.get_sourmash_internal(c)
        except sqlite3.OperationalError:
            info = {}

        if 'SqliteLineage' in info:
            if info['SqliteLineage'] != '1.0':
                raise IndexNotSupported

            table_name = 'sourmash_taxonomy'
        else:
            # legacy support for old taxonomy DB, pre sourmash_internal.
            try:
                c.execute('SELECT * FROM taxonomy LIMIT 1')
                table_name = 'taxonomy'
            except sqlite3.OperationalError:
                pass

        if table_name is None:
            raise ValueError("not a sqlite taxonomy database")

        return cls(conn, table_name=table_name)

    def _make_tup(self, row):
        "build a tuple of LineagePairs for this sqlite row"
        tup = [ lca_utils.LineagePair(n, r) for (n, r) in zip(taxlist(True), row) ]
        return tuple(tup)

    def __getitem__(self, ident):
        "Retrieve lineage for identifer"
        c = self.cursor
        c.execute(f'SELECT superkingdom, phylum, class, order_, family, genus, species, strain FROM {self.table_name} WHERE ident=?', (ident,))

        # retrieve names list...
        names = c.fetchone()
        if names:
            # ...and construct lineage tuple
            tup = self._make_tup(names)
            while tup and not tup[-1].name:
                tup = tup[:-1]

            return tup

        raise KeyError(ident)

    def __bool__(self):
        "Do we have any info?"
        return bool(len(self))

    def __len__(self):
        "Return number of rows"
        c = self.conn.cursor()
        c.execute(f'SELECT COUNT(DISTINCT ident) FROM {self.table_name}')
        nrows, = c.fetchone()
        return nrows

    def __iter__(self):
        "Return all identifiers"
        # create new cursor so as to allow other operations
        c = self.conn.cursor()
        c.execute(f'SELECT DISTINCT ident FROM {self.table_name}')

        for ident, in c:
            yield ident

    def items(self):
        "return all items in the sqlite database"
        c = self.conn.cursor()

        c.execute(f'SELECT DISTINCT ident, superkingdom, phylum, class, order_, family, genus, species, strain FROM {self.table_name}')

        for ident, *names in c:
            yield ident, self._make_tup(names)


class MultiLineageDB(abc.Mapping):
    "A wrapper for (dynamically) combining multiple lineage databases."

    # NTP: currently, later lineage databases will override earlier ones.
    # Do we want to report/summarize shadowed identifiers?

    def __init__(self):
        self.lineage_dbs = []

    @property
    def available_ranks(self):
        "build the union of available ranks across all databases"
        # CTB: do we need to worry about lineages of shadowed identifiers?
        x = set()
        for db in self.lineage_dbs:
            x.update(db.available_ranks)
        return x

    def add(self, db):
        "Add a new lineage database"
        self.lineage_dbs.insert(0, db)

    def __iter__(self):
        "Return all identifiers (once)"
        seen = set()
        for db in self.lineage_dbs:
            for k in db:
                if k not in seen:
                    seen.add(k)
                    yield k

    def items(self):
        "Return all (identifiers, lineage_tup), masking duplicate idents"
        seen = set()
        for db in self.lineage_dbs:
            for k, v in db.items():
                if k not in seen:
                    seen.add(k)
                    yield k, v

    def shadowed_identifiers(self):
        seen = set()
        dups = set()
        for db in self.lineage_dbs:
            for k, v in db.items():
                if k in seen:
                    dups.add(k)
                else:
                    seen.add(k)
        return seen

    def __getitem__(self, ident):
        "Return lineage tuple for first match to identifier."
        for db in self.lineage_dbs:
            if ident in db:
                return db[ident]

        # not found? KeyError!
        raise KeyError(ident)

    def __len__(self):
        "Return number of distinct identifiers. Currently iterates over all."
        # CTB: maybe we can make this unnecessary?
        x = set(self)
        return len(x)

    def __bool__(self):
        "True if any contained database has content."
        return any( bool(db) for db in self.lineage_dbs )

    def save(self, filename_or_fp, file_format):
        assert file_format in ('sql', 'csv')

        is_filename = False
        try:
            filename_or_fp.write
        except AttributeError:
            is_filename = True

        if file_format == 'sql':
            if not is_filename:
                raise ValueError("file format '{file_format}' requires a filename, not a file handle")
            self._save_sqlite(filename_or_fp)
        elif file_format == 'csv':
            # we need a file handle; open file.
            fp = filename_or_fp
            if is_filename:
                if filename_or_fp.endswith('.gz'):
                    fp = gzip.open(filename_or_fp, 'wt', newline="")
                else:
                    fp = open(filename_or_fp, 'w', newline="")

            try:
                self._save_csv(fp)
            finally:
                # close the file we opened!
                if is_filename:
                    fp.close()

    def _save_sqlite(self, filename, *, conn=None):
        from sourmash import sqlite_utils

        if conn is None:
            db = sqlite3.connect(filename)
        else:
            assert not filename
            db = conn

        cursor = db.cursor()
        try:
            sqlite_utils.add_sourmash_internal(cursor, 'SqliteLineage', '1.0')
        except sqlite3.OperationalError:
            raise ValueError("attempt to write a readonly database")

        try:
            # CTB: could add 'IF NOT EXIST' here; would need tests, too.
            cursor.execute("""

        CREATE TABLE sourmash_taxonomy (
            ident TEXT NOT NULL,
            superkingdom TEXT,
            phylum TEXT,
            class TEXT,
            order_ TEXT,
            family TEXT,
            genus TEXT,
            species TEXT,
            strain TEXT
        )
        """)
            did_create = True
        except sqlite3.OperationalError:
            # already exists?
            raise ValueError(f"taxonomy table already exists in '{filename}'")

        # follow up and create index
        cursor.execute("CREATE UNIQUE INDEX sourmash_taxonomy_ident ON sourmash_taxonomy(ident);")
        for ident, tax in self.items():
            x = [ident, *[ t.name for t in tax ]]

            # fill the taxonomy tuple with empty values until it's the
            # right length for the SQL statement -
            while len(x) < 9:
                x.append('')

            cursor.execute('INSERT INTO sourmash_taxonomy (ident, superkingdom, phylum, class, order_, family, genus, species, strain) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)', x)

        db.commit()

    def _save_csv(self, fp):
        headers = ['identifiers'] + list(taxlist(include_strain=True))
        w = csv.DictWriter(fp, fieldnames=headers)
        w.writeheader()

        for n, (ident, tax) in enumerate(self.items()):
            row = {}
            row['identifiers'] = ident

            # convert tax LineagePairs into dictionary
            for t in tax:
                row[t.rank] = t.name

            # add strain if needed
            if 'strain' not in row:
                row['strain'] = ''

            w.writerow(row)

    @classmethod
    def load(cls, locations, **kwargs):
        "Load one or more taxonomies from the given location(s)"
        force = kwargs.get('force', False)

        if isinstance(locations, str):
            raise TypeError("'locations' should be a list, not a string")

        tax_assign = cls()
        for location in locations:
            # try faster formats first
            loaded = False

            # sqlite db?
            try:
                this_tax_assign = LineageDB_Sqlite.load(location)
                loaded = True
            except ValueError:
                pass

            # CSV file?
            if not loaded:
                try:
                    this_tax_assign = LineageDB.load(location, **kwargs)
                    loaded = True
                except (ValueError, csv.Error) as exc:
                    # for the last loader, just pass along ValueError...
                    if not force:
                        raise ValueError(f"cannot read taxonomy assignments from '{location}': {str(exc)}")

            # nothing loaded, goodbye!
            if not loaded and not force:
                raise ValueError(f"cannot read taxonomy assignments from '{location}'")

            if loaded:
                tax_assign.add(this_tax_assign)

        return tax_assign


# strategy from: https://subscription.packtpub.com/book/programming/9781800207455/10/ch10lvl1sec01/using-dataclasses-to-simplify-working-with-csv-files
@dataclass
class GatherRow(): # all cols should match "gather_write_cols" in `search.py`
   # essential columns
   query_name: str
   name: str # match_name
   f_unique_weighted: float
   f_unique_to_query: float
   unique_intersect_bp: int
   remaining_bp: int
   query_md5: str
   query_filename: str
   # new essential cols: requires 4.4x
   query_bp: int
   ksize: int
   scaled: int

   # non-essential
   intersect_bp: int = None
   f_orig_query: float = None
   f_match: float = None
   average_abund: float = None
   median_abund: float = None
   std_abund: float = None
   filename: str = None
   md5: str = None
   f_match_orig: float = None
   gather_result_rank: str = None
   moltype: str = None
   query_n_hashes: int = None
   query_abundance: int = None
   query_containment_ani: float = None
   match_containment_ani: float = None
   average_containment_ani: float = None
   max_containment_ani: float = None
   potential_false_negative: bool = None
   n_unique_weighted_found: int = None
   sum_weighted_found: int = None
   total_weighted_hashes: int = None


@dataclass(frozen=True)
class QueryInfo():
   """Class for storing per-rank lineage information"""
   query_name: str
   query_md5: str
   query_filename: str
   query_bp: int
   query_n_hashes: int
   ksize: int
   scaled: int
   total_weighted_hashes: int = 0

   @property
   def total_weighted_bp(self):
       return self.total_weighted_hashes * self.scaled

@dataclass
class TaxResult():
    raw: GatherRow
    # can we get rid of these / just choose default ident hacking/slashing for future?
    keep_full_identifiers: bool = False
    keep_identifier_versions: bool = False

    query_name: str = field(init=False)
    query_info: QueryInfo = field(init=False)
    match_ident: str = field(init=False)
    lineageInfo: RankLineageInfo = RankLineageInfo()
    skipped_ident: bool = False
    missed_ident: bool = False
    match_lineage_attempted: bool = False

    def __post_init__(self):
        self.get_ident()
        self.query_name = self.raw.query_name
        if self.raw.query_n_hashes:
            self.raw.query_n_hashes = int(self.raw.query_n_hashes)
        if self.raw.total_weighted_hashes:
            self.raw.total_weighted_hashes = int(self.raw.total_weighted_hashes)
        else:
            self.raw.total_weighted_hashes = 0

        self.query_info = QueryInfo(query_name = self.raw.query_name,
                                  query_md5=self.raw.query_md5,
                                  query_filename = self.raw.query_filename,
                                  query_bp = int(self.raw.query_bp),
                                  query_n_hashes = self.raw.query_n_hashes,
                                  total_weighted_hashes = self.raw.total_weighted_hashes,
                                  ksize = int(self.raw.ksize),
                                  scaled = int(self.raw.scaled)
                                  )
        # cast and store the imp bits
        self.f_unique_to_query = float(self.raw.f_unique_to_query)
        self.f_unique_weighted = float(self.raw.f_unique_weighted)
        self.unique_intersect_bp = int(self.raw.unique_intersect_bp)

    def get_ident(self):
        # split identifiers = split on whitespace
        # keep identifiers = don't split .[12] from assembly accessions
        "Hack and slash identifiers."
        self.match_ident = self.raw.name
        if not self.keep_full_identifiers:
            self.match_ident = self.raw.name.split(' ')[0]
        else:
            #overrides version bc can't keep full without keeping version
            self.keep_identifier_versions = True
        if not self.keep_identifier_versions:
            self.match_ident = self.match_ident.split('.')[0]


    def get_match_lineage(self, tax_assignments, skip_idents=None, fail_on_missing_taxonomy=False):
        if skip_idents and self.match_ident in skip_idents:
            self.skipped_ident = True
        else:
            lin = tax_assignments.get(self.match_ident)
            if lin:
                self.lineageInfo = RankLineageInfo(lineage=lin)
            else:
                self.missed_ident=True
        self.match_lineage_attempted = True
        if self.missed_ident and fail_on_missing_taxonomy:
            raise ValueError(f"Error: ident '{self.match_ident}' is not in the taxonomy database. Failing, as requested via --fail-on-missing-taxonomy")
#            raise ValueError('Failing on missing taxonomy, as requested via --fail-on-missing-taxonomy.')

@dataclass
class SummarizedGatherResult():
#   """Class for storing summarized lineage information"""
    rank: str
    fraction: float
    lineage: RankLineageInfo
    f_weighted_at_rank: float
    bp_match_at_rank: int
    query_ani_at_rank: float

    def as_summary_dict(self, query_info=None, limit_float=False):
        sD = asdict(self)
        sD['lineage'] = self.lineage.display_lineage(null_as_unclassified=True)
        if query_info is not None:
            sD['query_name'] = query_info.query_name
            sD['query_md5'] = query_info.query_md5
            sD['query_filename'] = query_info.query_filename
            sD['total_weighted_hashes'] = query_info.total_weighted_hashes
        if limit_float:
            sD['fraction'] = f'{self.fraction:.3f}'
            sD['f_weighted_at_rank'] = f'{self.f_weighted_at_rank:.3f}'

        return(sD)

    def as_human_friendly_dict(self, query_info = None):
        sD = self.as_summary_dict(query_info=query_info, limit_float=True)
        sD['f_weighted_at_rank'] = f"{self.f_weighted_at_rank*100:>4.1f}%"
        if sD['query_ani_at_rank'] is not None:
            sD['query_ani_at_rank'] = f"{self.query_ani_at_rank*100:>3.1f}%"
        else:
            sD['query_ani_at_rank'] = '-    '
        return sD

    def as_kreport_dict(self, rankCode, total_weighted_hashes):
         sD = asdict(self)
         this_rank = self.lineage.lowest_rank
         this_sciname = self.lineage.lowest_lineage_name
         sD['ncbi_taxid'] = self.lineage.lowest_lineage_taxid
         sD['sci_name'] = this_sciname
         #sD['lineage'] = self.lineage.display_lineage(null_as_unclassified=True)
         sD['rank_code'] = rankCode[this_rank]
         sD['num_bp_assigned'] = 0
         # total percent containment, weighted to include abundance info
         sD['percent_containment'] = f'{self.f_weighted_at_rank * 100:.2f}'
         sD["num_bp_contained"] = int(self.f_weighted_at_rank * total_weighted_hashes)
         # the number of bp actually 'assigned' at this rank. Sourmash assigns everything
         # at genome level, but since kreport traditionally doesn't include 'strain' or genome,
         # it is reasonable to state that sourmash assigns at 'species' level for this.
         # can be modified later.
         lowest_assignment_rank = 'species'
         if this_rank == lowest_assignment_rank or this_sciname == "unclassified":
            sD["num_bp_assigned"] = sD["num_bp_contained"]

@dataclass
class ClassificationResult(SummarizedGatherResult):
#   """Class for storing summarized lineage information"""
    status: str


@dataclass
class QueryTaxResult():
    """Store all TaxResults for a query. Enable summarization."""
    query_info: QueryInfo # initialize with QueryInfo dataclass

    def __post_init__(self):
        self.query_name = self.query_info.query_name # for convenience
        self._init_taxresult_vars()
        self._init_summarization_vars()

    def _init_taxresult_vars(self):
        self.ranks = []
        self.raw_taxresults = []
        self.skipped_idents= set()
        self.missed_idents = set()
        self.n_missed = 0
        self.n_skipped = 0
        self.perfect_match = set()

    def _init_summarization_vars(self):
        self.sum_uniq_weighted = defaultdict(lambda: defaultdict(float))
        self.sum_uniq_to_query = defaultdict(lambda: defaultdict(float))
        self.sum_uniq_bp = defaultdict(lambda: defaultdict(int))
        self.summarized_ranks = []
        self._init_summarization_results()

    def _init_summarization_results(self):
        self.total_f_weighted = defaultdict(float) #0.0
        self.total_f_classified = defaultdict(float)#0.0
        self.total_bp_classified = defaultdict(int) #0
        self.summarized_lineage_results = defaultdict(list)

    def _init_classification_results(self):
        self.status = 'nomatch'
        self.classified_ranks = []
        self.classification_result = None
        self.krona_classification_result = None
        self.krona_unclassified_result = None
        self.krona_classification_header = []

    def is_compatible(self, taxresult):
        return taxresult.query_info == self.query_info

    @property
    def ascending_ranks(self):
        if not self.ranks:
            return []
        else:
            return self.ranks[::-1]

    def add_taxresult(self, taxresult):
        # check that all query parameters match
        if self.is_compatible(taxresult=taxresult):
            if not taxresult.match_lineage_attempted:
                raise ValueError("Error: Cannot add TaxResult. Please use get_match_lineage() to add taxonomic lineage information first.")
            if not self.ranks:
                self.ranks = taxresult.lineageInfo.ranks
            if taxresult.skipped_ident:
                self.n_skipped +=1
                self.skipped_idents.add(taxresult.match_ident)
            elif taxresult.missed_ident:
                self.n_missed +=1
                self.missed_idents.add(taxresult.match_ident)
            self.raw_taxresults.append(taxresult)
        else:
            raise ValueError("Error: Cannot add TaxResult: query information does not match.")

    def summarize_up_ranks(self, single_rank=None, force_resummarize=False):
        if self.summarized_ranks: # has already been summarized
            if force_resummarize:
                self._init_summarization_vars()
            else:
                raise ValueError("Error: already summarized using rank(s): '{', '.join(self.summarized_ranks)}'. Use 'force_resummarize=True' to reset and resummarize")
        # set ranks levels to summarize
        self.summarized_ranks = self.ascending_ranks
        if single_rank:
            if single_rank not in self.summarized_ranks:
                raise ValueError(f"Error: rank '{single_rank}' not in available ranks ({', '.join(self.summarized_ranks)})")
            self.summarized_ranks = [single_rank]
        notify(f"Starting summarization up rank(s): {', '.join(self.summarized_ranks)} ")
        for taxres in self.raw_taxresults:
            lininfo = taxres.lineageInfo
            if lininfo and lininfo.filled_lineage: # won't always have lineage to summarize (skipped idents, missed idents)
                # notify + track perfect matches
                if taxres.f_unique_to_query >= 1.0:
                    if taxres.match_ident not in self.perfect_match:
                        notify(f"WARNING: 100% match! Is query '{self.query_name}' identical to its database match, '{taxres.match_ident}'?")
                        self.perfect_match.add(taxres.match_ident)
                # add this taxresult to summary
                for rank in self.summarized_ranks:
                    if rank in lininfo.filled_ranks: # only store if this rank is filled.
                        lin_at_rank = lininfo.pop_to_rank(rank)
                        self.sum_uniq_weighted[rank][lin_at_rank] += taxres.f_unique_weighted
                        self.sum_uniq_to_query[rank][lin_at_rank] += taxres.f_unique_to_query
                        self.sum_uniq_bp[rank][lin_at_rank] += taxres.unique_intersect_bp
        # reset ranks levels to the ones that were actually summarized + that we can access for summarized result
        self.summarized_ranks = [x for x in self.summarized_ranks if x in self.sum_uniq_bp.keys()]
        if single_rank and single_rank not in self.summarized_ranks:
            raise ValueError(f"Error: rank '{single_rank}' was not available for any matching lineages.")

    def build_summarized_result(self, single_rank=None, force_resummarize=False):
        # just reset if we've already built summarized result (avoid adding to existing)? Or write in an error/force option?
        self._init_summarization_results()
        # if taxresults haven't been summarized, do that first
        if not self.summarized_ranks or force_resummarize:
            self.summarize_up_ranks(single_rank=single_rank, force_resummarize=force_resummarize)
        # catch potential error from running summarize_up_ranks separately and passing in different single_rank
        if single_rank and single_rank not in self.summarized_ranks:
            raise ValueError(f"Error: rank '{single_rank}' not in summarized rank(s), {','.join(self.summarized_ranks)}")
        # rank loop is currently done in __main__
        for rank in self.summarized_ranks[::-1]:  # reverse so that results are in descending order
            sum_uniq_to_query = self.sum_uniq_to_query[rank] #should be lineage: value
            # first, sort
            sorted_sum_uniq_to_query = list(sum_uniq_to_query.items())
            sorted_sum_uniq_to_query.sort(key = lambda x: -x[1])
            for lineage, f_unique in sorted_sum_uniq_to_query:
                query_ani = None
                if f_unique > 1:
                    raise ValueError(f"The tax summary of query '{self.query_name}' is {f_unique}, which is > 100% of the query!! This should not be possible. Please check that your input files come directly from a single gather run per query.")
                elif f_unique == 0: #no annotated results for this query. do we need to handle this differently now?
                    continue
                f_weighted_at_rank = self.sum_uniq_weighted[rank][lineage]
                bp_intersect_at_rank = self.sum_uniq_bp[rank][lineage]

                # NTP Note: These change by rank ONLY when doing best_only (selecting top hit at that particular rank)
                # now that I pulled best_only into separate fn, these don't need to be dicts...
                self.total_f_classified[rank] += f_unique
                self.total_f_weighted[rank] += f_weighted_at_rank
                self.total_bp_classified[rank] += bp_intersect_at_rank

                query_ani = containment_to_distance(f_unique, self.query_info.ksize, self.query_info.scaled,
                                                    n_unique_kmers=self.query_info.query_hashes,
                                                    sequence_len_bp=self.query_info.query_bp).ani
                sres = SummarizedGatherResult(lineage=lineage, rank=rank,
                                              f_weighted_at_rank=f_weighted_at_rank, fraction=f_unique,
                                              bp_match_at_rank=bp_intersect_at_rank, query_ani_at_rank=query_ani)
                self.summarized_lineage_results[rank].append(sres)

            # record unclassified
            lineage = ()
            query_ani = None
            f_unique = 1.0 - self.total_f_classified[rank]
            if f_unique > 0:
                f_weighted_at_rank = 1.0 - self.total_f_weighted[rank]
                bp_intersect_at_rank = self.query_info.query_bp - self.total_bp_classified[rank]
                sres = SummarizedGatherResult(lineage=lineage, rank=rank, f_weighted_at_rank=f_weighted_at_rank,
                                              fraction=f_unique, bp_match_at_rank=bp_intersect_at_rank, query_ani_at_rank=query_ani)
                self.summarized_lineage_results[rank].append(sres)

    def build_classification_result(self, rank=None, ani_threshold=None, containment_threshold=0.1, force_resummarize=False):
        if containment_threshold and not 0 <= containment_threshold <= 1:
            raise ValueError(f"Containment threshold must be between 0 and 1 (input value: {containment_threshold}).")
        if ani_threshold and not 0 <= ani_threshold <= 1:
            raise ValueError(f"ANI threshold must be between 0 and 1 (input value: {ani_threshold}).")
        self._init_classification_results() # init some fields
        #for sg in self.summarized_lineage_results:
        if not self.summarized_ranks or force_resummarize:
            self.summarize_up_ranks(single_rank=rank, force_resummarize=force_resummarize)
        # catch potential error from running summarize_up_ranks separately and passing in different single_rank
        self.classified_ranks = self.summarized_ranks
        # if a rank is provided, we need to classify ONLY using that rank
        if rank:
            if rank not in self.summarized_ranks:
                raise ValueError(f"Error: rank '{rank}' not in summarized rank(s), {','.join(self.summarized_ranks)}")
            else:
                self.classified_ranks = [rank]
        # CLASSIFY using summarization--> best only result. Best way = use ANI or containment threshold
        for this_rank in self.classified_ranks: # ascending order or just single rank
            # reset for this rank
            f_weighted=0.0
            f_unique_at_rank=0.0
            bp_intersect_at_rank=0
            #sum_uniq_to_query = self.sum_uniq_to_query[this_rank]
            sum_uniq_weighted = self.sum_uniq_weighted[this_rank]  ##SHOULD WE BE USING WEIGHTED HERE? I THINK YES, but this is a change from before.
            # sort the results and grab best
            sorted_sum_uniq_weighted = list(sum_uniq_weighted.items())
            sorted_sum_uniq_weighted.sort(key = lambda x: -x[1])
            # best only
            this_lineage, f_weighted = sorted_sum_uniq_weighted[0]
            if f_weighted > 1:
                raise ValueError(f"The tax classification of query '{self.query_name}' is {f_weighted}, which is > 100% of the query!! This should not be possible. Please check that your input files come directly from a single gather run per query.")
            elif f_weighted == 0: #no annotated results for this query. do we need to handle this differently?
                continue
            else:
                status="below_threshold"
            f_unique_at_rank = self.sum_uniq_to_query[this_rank][this_lineage]
            query_ani = containment_to_distance(f_unique_at_rank, self.query_info.ksize, self.query_info.scaled,
                                                n_unique_kmers=self.query_info.query_hashes,
                                                sequence_len_bp=self.query_info.query_bp).ani
            # set classification status based on thresholds
            if ani_threshold: # if provided, just use ani thresh
                if query_ani and query_ani >= ani_threshold:
                    status = 'match'
                else:
                    status = 'below_threshold'
            elif f_weighted >= containment_threshold:
                status = 'match'
            # determine whether to move on to a higher tax rank (if avail)
            if status == 'match':
                break

        if this_rank == "superkingdom" and status == "nomatch":
            status="below_threshold"

        # what happens when we don't have a gather match at all?
        bp_intersect_at_rank = self.sum_uniq_bp[this_rank][this_lineage]

        classif = ClassificationResult(status=status, rank=this_rank, fraction=f_unique_at_rank, lineage=this_lineage,
                                       f_weighted_at_rank=f_weighted, bp_match_at_rank=bp_intersect_at_rank,
                                       query_ani_at_rank= query_ani)
        self.classification_result = classif

        if rank is not None:
            lin_as_list = this_lineage.display_lineage().split(';')
            krona_classification = (f_weighted, *lin_as_list)
            self.krona_classification_result = (krona_classification)
            # handle unclassified - do we want/need this?
            unclassified_fraction= 1.0-f_weighted
            len_unclassified_lin = len(lin_as_list)
            unclassifed_lin = ["unclassified"]*(len_unclassified_lin)
            self.krona_unclassified_result = (unclassified_fraction, *unclassifed_lin)
            self.krona_classification_header=self.make_krona_header(min_rank=rank)

    def make_krona_header(self, min_rank):
        "make header for krona output"
        if min_rank not in self.summarized_ranks:
            raise ValueError(f"Rank '{min_rank}' not present in summarized ranks.")
        else:
            rank_index = self.ranks.index(min_rank)
        return ["fraction"] + list(self.ranks[:rank_index+1])

    def make_human_summary(self, display_rank, classification=False):
        results = []
        if classification:
            if not self.classification_result:
                raise ValueError("query not classified yet.")
            display_rank_results = [self.classification_result]
        else:
            if not self.summarized_lineage_results:
                raise ValueError("lineages not summarized yet.")
            display_rank_results = self.summarized_lineage_results[display_rank]
            display_rank_results.sort(key=lambda res: -res.f_weighted_at_rank)

        for res in display_rank_results:
            results.append(res.as_human_friendly_dict())
        return results

    def make_full_summary(self, classification=False, limit_float=False):
        results = []
        if classification:
            header= ["query_name", "status", "rank", "fraction", "lineage",
                     "query_md5", "query_filename", "f_weighted_at_rank",
                     "bp_match_at_rank", "query_ani_at_rank"]
            if not self.classification_result:
                raise ValueError("query not classified yet.")
            rD = self.classification_result.as_summary_dict(query_info = self.query_info, limit_float=limit_float)
            results.append(rD)
        else:
            header= ["query_name", "rank", "fraction", "lineage", "query_md5",
                     "query_filename", "f_weighted_at_rank", "bp_match_at_rank",
                     "query_ani_at_rank", "total_weighted_hashes"]
            if not self.summarized_lineage_results:
                raise ValueError("lineages not summarized yet.")

            for rank in self.summarized_ranks[::-1]:
                rank_results = self.summarized_lineage_results[rank]
                rank_results.sort(key=lambda res: -res.f_weighted_at_rank)
                for res in rank_results:
                    results.append(res.as_summary_dict(query_info=self.query_info, limit_float=limit_float))
        return header, results

    def make_kreport_results(self):
        rankCode = { "superkingdom": "D", "kingdom": "K", "phylum": "P", "class": "C",
                        "order": "O", "family":"F", "genus": "G", "species": "S", "unclassified": "U"}
        if self.query_info.total_weighted_hashes == 0:
            raise ValueError("ERROR: cannot produce 'kreport' format from gather results before sourmash v4.5.0")
        required_ranks = set(rankCode.keys).pop('unclassified')
        if not set(required_ranks).issubset(set(self.ranks)):
            raise ValueError("ERROR: cannot produce 'kreport' format from ranks {', '.join(self.ranks)}")
        kreport_results = []
        unclassified_recorded=False
        # want to order results descending by rank
        for rank in self.ranks:
            if rank == 'strain': # no code for strain, can't include in this output afaik
                continue
            rank_results = self.summarized_lineage_results[rank]
            for res in rank_results:
                kresD = res.as_kreport_dict(rankCode, self.query_info.total_weighted_hashes)
                if kresD['sci_name'] == "unclassified":
                    # SummarizedGatherResults have an unclassified lineage at every rank, to facilitate reporting at a specific rank.
                    # Here, we only need to report it once, since it will be the same fraction for all ranks
                    if unclassified_recorded:
                        continue
                    else:
                        unclassified_recorded = True
                kreport_results.append(kresD)
        return(kreport_results)
