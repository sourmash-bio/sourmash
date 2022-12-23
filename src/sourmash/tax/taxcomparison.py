"""
Taxonomic Information Classes
"""
from dataclasses import dataclass, field
from itertools import zip_longest
from collections import defaultdict

from sourmash.logging import notify

from sourmash.lca import lca_utils

from sourmash.distance_utils import containment_to_distance

#from sourmash.tax.tax_utils import find_match_lineage

@dataclass(frozen=True)
class LineagePair():
    """Class for storing per-rank lineage information"""
    rank: str=None
    name: str=None
#    taxid: int=None # would actually be easier just to change it

    def is_empty(self):
        return any(self.name is None, self.rank is None)

@dataclass(frozen=True)
class LineageTuple(LineagePair):
    """Class for storing per-rank lineage information"""
    taxid: int = None # taxid allowed to be empty

    # do we want to force taxid to be int? Seems safe ish, prevents regular ranks from being input here..
    def __post_init__(self):
        if self.taxid is not None and not isinstance(self.taxid, int):
            raise TypeError('taxid not an int')


# storing lineage as list make manipulation easier within the class. BUT, do I only want to be returning
# lineage/filled_lineage as tuples, instead?
@dataclass
class BaseLineageInfo:
    ranks: list[str] = field(default_factory=lambda: [])
    lineage: list = field(default_factory=lambda: [LineageTuple()]) #list of LineageTuples/LineagePairs
    lineage_str: str = None # ';'- or ','-separated str of lineage names
    lineage_dict: dict = None # dict of rank: name
    ident: str = None

    def __post_init__(self):
        "Initialize according to passed values"
        if self.lineage != [LineageTuple()]:
            self.init_from_lineage()
        elif self.lineage_str is not None and self.ranks:
            self.make_lineage()
        elif self.lineage_dict is not None and self.ranks:
            self.init_from_lineage_dict()
        elif self.ranks:
            self.init_empty()
        else:
            raise ValueError("Cannot initialize BaseLineageInfo. Please provide lineage or rank info.")

    def __eq__(self, other): # ignore lineage_str
        return all([self.ranks == other.ranks, self.lineage==other.lineage])

    @property
    def taxlist(self):
        return self.ranks
    
    @property
    def ascending_taxlist(self):
        return self.ranks[::-1]
    
    @property
    def lowest_rank(self):
        return self.filled_ranks[-1]

    def lineageD(self):
        return {lin_tup.rank: lin_tup.name for lin_tup in self.lineage}

    def rank_index(self, rank):
        return self.ranks.index(rank)

    @property
    def filled_lineage(self):
        if not self.filled_ranks:
            return []
        # return lineage down to lowest non-empty rank. Preserves missing ranks above.
        # Would we prefer this to be the default returned by lineage??
        lowest_filled_rank_idx = self.rank_index(self.filled_ranks[-1])
        return self.lineage[:lowest_filled_rank_idx+1]
    
    def init_empty(self):
        'initialize empty genome lineage'
        if self.lineage and self.lineage != [LineageTuple()]:
            raise ValueError("lineage not empty")
        self.lineage = []
        for rank in self.ranks:
            self.lineage.append(LineageTuple(rank=rank))
        self.filled_ranks = []
       
    def init_from_lineage(self):
        'initialize from lineage tuples, allowing empty ranks and reordering if necessary'
        # first, initialize_empty
        new_lineage = []
        # check this is a list of lineage tuples:
        for lin_tup in self.lineage:
            if not isinstance(lin_tup, (LineageTuple, LineagePair, lca_utils.LineagePair)):
                raise ValueError(f"{lin_tup} is not LineageTuple or LineagePair")
        if self.ranks:
            # build empty lineage
            for rank in self.ranks:
                new_lineage.append(LineageTuple(rank=rank))
            # now add input tuples in correct spots. This corrects for order and allows empty values.
            for lin_tup in self.lineage:
                # find index for this rank
                if lin_tup.rank: # skip this tuple if rank is None or "" (empty lineage tuple. is this needed?)
                    try:
                        rank_idx = self.rank_index(lin_tup.rank)
                    except ValueError as e:
                        raise ValueError(f"Rank '{lin_tup.rank}' not present in {', '.join(self.ranks)}") from e
                    # make sure we're adding lineagetuples, not lineagePairs for consistency
                    if not isinstance(lin_tup, LineageTuple):
                        new_lineage[rank_idx] =  LineageTuple(rank=lin_tup.rank, name=lin_tup.name)
                    else:
                        new_lineage[rank_idx] =  lin_tup
            self.lineage = new_lineage
        else:
            # if ranks not provided, check that all ranks are provided; build ranks from lineage tuples
            self.ranks=[]
            for lin_tup in self.lineage:
                if not lin_tup.rank:
                    raise ValueError(f"Missing Rank for Lineage Tuple {lin_tup} and ranks not provided!")
                self.ranks.append(lin_tup.rank)
        # build list of filled ranks
        self.filled_ranks = [a.rank for a in self.lineage if a.name]

    def init_from_lineage_dict(self):
        'initialize from lineage dict, e.g. from gather csv, allowing empty ranks and reordering if necessary'
        if not isinstance(self.lineage_dict, (dict)):
            raise ValueError(f"{self.lineage_dict} is not dictionary")
        # first, initialize_empty
        new_lineage = []
        # build empty lineage
        for rank in self.ranks:
            new_lineage.append(LineageTuple(rank=rank))
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
            new_lineage[rank_idx] =  LineageTuple(rank=rank, name=name, taxid=taxid)
        self.lineage = new_lineage
        # build list of filled ranks
        self.filled_ranks = [a.rank for a in self.lineage if a.name]
        
    def make_lineage(self):
        "Turn a ; or ,-separated set of lineages into a tuple of LineageTuple objs."
        if not self.ranks:
            raise ValueError(f"Must provide ordered ranks for {self.lineage_str}")
        new_lin = self.lineage_str.split(';')
        if len(new_lin) == 1:
            new_lin = self.lineage_str.split(',')
        new_lin = [ LineageTuple(rank=rank, name=n) for (rank, n) in zip_longest(self.ranks, new_lin) ]
        self.lineage=new_lin
        self.filled_ranks = [a.rank for a in self.lineage if a.name]
    
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
        # don't turn None into str(None)
        if truncate_empty:
            zipped = [a.taxid for a in self.filled_lineage]
        else:
            zipped = [a.taxid for a in self.lineage]
        # replace None with empty string (""); cast taxids to str
        zipped = ['' if x is None else str(x) for x in zipped]
        
        return zipped

    def display_lineage(self, truncate_empty=True):
        # default truncate empty??
        "Return lineage names as ';'-separated list"
        return ";".join(self.zip_lineage(truncate_empty=truncate_empty))

    def display_taxid(self, truncate_empty=True):
        "Return lineage taxids as ';'-separated list"
        return ";".join(self.zip_taxid(truncate_empty=truncate_empty))

    def is_lineage_match(self, other, rank):
        """
        check to see if two lineages are a match down to given rank.
        """
        if not other.ranks == self.ranks: # check same ranks
            raise ValueError("Cannot compare lineages from taxonomies with different ranks.")
        if rank not in self.ranks: # rank is available
            raise ValueError(f"Desired Rank {rank} not available for this lineage")
        # always return false if rank is not filled in either of the two lineages
        if rank in self.filled_ranks and rank in other.filled_ranks:
            rank_idx = self.rank_index(rank)
            a_lin = self.lineage[:rank_idx+1]
            b_lin = other.lineage[:rank_idx+1]
            if a_lin == b_lin:
                return 1
        return 0

    def pop_to_rank(self, rank):
        # current behavior: removes information, same as original pop_to_rank. Is that desired here if we can do it another way?)
        "Remove lineage tuples from given lineage `lin` until `rank` is reached."
        if rank not in self.ranks:
            raise ValueError(f"Desired Rank {rank} not available for this lineage")
        # are we already above rank?
        if rank not in self.filled_ranks:
            return self
        # if not, pop lineage to this rank, replacing lower with empty lineagetuples
        while self.filled_ranks and self.filled_ranks[-1] != rank:
            # remove rank from filled ranks
            this_rank = self.filled_ranks[-1]
            self.filled_ranks.pop()
            # replace LineageTuple at this rank with empty rank lineage tuple
            rank_idx = self.rank_index(this_rank)
            self.lineage[rank_idx] = LineageTuple(rank=this_rank)
        
        return self
    
    def lineage_at_rank(self, rank):
        # non-descructive pop_to_rank. Returns tuple of lineagetuples
        "Remove lineage tuples from given lineage `lin` until `rank` is reached."
        if rank not in self.ranks:
            raise ValueError(f"Desired Rank {rank} not available for this lineage")
        # are we already above rank?
        if rank not in self.filled_ranks:
            return self.filled_lineage
        # if not, return lineage tuples down to desired rank
        rank_idx = self.rank_index(rank)
        return self.filled_lineage[:rank_idx+1]


@dataclass
class RankLineageInfo(BaseLineageInfo):
    """Class for storing multi-rank lineage information"""
    ranks: list = field(default_factory=lambda: ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain'])

    def __post_init__(self):
        "Initialize according to passed values"
        if self.lineage != [LineageTuple()]:
            self.init_from_lineage()
        elif self.lineage_str is not None:
            self.make_lineage()
        elif self.lineage_dict is not None:
            self.init_from_lineage_dict()
        else:
            self.init_empty()

    def __eq__(self, other):
        return all([self.ranks == other.ranks, self.lineage==other.lineage])


@dataclass
class LINSLineageInfo(BaseLineageInfo):
    """Class for storing multi-rank lineage information"""
    num_positions: int = None
    ## WHAT special considerations do we have here?
    def __post_init__(self):
        "Initialize according to passed values"
        if self.lineage != [LineageTuple()]:
            self.init_from_lineage()
        elif self.lineage_str is not None and self.ranks:
            self.make_lineage()
        elif self.ranks:
            self.init_empty()
        else:
            raise ValueError("Cannot initialize LINSLineageInfo. Please provide lineage or rank info.")
    

def build_tree(assignments, initial=None):
    """
    Builds a tree of dictionaries from lists of LineageInfo or LineagePair objects
    in 'assignments'.  This tree can then be used to find lowest common
    ancestor agreements/confusion.
    """
    if initial is None:
        tree = {}
    else:
        tree = initial

    if not assignments:
        raise ValueError("empty assignment passed to build_tree")
    
    if not isinstance(assignments, list):
        assignments = [assignments]

    for assignment in assignments:
        node = tree

        if isinstance(assignment, (BaseLineageInfo, RankLineageInfo, LINSLineageInfo)):
            assignment = assignment.filled_lineage

        for lineage_tup in assignment:
            if lineage_tup.name:
                child = node.get(lineage_tup, {})
                node[lineage_tup] = child

                # shift -> down in tree
                node = child

    return tree


# strategy from: https://subscription.packtpub.com/book/programming/9781800207455/10/ch10lvl1sec01/using-dataclasses-to-simplify-working-with-csv-files
@dataclass
class GatherRow(): # should match "gather_write_cols" in `search.py`
   # essential columns
   query_name: str
   name: str # match_name
   f_unique_weighted: float
   f_unique_to_query: float
   unique_intersect_bp: int
   remaining_bp: int
   query_md5: str
   query_filename: str
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
   query_bp: int = None
   ksize: int = None
   moltype: str = None
   scaled: int = None
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


#goal: associate all GatherRows from a particular query. Use query name to grab tax info from database.
#add functions to summarize all of these at each taxonomic rank; report (general fn's from tax_utils)
@dataclass(frozen=True)
class QueryInfo():
   """Class for storing per-rank lineage information"""
   query_name: str
   query_md5: str
   query_filename: str
   query_bp: int
   query_hashes: int
   total_weighted_hashes: int
   ksize: int
   scaled: int

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

    def __post_init__(self):
        self.get_ident()
        self.query_name = self.raw.query_name

        # before 4.x, these weren't in gather output. 
        # Get query_bp; only change type here if they exist
        if not self.raw.query_bp:
            raise ValueError("Error: Gather Results too old. Please regenerate with sourmash > 4.4.")
            #self.raw.query_bp = int(self.raw.unique_intersect_bp) + int(self.raw.remaining_bp)
        if self.raw.query_n_hashes:
            self.raw.query_n_hashes = int(self.raw.query_n_hashes)
        if self.raw.total_weighted_hashes:
            self.raw.total_weighted_hashes = int(self.raw.total_weighted_hashes)
        if self.raw.scaled:
            self.raw.scaled = int(self.raw.scaled)

        self.query_info = QueryInfo(query_name = self.raw.query_name,
                                  query_md5=self.raw.query_md5,
                                  query_filename = self.raw.query_filename,
                                  query_bp = int(self.raw.query_bp),
                                  query_hashes = self.raw.query_n_hashes,
                                  total_weighted_hashes = self.raw.total_weighted_hashes,
                                  ksize = self.raw.ksize,
                                  scaled = self.raw.scaled
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
        if self.missed_ident and fail_on_missing_taxonomy:
            raise ValueError(f"ident {self.match_ident} is not in the taxonomy database.")


# desired behavior for this: initialize once for each query,
# then add results (TaxResult) as you go, to build full representation of all results for a single query
# once all results added, write/use summarize_at_rank(self, rank) method to summarize?
# finally, add methods: write_classification_result, write_kreport, etc??


@dataclass
class SummarizedGatherResult():
#   """Class for storing summarized lineage information"""
   query_name: str
   query_md5: str
   query_filename: str
   rank: str
   fraction: float
   f_weighted_at_rank: float
   lineage: tuple
   bp_match_at_rank: int
   query_ani: float

@dataclass
class QueryTaxResult(): 
    """Store all TaxResult for a query. Enable summarization."""
    query_info: QueryInfo # initialize with QueryInfo dataclass
    query_name: str = field(init=False)
    raw_taxresults: list = field(default_factory=lambda: [])
    skipped_idents: set = field(default_factory=lambda: set()) 
    missed_idents: set = field(default_factory=lambda: set())
    sum_uniq_weighted: dict = field(default_factory=lambda: defaultdict(lambda: defaultdict(float)))
    sum_uniq_to_query: dict = field(default_factory=lambda: defaultdict(lambda: defaultdict(float)))
    sum_uniq_bp: dict = field(default_factory=lambda: defaultdict(lambda: defaultdict(int)))
    estimate_query_ani: bool = False
    perfect_match: set = field(default_factory=lambda: set())
    total_f_weighted: float = 0.0
    total_f_classified: float = 0.0
    total_bp_classified: int = 0
    summarized_lineage_results: dict = field(default_factory=lambda: defaultdict(list))
    best_only: bool = False
    ranks: list = field(default_factory=lambda: [])
    ascending_ranks: list = field(default_factory=lambda: [])
    summarized_ranks: list = field(default_factory=lambda: [])

    def __post_init__(self):
        # pull this out for convenience?
        self.query_name = self.query_info.query_name
        if self.estimate_query_ani and (not self.query_info.ksize or not self.query_info.scaled):
            self.estimate_query_ani = False
            notify("WARNING: Please run gather with sourmash >= 4.4 to estimate query ANI at rank. Continuing without ANI...")
    
    def is_compatible(self, taxresult):
        return taxresult.query_info == self.query_info
    
    def add_taxresult(self, taxresult):
        # check that all query parameters match
        if self.is_compatible(taxresult=taxresult):
            if not self.ranks:
                self.ranks = taxresult.lineageInfo.ranks
                self.ascending_ranks = taxresult.lineageInfo.ranks[::-1]
            self.raw_taxresults.append(taxresult)
            if taxresult.skipped_ident:
                self.skipped_idents.add(taxresult.match_ident)
            elif taxresult.missed_ident:
                self.missed_idents.add(taxresult.match_ident)
        else:
            #notify("this QueryInfo: {self.query_info}")
            #notify("TaxResult QueryInfo: {taxresult.query_info}")
            #notify(f"Cannot add TaxResult as query info does not match")
            raise ValueError("Cannot add TaxResult {taxresult}: Query Info does not match.")
    
    def summarize_up_ranks(self, single_rank=None):
        for taxres in self.raw_taxresults:
            lininfo = taxres.lineageInfo
            if lininfo: # we should have lineage if it wasn't in skip_idents
                if taxres.f_unique_to_query >= 1.0:
                    if taxres.match_ident not in self.perfect_match:
                        notify(f'WARNING: 100% match! Is query "{self.query_name}" identical to its database match, {taxres.match_ident}?')
                        self.perfect_match.add(taxres.match_ident)
                self.summarized_ranks = lininfo.ascending_taxlist
                if single_rank:
                    if single_rank not in self.summarized_ranks:
                        raise ValueError(f"Error: rank {single_rank} not in available ranks, {','.join(self.summarized_ranks)}")
                    self.summarized_ranks = [single_rank]
                for rank in self.summarized_ranks:
                    lin_at_rank = tuple(lininfo.lineage_at_rank(rank))
                    self.sum_uniq_weighted[rank][lin_at_rank] += taxres.f_unique_weighted
                    self.sum_uniq_to_query[rank][lin_at_rank] += taxres.f_unique_to_query
                    self.sum_uniq_bp[rank][lin_at_rank] += taxres.unique_intersect_bp
    
    def build_summarized_result(self, single_rank=None, best_only=False):
        if not self.sum_uniq_weighted: # hasn't been summarized, do that first
            self.summarize_up_ranks(single_rank=single_rank)
        # catch potential error from running summarize_up_ranks separately and passing in different single_rank
        if single_rank and single_rank not in self.summarized_ranks:
            self.best_only = True
            raise ValueError(f"Error: rank {single_rank} not in summarized ranks, {','.join(self.summarized_ranks)}")
        # rank loop is currently done in __main__
        for rank in self.summarized_ranks:  # ascending ranks or single rank
            sum_uniq_to_query = self.sum_uniq_to_query[rank] #should be lineage: value
            # first, sort
            #sorted_sum_uniq_to_query = sorted(sum_uniq_to_query.items(), lambda kv: kv[1], reverse=True)
            sorted_sum_uniq_to_query = list(sum_uniq_to_query.items())
            sorted_sum_uniq_to_query.sort(key = lambda x: -x[1])
            if best_only: # consider first match only
                sorted_sum_uniq_to_query = [sorted_sum_uniq_to_query[0]] 
            for lineage, f_unique in sorted_sum_uniq_to_query:
                query_ani = None
                if f_unique > 1:
                    raise ValueError(f"The tax summary of query '{self.query_name}' is {f_unique}, which is > 100% of the query!! This should not be possible. Please check that your input files come directly from a single gather run per query.")
                elif f_unique == 0: #no annotated results for this query. do we need to handle this differently now?
                    continue
                f_weighted_at_rank = self.sum_uniq_weighted[rank][lineage]
                bp_intersect_at_rank = self.sum_uniq_bp[rank][lineage]
                
                self.total_f_classified += f_unique
                self.total_f_weighted += f_weighted_at_rank
                self.total_bp_classified += bp_intersect_at_rank
                
                if self.estimate_query_ani:
                    query_ani = containment_to_distance(f_unique, self.query_info.ksize, self.query_info.scaled,
                                                        n_unique_kmers=self.query_info.query_hashes, 
                                                        sequence_len_bp=self.query_info.query_bp).ani
                sres = SummarizedGatherResult(query_name=self.query_name, query_md5=self.query_info.query_md5, 
                                              query_filename=self.query_info.query_filename,
                                              lineage=lineage, rank=rank,
                                              f_weighted_at_rank=f_weighted_at_rank, fraction=f_unique,
                                              bp_match_at_rank=bp_intersect_at_rank, query_ani=query_ani)
                self.summarized_lineage_results[rank].append(sres)
                #self.summarized_lineage_results.append(sres)

            # record unclassified
            lineage = ()
            query_ani = None
            f_unique = 1.0 - self.total_f_classified
            if f_unique > 0:
                f_weighted_at_rank = 1.0 - self.total_f_weighted
                bp_intersect_at_rank = self.query_info.query_bp - self.total_bp_classified
                sres = SummarizedGatherResult(query_name=self.query_name, query_md5=self.query_info.query_md5, 
                                              query_filename=self.query_info.query_filename,
                                              lineage=lineage, rank=rank,
                                              f_weighted_at_rank=f_weighted_at_rank, fraction=f_unique,
                                              bp_match_at_rank=bp_intersect_at_rank, query_ani=query_ani)
                self.summarized_lineage_results[rank].append(sres)
                #self.summarized_lineage_results.append(sres)
