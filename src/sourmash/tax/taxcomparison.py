"""
Taxonomic Information Classes
"""
import csv
from dataclasses import dataclass, field, asdict
from itertools import zip_longest
from collections import defaultdict

from sourmash.logging import notify

from sourmash.lca import lca_utils

from sourmash.distance_utils import containment_to_distance

#from sourmash.tax.tax_utils import find_match_lineage

@dataclass(frozen=True, order=True)
class LineagePair():
    """Class for storing per-rank lineage information"""
    rank: str=None
    name: str=None
#    taxid: int=None # would actually be easier just to change it

    def is_empty(self):
        return any(self.name is None, self.rank is None)

@dataclass(frozen=True, order=True)
class LineageTuple(LineagePair):
    """Class for storing per-rank lineage information"""
    taxid: int = None # taxid allowed to be empty

    # do we want to force taxid to be int? Seems safe ish, prevents regular ranks from being input here..
    def __post_init__(self):
        if self.taxid is not None and not isinstance(self.taxid, int):
            raise TypeError('taxid not an int')


# storing lineage as list make manipulation easier within the class. BUT, do I only want to be returning
# lineage/filled_lineage as tuples, instead?
@dataclass#(frozen=True, order=True)
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

    @property
    def lowest_lineage_name(self, null_as_unclassified = False):
        if not self.filled_ranks:
            if null_as_unclassified:
                return "unclassified"
            else:
                return ""
        return self.filled_lineage[-1].name
    
    @property
    def lowest_lineage_taxid(self):
        if not self.filled_ranks:
            return ""
        return self.filled_lineage[-1].taxid
    
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
        "Turn a ; or ,-separated set of lineages into a list of LineageTuple objs."
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

    def display_lineage(self, truncate_empty=True, null_as_unclassified=False):
        # default truncate empty??
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

    # Probably need to remove this/replace with `lineage_at_rank` if we need to make this class 'frozen'
    # Or make/return new instance instead of self
    def pop_to_rank(self, rank):
        # current behavior: removes information, same as original pop_to_rank. Is that desired here if we can do it another way?)
        "Remove lineage tuples from given lineage `lin` until `rank` is reached."
        if rank not in self.ranks:
            raise ValueError(f"Desired Rank '{rank}' not available for this lineage")
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
            raise ValueError(f"Desired Rank '{rank}' not available for this lineage")
        # are we already above rank?
    # if rank is not provided, we only return the filled lineage, to follow original pop_to_rank behavior.
    # Do we want to alter this behavior?
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
                                  query_hashes = self.raw.query_n_hashes,
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
    query_name: str
    rank: str
    fraction: float
    lineage: tuple
    query_md5: str
    query_filename: str
    f_weighted_at_rank: float
    bp_match_at_rank: int
    query_ani_at_rank: float
    total_weighted_hashes: int

    def as_float_limited_dict(self):
        sD = asdict(self)
        sD['fraction'] = f'{self.fraction:.3f}'
        sD['f_weighted_at_rank'] = f'{self.f_weighted_at_rank:.3f}'
        return sD

    def as_human_friendly_dict(self):
        sD = asdict(self)
        sD['fraction'] = f'{self.fraction:.3f}'
        sD['f_weighted_at_rank'] = f"{self.f_weighted_at_rank*100:>4.1f}%"
        if sD['query_ani_at_rank'] is not None:
            sD['query_ani_at_rank'] = f"{self.query_ani_at_rank*100:>3.1f}%"
        else:
            sD['query_ani_at_rank'] = '-    '
            #cD['lineage'] = display_lineage(res.lineage)
            sD['lineage'] = self.lineage.display_lineage(null_as_unclassified=True)
        return sD
    
    def as_kreport_dict(self, rankCode):
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
         sD["num_bp_contained"] = int(self.f_weighted_at_rank * self.total_weighted_hashes)
         # the number of bp actually 'assigned' at this rank. Sourmash assigns everything 
         # at genome level, but since kreport traditionally doesn't include 'strain' or genome,
         # it is reasonable to state that sourmash assigns at 'species' level for this.
         # can be modified later.
         lowest_assignment_rank = 'species'
         if this_rank == lowest_assignment_rank or this_sciname == "unclassified":
            sD["num_bp_assigned"] = sD["num_bp_contained"]

@dataclass
class ClassificationResult():
#   """Class for storing summarized lineage information"""
    query_name: str
    status: str
    rank: str
    fraction: float
    lineage: tuple
    query_md5: str
    query_filename: str
    f_weighted_at_rank: float
    bp_match_at_rank: int
    query_ani_at_rank: float

    def as_float_limited_dict(self):
        cD = asdict(self)
        cD['fraction'] = f'{self.fraction:.3f}'
        cD['f_weighted_at_rank'] = f'{self.f_weighted_at_rank:.3f}'
        return cD



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
        self.best_only_result = False
        self.total_f_weighted = defaultdict(float) #0.0
        self.total_f_classified = defaultdict(float)#0.0
        self.total_bp_classified = defaultdict(int) #0
        self.summarized_lineage_results = defaultdict(list)

    def _init_classification_results(self):
        self.best_only_result = True
        self.status = 'nomatch'
        self.classification_result = None
        self.krona_classification_result = None
    
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
                raise ValueError(f"Error: rank {single_rank} not in available ranks, {','.join(self.summarized_ranks)}")
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
                    lin_at_rank = tuple(lininfo.lineage_at_rank(rank))
                    self.sum_uniq_weighted[rank][lin_at_rank] += taxres.f_unique_weighted
                    self.sum_uniq_to_query[rank][lin_at_rank] += taxres.f_unique_to_query
                    self.sum_uniq_bp[rank][lin_at_rank] += taxres.unique_intersect_bp
     

    def build_summarized_result(self, single_rank=None, best_only=False, force_resummarize=False):
        # just reset if we've already built summarized result (avoid adding to existing)? Or write in an error/force option?
        self._init_summarization_results()
        # if taxresults haven't been summarized, do that first
        if not self.summarized_ranks or force_resummarize:
            self.summarize_up_ranks(single_rank=single_rank, force_resummarize=force_resummarize)
        # catch potential error from running summarize_up_ranks separately and passing in different single_rank
        if single_rank and single_rank not in self.summarized_ranks:
            raise ValueError(f"Error: rank '{single_rank}' not in summarized rank(s), {','.join(self.summarized_ranks)}")
        if best_only:
            self.best_only_result = True
        # rank loop is currently done in __main__
        for rank in self.summarized_ranks[::-1]:  # reverse so that results are in descending order
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
                
                # These change by rank ONLY when doing best_only (selecting top hit at that particular rank)
                self.total_f_classified[rank] += f_unique
                self.total_f_weighted[rank] += f_weighted_at_rank
                self.total_bp_classified[rank] += bp_intersect_at_rank
                
#                if self.estimate_query_ani:
                query_ani = containment_to_distance(f_unique, self.query_info.ksize, self.query_info.scaled,
                                                    n_unique_kmers=self.query_info.query_hashes,
                                                    sequence_len_bp=self.query_info.query_bp).ani
                sres = SummarizedGatherResult(query_name=self.query_name, query_md5=self.query_info.query_md5, 
                                              query_filename=self.query_info.query_filename,
                                              lineage=lineage, rank=rank,
                                              f_weighted_at_rank=f_weighted_at_rank, fraction=f_unique,
                                              bp_match_at_rank=bp_intersect_at_rank, query_ani_at_rank=query_ani,
                                              total_weighted_hashes=self.query_info.total_weighted_hashes)
                self.summarized_lineage_results[rank].append(sres)

            # record unclassified
            lineage = ()
            query_ani = None
            f_unique = 1.0 - self.total_f_classified[rank]
            if f_unique > 0:
                f_weighted_at_rank = 1.0 - self.total_f_weighted[rank]
                bp_intersect_at_rank = self.query_info.query_bp - self.total_bp_classified[rank]
                sres = SummarizedGatherResult(query_name=self.query_name, query_md5=self.query_info.query_md5, 
                                              query_filename=self.query_info.query_filename,
                                              lineage=lineage, rank=rank,
                                              f_weighted_at_rank=f_weighted_at_rank, fraction=f_unique,
                                              bp_match_at_rank=bp_intersect_at_rank, query_ani_at_rank=query_ani,
                                              total_weighted_hashes=self.query_info.total_weighted_hashes)
                self.summarized_lineage_results[rank].append(sres)

    def build_classification_result(self, rank=None, ani_threshold=None, containment_threshold=None, force_resummarize=False):
        self._init_classification_results() # init some fields
        #for sg in self.summarized_lineage_results:
        if not self.summarized_ranks or force_resummarize:
            self.summarize_up_ranks(single_rank=rank, force_resummarize=force_resummarize)
        # catch potential error from running summarize_up_ranks separately and passing in different single_rank
        if rank and rank not in self.summarized_ranks:
            raise ValueError(f"Error: rank '{rank}' not in summarized rank(s), {','.join(self.summarized_ranks)}")
        # CLASSIFY using summarization--> best only result. Best way = use ANI or containment threshold
        for this_rank in self.summarized_ranks: # ascending order or just single rank
            # reset for this rank
            f_weighted=0.0
            f_unique_at_rank=0.0
            bp_intersect_at_rank=0
            #sum_uniq_to_query = self.sum_uniq_to_query[this_rank] 
            sum_uniq_weighted = self.sum_uniq_weighted[this_rank]  ##SHOULD WE BE USING WEIGHTED HERE? I THINK MAYBE YES?
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
            if ani_threshold and query_ani and query_ani >= ani_threshold:
                status = 'match'
            elif containment_threshold and f_weighted >= containment_threshold:
                status = 'match'
            # determine whether to move on to a higher tax rank (if avail)
            if status == 'match':
                break
        
        if this_rank == "superkingdom" and status == "nomatch":
            status="below_threshold"

        # what happens when we don't have a gather match at all?
        bp_intersect_at_rank = self.sum_uniq_bp[this_rank][this_lineage]

        classif = ClassificationResult(query_name=self.query_name, status=status, rank=this_rank,
                                       fraction=f_unique_at_rank, lineage=this_lineage,
                                       query_md5=self.query_info.query_md5, query_filename=self.query_info.query_filename,
                                       f_weighted_at_rank=f_weighted, bp_match_at_rank=bp_intersect_at_rank,
                                       query_ani_at_rank= query_ani)
        # use same formatting as summarized gather result for easier writing
        self.classification_result = classif

        def write_krona_classification_result(self, rank, *, ani_threshold=None, containment_threshold=None):
            if not self.classification_result:
                self.build_classification_result(rank=rank, ani_threshold=ani_threshold, containment_threshold=containment_threshold)
            lin = self.classification_result.lineage
            #if lin == ():
            #    return (f_weighted, "unclassified")
            lin_list = lin.display_lineage().split(';') #display_lineage(sg.lineage).split(';')
            return (f_weighted, *lin_list)
        
        def make_krona_header(self, min_rank):
            "make header for krona output"
            if min_rank not in self.ranks:
                raise ValueError(f"Rank {min_rank} not present in available ranks!")
            else: 
                rank_index = self.ranks.index(min_rank)
            return ["fraction"] + self.ranks[:rank_index+1]
        
        def write_human_summary(self, display_rank, csv_fp):
            results = []
            display_rank_results = self.summarized_lineage_results[display_rank]
            display_rank_results.sort(key=lambda res: -res.f_weighted_at_rank)
            for res in display_rank_results:
                rD = res.as_human_friendly_dict()
                rD['query_name'] = self.query_name
                rD['query_md5'] = self.query_info.query_md5
                rD['query_filename'] = self.query_info.query_filename
                rD['total_weighted_hashes'] = self.query_info.total_weighted_hashes
                results.append(rD)
            return results
        
        def write_kreport(self, csv_fp, sep='\t'):
            columns = ["percent_containment", "num_bp_contained", "num_bp_assigned", "rank_code", "ncbi_taxid", "sci_name"]
            w = csv.DictWriter(csv_fp, columns, delimiter=sep)
            rankCode = { "superkingdom": "D", "kingdom": "K", "phylum": "P", "class": "C",
                         "order": "O", "family":"F", "genus": "G", "species": "S", "unclassified": "U"}
            if self.query_info.total_weighted_hashes == 0:
                raise ValueError("ERROR: cannot produce 'kreport' format from gather results before sourmash v4.5.0")
            required_ranks = set(rankCode.keys).pop('unclassified')
            if not set(required_ranks).issubset(set(self.ranks)):
                raise ValueError("ERROR: cannot produce 'kreport' format from ranks {', '.join(self.ranks)}")
            unclassified_written=False
            # want to order results descending by rank
            for rank in self.ranks:
                if rank == 'strain': # no code for strain, can't include in this output afaik
                    continue
                rank_results = self.summarized_lineage_results[rank]
                for res in rank_results:
                    kresD = res.as_kreport_dict(rankCode)
                    if kresD['sci_name'] == "unclassified":
                        # SummarizedGatherResults have an unclassified lineage at every rank, to facilitate reporting at a specific rank.
                        # Here, we only need to report it once, since it will be the same fraction for all ranks
                        if unclassified_written:
                            continue
                        else:
                            unclassified_written = True
                    w.writerow(kresD)
            


# utility functions for testing
def make_mini_taxonomy(tax_info=None):
    from sourmash.lca import lca_utils
     #usage: pass in list of tuples: [(name, lineage)]
    if not tax_info:
         tax_info = [("gA", "a;b;c")]
    taxD = {}
    for (name,lin) in tax_info:
        taxD[name] = lca_utils.make_lineage(lin)
    return taxD


def make_GatherRow(gather_dict=None):
    """Load artificial gather row (dict) into GatherRow class"""
    # default contains just the essential cols
    gatherD = {'query_name': 'q1',
               'query_md5': 'md5',
               'query_filename': 'query_fn',
               'name': 'gA',
               'f_unique_weighted': 0.2,
               'f_unique_to_query': 0.1,
               'query_bp':100,
               'unique_intersect_bp': 20,
               'remaining_bp': 1,
               'ksize': 31,
               'scaled': 1}
    if gather_dict is not None:
        gatherD.update(gather_dict)
    gatherRaw = GatherRow(**gatherD)
    return gatherRaw


def make_TaxResult(gather_dict=None, taxD=None, keep_full_ident=False, keep_ident_version=False, skip_idents=None):
    """Make TaxResult from artificial gather row (dict)"""
    gRow = make_GatherRow(gather_dict)
    taxres = TaxResult(raw=gRow, keep_full_identifiers=keep_full_ident, keep_identifier_versions=keep_ident_version)
    if taxD is not None:
        taxres.get_match_lineage(tax_assignments=taxD, skip_idents=skip_idents)
    return taxres


def make_QueryTaxResults(gather_info, taxD=None, single_query=False, keep_full_ident=False, keep_ident_version=False,
                        skip_idents=None, summarize=False, best_only=False):
    """Make QueryTaxResult(s) from artificial gather information, formatted as list of gather rows (dicts)"""
    gather_results = {}
    this_querytaxres = None
    for gather_infoD in gather_info:
        taxres = make_TaxResult(gather_infoD, taxD=taxD,  keep_full_ident=keep_full_ident,
                                keep_ident_version=keep_ident_version, skip_idents=skip_idents)
        query_name = taxres.query_name
        # add to matching QueryTaxResult or create new one
        if not this_querytaxres or not this_querytaxres.is_compatible(taxres):
            # get existing or initialize new
            this_querytaxres = gather_results.get(query_name, QueryTaxResult(taxres.query_info))
        this_querytaxres.add_taxresult(taxres)
#        print('missed_ident?', taxres.missed_ident)
        gather_results[query_name] = this_querytaxres
    if summarize:
        for query_name, qres in gather_results.items():
            qres.build_summarized_result(best_only=best_only)
    # for convenience: If working with single query, just return that QueryTaxResult.
    if single_query:
        if len(gather_results.keys()) > 1:
            raise ValueError("You passed in results for more than one query")
        else:
            return next(iter(gather_results.values()))
    return gather_results
