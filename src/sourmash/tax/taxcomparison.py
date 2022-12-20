"""
Taxonomic Information Classes
"""
from dataclasses import dataclass, field
from itertools import zip_longest

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

    def __post_init__(self):
        "Initialize according to passed values"
        if self.lineage != [LineageTuple()]:
            self.init_from_lineage()
        elif self.lineage_str is not None and self.ranks:
            self.make_lineage()
        elif self.ranks:
            self.init_empty()
        else:
            raise ValueError("Cannot initialize BaseLineageInfo. Please provide lineage or rank info.")

    def __eq__(self, other): # ignore lineage_str
        return all([self.ranks == other.ranks, self.lineage==other.lineage])

    def taxlist(self):
        return self.ranks
    
    def ascending_taxlist(self):
        return self.ranks[::-1]

    def lineageD(self):
        return {lin_tup.rank: lin_tup.name for lin_tup in self.lineage}

    def rank_index(self, rank):
        return self.ranks.index(rank)

    @property
    def filled_lineage(self):
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
            if not isinstance(lin_tup, (LineageTuple, LineagePair)):
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
            zipped = [str(a.taxid) for a in self.filled_lineage if a.taxid is not None]
        else:
            zipped = [str(a.taxid) for a in self.lineage if a.taxid is not None]
        # replace None with empty string ("")
        if None in zipped:
            zipped = ['' if x is None else x for x in zipped]
        
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


@dataclass
class RankLineageInfo(BaseLineageInfo):
    """Class for storing multi-rank lineage information"""
    ranks: list = field(default_factory=lambda: ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'])
    include_strain: bool = False

    def __post_init__(self):
        "Initialize according to passed values"
        if self.include_strain:
            self.ranks.append("strain")
        if self.lineage != [LineageTuple()]:
            self.init_from_lineage()
        elif self.lineage_str is not None:
            self.make_lineage()
        else:
            self.init_empty()

    def __eq__(self, other): # ignore lineage_str
        return all([self.ranks == other.ranks, self.lineage==other.lineage])

@dataclass
class LINSLineageInfo(BaseLineageInfo):
    """Class for storing multi-rank lineage information"""
    ranks: list
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


#@dataclass
#class QueryInfo: # prob don't need this if we just have iall info in Base gather result
#    res: list = None

#@dataclass
#class BaseGatherResult:
#    res: list = None

#@dataclass
#class SummarizedGatherResult(BaseGatherResult):
#    res: list = None

#@dataclass
#class ClassificationResult(BaseGatherResult):
#    res: list = None

