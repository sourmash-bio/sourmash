"""
Taxonomic Information Classes
"""
from dataclasses import dataclass

@dataclass
class LineagePair:
    """Class for storing per-rank lineage information"""
    name: str = None
    rank: str = None

    def is_empty(self):
        return any(self.name is None, self.rank is None)

@dataclass
class LineageTuple(LineagePair):
    """Class for storing per-rank lineage information"""
    taxid: int = None # taxid allowed to be empty


@dataclass
class BaseLineageInfo:
    ranks: list = None
    lineage: tuple = None # tuple of LineageTuples/LineagePairs
    lineage_str: str = None # ';'- or ','-separated str of lineage names
    
    def __post_init__(self):
        "Initialize according to passed values"
        if self.lineage is None:
            if self.ranks is not None:
                if self.lineage_str is not None:
                        self.lineage = self.make_lineage(self.lineage_str)
                else:
                    self.lineage= self.init_empty()
            else:
                raise ValueError("Must provide ordered ranks for lineage_str")
        else:
            if self.ranks is not None:
                self.validate_lineage()
            else:
                self.ranks = [a.rank for a in self.lineage]

    def taxlist(self):
        return self.ranks
    
    def ascending_taxlist(self):
        return self.ranks[::-1]
    
    def init_empty(self):
        'initialize empty genome lineage'
        if self.lineage != None:
            raise ValueError("lineage not empty")
        for rank in self.ranks:
            self.lineage.append(LineageTuple(rank=rank))
    
    def validate_lineage(self):
        "Check if all lineage ranks are in allowed ranks; check that they are in order"
        for taxrank, lin in zip(self.ranks, self.lineage):
            if lin.rank != taxrank:
                raise ValueError(f'incomplete lineage at {taxrank} - is {lin.rank} instead')
            if lin.rank not in self.ranks:
                raise ValueError("Error: Lineage not valid. Rank {rank} not in set ranks: {self.ranks}")

    def make_lineage(self, lin):
        "Turn a ; or ,-separated set of lineages into a tuple of LineageTuple objs."
        new_lin = lin.split(';')
        if len(new_lin) == 1:
            new_lin = lin.split(',')
        new_lin = [ LineageTuple(rank=rank, name=n) for (rank, n) in zip(self.ranks, new_lin) ]
        self.lineage=tuple(new_lin)
    
    def zip_lineage(self, truncate_empty=False):
        """
        Return lineage names as a list
        """
        zipped = [a.name for a in self.lineage]
        # eliminate empty if so requested
        if truncate_empty:
            empty = ""
            last_lineage_name = zipped[-1]
            while zipped and last_lineage_name == empty:
                zipped.pop(-1)
                if zipped:
                    last_lineage_name = zipped[-1]
        return zipped

    def zip_taxid(self, truncate_empty=False):
        """
        Return taxids as a list
        """
        zipped = [a.taxid for a in self.lineage]
        # eliminate empty if so requested
        if truncate_empty:
            empty = ""
            last_lineage_taxid = zipped[-1]
            while zipped and last_lineage_taxid == empty:
                zipped.pop(-1)
                if zipped:
                    last_lineage_taxid = zipped[-1]
        return zipped

    def display_lineage(self, truncate_empty=False):
        "Return lineage names as ';'-separated list"
        return ";".join(self.zip_lineage(truncate_empty=truncate_empty))

    def display_taxid(self, truncate_empty=False):
        "Return lineage taxids as ';'-separated list"
        return ";".join(self.zip_taxid(truncate_empty=truncate_empty))

    def is_lineage_match(self, other, rank):
        """
        check to see if two lineages are a match down to given rank.
        """
        if not other.ranks == self.ranks: # check same ranks
            raise ValueError("Cannot compare lineages from taxonomies with different ranks.")
        if rank not in self.ranks: # rank is available
            raise ValueError("Desired Rank {rank} not available for this lineage")
        for a, b in zip(self.lineage, other.lineage):
            assert a.rank == b.rank
            if a.rank == rank:
                if a == b:
                    return 1
            if a != b:
                return 0

        return 0

    def pop_to_rank(self, rank):
        "Remove lineage tuples from given lineage `lin` until `rank` is reached."
        if rank not in self.ranks:
            raise ValueError("Desired Rank {rank} not available for this lineage")

        before_rank = []
        for tax_rank in self.ranks:
            if tax_rank != rank:
                before_rank.append(tax_rank)
            else:
                break        
        # are we already above rank?
        if self.lineage and self.lineage[-1].rank in before_rank:
            return tuple(self.lineage)
        
        # if not, get lineage at this rank
        while self.lineage and self.lineage[-1].rank != rank:
            self.lineage.pop()
        
        return tuple(self.lineage)
    

@dataclass
class LineageInfoRanks(BaseLineageInfo):
    """Class for storing multi-rank lineage information"""
    ranks: list = ['superkingdom', 'phylum', 'class', 'order', 'family','genus', 'species']
    include_strain: bool = False

    def __post_init__(self):
        "Initialize according to passed values"
        if self.include_strain:
            self.ranks.append("strain")
        if self.lineage is None:
            if self.lineage_str is not None:
                    self.lineage = self.make_lineage(self.lineage_str)
            else:
                self.lineage= self.init_empty()
        else:
            self.validate_lineage()
    

@dataclass
class LineageInfoLINS(BaseLineageInfo):
    """Class for storing multi-rank lineage information"""
    ranks: list
    ## WHAT special considerations do we have here?

    def __post_init__(self):
        "Initialize according to passed values"
        if self.lineage is None:
            if self.ranks is not None:
                if self.lineage_str is not None:
                        self.lineage = self.make_lineage(self.lineage_str)
                else:
                    self.lineage= self.init_empty()
            else:
                raise ValueError("Must provide ordered ranks for lineage_str")
        else:
            if self.ranks is not None:
                self.validate_lineage()
            else:
                self.ranks = [a.rank for a in self.lineage]
    
    def zip_taxid(self, truncate_empty=False):
        raise NotImplementedError
    
    def display_taxid(self, truncate_empty=False):
        raise NotImplementedError


# not sure where to go here.. we already have MultiLineagesDB .. can we use/mod that instead?
#@dataclass
#class MultiLineages:
#    # NOTE: explicitly allow any ranks so this will worth with LINS
#    """Class for manipulating groups of LineageInfo"""
#    lineages:  # list of LineageInfo??
#
#    def build_tree(self):
#        return self
#    def find_lca(self):
#        return self


@dataclass 
class QueryInfo: # prob don't need this if we just have iall info in Base gather result
    res: list = None

@dataclass 
class BaseGatherResult:
    res: list = None

@dataclass 
class SummarizedGatherResult(BaseGatherResult):
    res: list = None

@dataclass 
class ClassificationResult(BaseGatherResult):
    res: list = None

