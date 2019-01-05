"An Abstract Base Class for collections of signatures."

from abc import ABCMeta, abstractmethod
from collections import namedtuple

# @CTB copied out of search.py to deal with import order issues, #willfix
SearchResult = namedtuple('SearchResult',
                          'similarity, match_sig, md5, filename, name')


# compatible with Python 2 *and* 3:
ABC = ABCMeta("ABC", (object,), {"__slots__": ()})


class Index(ABC):
    @abstractmethod
    def find(self, search_fn, *args, **kwargs):
        """ """

    @abstractmethod
    def search(self, signature, *args, **kwargs):
        """ """

    @abstractmethod
    def gather(self, signature, *args, **kwargs):
        """ """

    @abstractmethod
    def insert(self, node):
        """ """

    @abstractmethod
    def save(self, path, storage=None, sparseness=0.0, structure_only=False):
        """ """

    @classmethod
    @abstractmethod
    def load(cls, location, leaf_loader=None, storage=None, print_version_warning=True):
        """ """


class LinearIndex(Index):
    def __init__(self):
        self.signatures = set()

    def insert(self, node):
        self.signatures.add(node)

    def find(self, search_fn, *args, **kwargs):
        matches = []

        for node in self.signatures:
            if search_fn(node, *args):
                matches.append(node)
        return matches

    def search(self, query, *args, **kwargs):
        """@@

        Note, the "best only" hint is ignored by LinearIndex.
        """

        # check arguments
        if 'threshold' not in kwargs:
            raise TypeError("'search' requires 'threshold'")
        threshold = kwargs['threshold']

        do_containment = kwargs.get('do_containment', False)
        ignore_abundance = kwargs.get('ignore_abundance', False)

        # configure search - containment? ignore abundance?
        if do_containment:
            query_match = lambda x: query.contained_by(x, downsample=True)
        else:
            query_match = lambda x: query.similarity(
                x, downsample=True, ignore_abundance=ignore_abundance)

        # do the actual search:
        matches = []

        for ss in self.signatures:
            similarity = query_match(ss)
            if similarity >= threshold:
                # @CTB: check duplicates via md5sum - here or ??
                sr = SearchResult(similarity=similarity,
                                  match_sig=ss,
                                  md5=ss.md5sum(),
                                  filename = None,
                                  name=ss.name())
                matches.append(sr)

        # @CTB sort here or ??
        return matches

    def gather(self, signature, *args, **kwargs):
        pass

    def save(self, path):
        pass

    @classmethod
    def load(cls, location):
        pass
