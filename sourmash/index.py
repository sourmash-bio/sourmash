from abc import ABCMeta, abstractmethod


# compatible with Python 2 *and* 3:
ABC = ABCMeta("ABC", (object,), {"__slots__": ()})


class Index(ABC):
    @abstractmethod
    def find(self, search_fn, *args, **kwargs):
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

    def save(self, path):
        pass

    @classmethod
    def load(cls, location):
        pass
