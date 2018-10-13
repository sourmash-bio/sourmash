from abc import ABCMeta, abstractmethod


# compatible with Python 2 *and* 3:
ABC = ABCMeta('ABC', (object,), {'__slots__': ()})


class Index(ABC):

    @abstractmethod
    def find(self, search_fn, *args, **kwargs):
        ''' '''

    @abstractmethod
    def insert(self, node):
        ''' '''

    @abstractmethod
    def save(self, path, storage=None, sparseness=0.0, structure_only=False):
        ''' '''

    @classmethod
    @abstractmethod
    def load(cls, location, leaf_loader=None, storage=None, print_version_warning=True):
        ''' '''
