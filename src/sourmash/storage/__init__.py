from fnmatch import fnmatch

from sourmash.logging import notify


class Handler:
    def __init__(self, name, matcher_fn, priority, loader_fn):
        self.name = name
        self.matcher_fn = matcher_fn
        self.priority = priority
        self.loader_fn = loader_fn

    def __cmp__(self, other):
        return cmp(self.priority, other.priority)

class GlobHandler(Handler):
    def __init__(self, name, pattern, priority, loader_fn):
        def matcher_fn(filename):
            if fnmatch(filename, pattern):
                return True
            return False

        super().__init__(name, matcher_fn, priority, loader_fn)

handlers = []
def add_handler(name, matcher_fn, priority, loader_fn):
    h = Handler(name, matcher_fn, priority, loader_fn)
    add_handler_obj(h)


def add_handler_obj(h):
    handlers.append(h)


def load(filename, traverse_yield_all, cache_size=None):
    matching_handlers = [ h for h in handlers if h.matcher_fn(filename) ]
    matching_handlers.sort()

    for h in matching_handlers:
        notify(f'trying: {h.name}, {h.priority}')
        try:
            db, dbtype = h.loader_fn(filename, traverse_yield_all,
                                     cache_size=cache_size)
            if db:
                notify(f'success!')
                return db, dbtype
        except:
            # @CTB do some kind of logging here!
            import traceback
            traceback.print_exc()
            notify(f'loader {h.name} failed; trying next')

    return None, None

