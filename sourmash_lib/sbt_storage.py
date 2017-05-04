from __future__ import print_function, unicode_literals, division

from io import BytesIO
import operator
import os
import tarfile

from cachetools import LFUCache, cachedmethod


DEFAULT_CACHESIZE = 256


class Storage(object):

    def __init__(self, cachesize=DEFAULT_CACHESIZE):
        self._cache = LFUCache(maxsize=cachesize)

    def save(self, path, content):
        raise NotImplementedError

    @cachedmethod(operator.attrgetter('_cache'))
    def load(self, path):
        raise NotImplementedError

    def init_args(self):
        return {}

    def __enter__(self):
        return self

    def __exit__(self, typ, value, traceback):
        pass


class FSStorage(Storage):

    def __init__(self, path, cachesize=DEFAULT_CACHESIZE):
        super(FSStorage, self).__init__(cachesize=cachesize)

        self.path = path
        if not os.path.exists(path):
            os.makedirs(path)

    def init_args(self):
        return {'path': self.path}

    def save(self, path, content):
        with open(os.path.join(self.path, path), 'wb') as f:
            f.write(content)

        return path

    @cachedmethod(operator.attrgetter('_cache'))
    def load(self, path):
        out = BytesIO()
        with open(os.path.join(self.path, path), 'rb') as f:
            out.write(f.read())

        return out.getvalue()


class TarStorage(Storage):

    def __init__(self, path=None, cachesize=DEFAULT_CACHESIZE):
        super(TarStorage, self).__init__(cachesize=cachesize)

        # TODO: leave it open, or close/open every time?

        if path is None:
            # TODO: Open a temporary file?
            pass

        self.path = os.path.abspath(path)

        dirname = os.path.dirname(self.path)
        if not os.path.exists(dirname):
            os.makedirs(dirname)

        if os.path.exists(self.path):
            self.tarfile = tarfile.open(path, 'r')
        else:
            self.tarfile = tarfile.open(path, 'w:gz')

    def save(self, path, content):
        info = tarfile.TarInfo(path)
        info.size = len(content)

        # TODO: check tarfile mode, if read-only reopen as writable
        self.tarfile.addfile(info, BytesIO(content))

        return path

    @cachedmethod(operator.attrgetter('_cache'))
    def load(self, path):
        content = self.tarfile.getmember(path)
        f = self.tarfile.extractfile(content)
        return f.read()

    def init_args(self):
        return {'path': self.path}

    def __exit__(self, type, value, traceback):
        self.tarfile.close()


class IPFSStorage(Storage):

    def __init__(self, cachesize=DEFAULT_CACHESIZE, **kwargs):
        super(IPFSStorage, self).__init__(cachesize=cachesize)

        import ipfsapi
        self.ipfs_args = kwargs
        self.api = ipfsapi.connect(**self.ipfs_args)

    def save(self, path, content):
        # api.add_bytes(b"Mary had a little lamb")
        return self.api.add_bytes(content)

        # TODO: the above solution is quick and dirty.
        # we actually want something more organized,
        # like putting all the generated objects inside the same dir.
        # Check this call using the files API for an example.
        # api.files_write("/test/file", io.BytesIO(b"hi"), create=True)

    @cachedmethod(operator.attrgetter('_cache'))
    def load(self, path):
        return self.api.cat(path)

    def init_args(self):
        return self.ipfs_args

    def __exit__(self, type, value, traceback):
        # TODO: do nothing for now,
        # but we actually want something more organized,
        # like putting all the generated objects inside the same dir.
        # Use the files API,
        # add files without flush(),
        # and then flush it here?
        pass


class RedisStorage(Storage):

    def __init__(self, cachesize=DEFAULT_CACHESIZE, **kwargs):
        super(RedisStorage, self).__init__(cachesize=cachesize)

        import redis
        self.redis_args = kwargs
        self.conn = redis.Redis(**self.redis_args)

    def save(self, path, content):
        self.conn.set(path, content)
        return path

    @cachedmethod(operator.attrgetter('_cache'))
    def load(self, path):
        return self.conn.get(path)

    def init_args(self):
        # TODO: do we want to remove stuff like password from here?
        return self.redis_args

    def __exit__(self, type, value, traceback):
        pass
