from __future__ import print_function, unicode_literals, division

import abc
from io import BytesIO
import os
import sys
import tarfile
import zipfile


class Storage(abc.ABCMeta(str('ABC'), (object,), {'__slots__': ()})):
    # this weird baseclass is compatible with Python 2 *and* 3,
    # we can remove once we support only py3.4+

    @abc.abstractmethod
    def save(self, path, content):
        pass

    @abc.abstractmethod
    def load(self, path):
        pass

    def init_args(self):
        return {}

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def close(self):
        pass

    def can_open(self, location):
        return False


class FSStorage(Storage):

    def __init__(self, location, subdir):
        self.location = location
        self.subdir = subdir

        fullpath = os.path.join(location, subdir)
        if not os.path.exists(fullpath):
            os.makedirs(fullpath)

    def init_args(self):
        return {'path': self.subdir}

    def save(self, path, content):
        "Save a node/leaf."
        with open(os.path.join(self.location, self.subdir, path), 'wb') as f:
            f.write(content)

        return path

    def load(self, path):
        out = BytesIO()
        with open(os.path.join(self.location, self.subdir, path), 'rb') as f:
            out.write(f.read())

        return out.getvalue()


class TarStorage(Storage):

    def __init__(self, path=None):
        # TODO: leave it open, or close/open every time?

        if path is None:
            # TODO: Open a temporary file?
            pass                          # CTB: should raise an exception, no?

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

    def load(self, path):
        content = self.tarfile.getmember(path)
        f = self.tarfile.extractfile(content)
        return f.read()

    def init_args(self):
        return {'path': self.path}

    def __exit__(self, type, value, traceback):
        self.tarfile.close()


class ZipStorage(Storage):

    def __init__(self, path=None):
        # TODO: leave it open, or close/open every time?

        if path is None:
            # TODO: Open a temporary file?
            pass

        self.path = os.path.abspath(path)

        dirname = os.path.dirname(self.path)
        if not os.path.exists(dirname):
            os.makedirs(dirname)

        # Turns out we can't delete/modify an entry in a zipfile easily.
        # For now, if the file already exists open it in read mode,
        # otherwise open in write mode.
        # This causes issues when calling the `sourmash index` command
        # many times in a row, because it will try to write in a read-only file...
        # Opening in append mode is an alternative but is misleading, because
        # duplicated entries are written again in the file, generating
        # potentially large zip files...
        # More info: https://bugs.python.org/issue6818
        if os.path.exists(self.path):
            self.zipfile = zipfile.ZipFile(path, 'r')
        else:
            self.zipfile = zipfile.ZipFile(path, mode='w',
                                           compression=zipfile.ZIP_STORED)

        self.subdir = None
        subdirs = [f for f in self.zipfile.namelist() if f.endswith("/")]
        if len(subdirs) == 1:
            self.subdir = subdirs[0]

    def save(self, path, content):
        # Sigh. SIGH. When we are Python 3.6+, remove this mode madness.
        mode = "w"
        if sys.version_info[:3] < (3, 6):
            mode = "U"
        try:
            with self.zipfile.open(path, mode=mode) as entry:
                entry.write(content)
        except KeyError:
            self.zipfile.writestr(path, content)

        # After 3.6+, we can do just this:
        #with self.zipfile.open(path, mode='w') as entry:
        #    entry.write(content)

        return path

    def load(self, path):
        try:
            return self.zipfile.read(path)
        except KeyError:
            path = os.path.join(self.subdir, path)
            return self.zipfile.read(path)

    def init_args(self):
        return {'path': self.path}

    def __exit__(self, type, value, traceback):
        self.close()

    def close(self):
        self.zipfile.close()

    @staticmethod
    def can_open(location):
        return zipfile.is_zipfile(location)

    def list_sbts(self):
        return [f for f in self.zipfile.namelist() if f.endswith(".sbt.json")]


class IPFSStorage(Storage):

    def __init__(self, pin_on_add=True, **kwargs):
        import ipfshttpclient
        self.ipfs_args = kwargs
        self.pin_on_add = pin_on_add
        self.api = ipfshttpclient.connect(**self.ipfs_args)

    def save(self, path, content):
        # api.add_bytes(b"Mary had a little lamb")
        new_obj = self.api.add_bytes(content)
        if self.pin_on_add:
            self.api.pin.add(new_obj)
        return new_obj

        # TODO: the above solution is quick and dirty.
        # we actually want something more organized,
        # like putting all the generated objects inside the same dir.
        # Check this call using the files API for an example.
        # api.files_write("/test/file", io.BytesIO(b"hi"), create=True)

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

    def __init__(self, **kwargs):
        import redis
        self.redis_args = kwargs
        self.conn = redis.Redis(**self.redis_args)

    def save(self, path, content):
        if not isinstance(content, bytes):
            content = bytes(content)
        self.conn.set(path, content)
        return path

    def load(self, path):
        return self.conn.get(path)

    def init_args(self):
        # TODO: do we want to remove stuff like password from here?
        return self.redis_args

    def __exit__(self, type, value, traceback):
        pass
