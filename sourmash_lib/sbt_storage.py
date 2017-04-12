from __future__ import print_function, unicode_literals, division

from io import BytesIO


class Storage(object):

    def save(self, path, content):
        with open(path, 'wb') as f:
            f.write(content)

    def load(self, path):
        out = BytesIO()
        with open(path, 'rb') as f:
            out.write(f.read())

        return out.getvalue()


class FSStorage(Storage):
    pass
