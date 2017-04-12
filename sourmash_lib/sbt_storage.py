from __future__ import print_function, unicode_literals, division

import khmer


class Storage(object):

    def save(self, path, content):
        # TODO: this works for khmer nodetables only
        content.save(path)

    def load(self, path):
        return khmer.load_nodegraph(path)


class FSStorage(Storage):
    pass
