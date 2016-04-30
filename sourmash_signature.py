#! /usr/bin/env python
"""
Save and load MinHash sketches in a YAML format, along with some metadata.
"""
import yaml
import hashlib
import sourmash_lib


class SourmashSignature(object):
    "Main class for signature information."

    def __init__(self, email, estimator, name='', filename=''):
        self.d = {}
        self.d['class'] = 'sourmash_signature'
        self.d['type'] = 'mrnaseq'
        self.d['email'] = email
        if name:
            self.d['name'] = name
        if filename:
            self.d['filename'] = filename

        self.estimator = estimator

    def md5sum(self):
        "Calculate md5 hash of the bottom sketch, specifically."
        # @CTB: should include ksize, prime.
        m = hashlib.md5()
        for k in self.estimator.mh.get_mins():
            m.update(str(k).encode('utf-8'))
        return m.hexdigest()

    def name(self):
        "Return as nice a name as possible, defaulting to md5 prefix."
        # @CTB convert to printable or something.
        if 'name' in self.d:
            return self.d.get('name')
        elif 'filename' in self.d:
            return self.d.get('filename')
        else:
            return self.md5sum()[:8]

    def save(self):
        "Return metadata and a dictionary containing the sketch info."
        e = dict(self.d)
        estimator = self.estimator

        sketch = {}
        sketch['ksize'] = int(estimator.ksize)
        sketch['prime'] = estimator.p
        sketch['mins'] = list(map(int, estimator.mh.get_mins()))
        sketch['md5sum'] = self.md5sum()
        e['signature'] = sketch

        return self.d.get('email'), self.d.get('name'), \
            self.d.get('filename'), sketch

    def similarity(self, other):
        "Compute similarity with the stored MinHash."
        return self.estimator.similarity(other.estimator)


class SourmashCompositeSignature(SourmashSignature):

    def md5sum(self):
        "Calculate md5 hash of the bottom sketch, specifically."
        # @CTB: should include ksize, prime.
        m = hashlib.md5()
        for sk in self.estimator.sketches:
            for k in sk.mh.get_mins():
                m.update(str(k).encode('utf-8'))
        return m.hexdigest()

    def save(self):
        "Return metadata and a dictionary containing the sketch info."
        e = dict(self.d)
        estimator = self.estimator

        sketch = {}
        sketch['ksize'] = int(estimator.ksize)
        sketch['prime'] = estimator.p
        sketch['prefixsize'] = estimator.prefixsize
        sketch['type'] = 'composition'

        x = []
        for i, sk in enumerate(estimator.sketches):
            d = {}
            d['num'] = i
            d['mins'] = list(map(int, sk.mh.get_mins()))
            x.append(d)
        sketch['subsketches'] = x

        sketch['md5sum'] = self.md5sum()
        e['signature'] = sketch

        return self.d.get('email'), self.d.get('name'), \
            self.d.get('filename'), sketch


def load_signatures(data, select_ksize=None, ignore_md5sum=False):
    """Load a YAML file with signatures into classes.

    Returns list of SourmashSignature objects.
    """

    # record header
    d = yaml.safe_load(data)
    if d.get('class') != 'sourmash_signature':
        raise Exception("incorrect class: %s" % d.get('class'))
    email = d['email']

    name = ''
    if 'name' in d:
        name = d['name']

    filename = ''
    if 'filename' in d:
        filename = d['filename']

    # one (old format) or more (new) signatures
    if 'signature' in d:          # old format
        assert d['version'] == '0.1'
        sketch = d['signature']
        sig = _load_one_signature(sketch, email, name, filename, ignore_md5sum)

        return [sig]
    elif 'signatures' in d:       # new format
        assert d['version'] == '0.2'

        siglist = []
        for sketch in d['signatures']:
            sig = _load_one_signature(sketch, email, name, filename,
                                      ignore_md5sum)
            if not select_ksize or select_ksize == sig.estimator.ksize:
                siglist.append(sig)
        return siglist


def _load_one_signature(sketch, email, name, filename, ignore_md5sum=False):
    """Helper function to unpack and check one signature block only."""
    ksize = sketch['ksize']
    prime = sketch['prime']
    if sketch.get('type') == 'composition':
        prefixsize = sketch['prefixsize']
        n = int(sketch['subsketches']['num'])
        e = sourmash_lib.CompositionSketch(ksize=ksize, max_prime=prime,
                                           n=n, prefixsize=prefixsize)

        for item in sketch['subsketches']:
            n = item['num']
            mins = item['mins']
            n = int(n)
            for m in map(int, mins):
                e.sketches[n].mh.add_hash(m)

        sig = SourmashCompositeSignature(email, e)
    else:
        mins = list(map(int, sketch['mins']))
        n = len(mins)
        e = sourmash_lib.Estimators(ksize=ksize, max_prime=prime, n=n)
        for m in mins:
            e.mh.add_hash(m)

        sig = SourmashSignature(email, e)

    if not ignore_md5sum:
        md5sum = sketch['md5sum']
        if md5sum != sig.md5sum():
            raise Exception('error loading - md5 of estimator does not match')

    if name:
        sig.d['name'] = name
    if filename:
        sig.d['filename'] = filename

    return sig


def save_signatures(siglist):
    """Save multiple signatures into a YAML string."""
    top_records = {}
    for sig in siglist:
        email, name, filename, sketch = sig.save()
        if not email:
            raise Exception('email must be non-unique')
        k = (email, name, filename)

        x = top_records.get(k, [])
        x.append(sketch)
        top_records[k] = x

    if len(top_records) > 1:  # not yet tested
        raise Exception("no support for multiple email/name/filename yet")

    for (email, name, filename), sketches in top_records.items():
        record = {}
        record['email'] = email
        if name:
            record['name'] = name
        if filename:
            record['filename'] = filename
        record['signatures'] = sketches

        record['version'] = '0.2'
        record['class'] = 'sourmash_signature'
        record['type'] = 'mrnaseq'

        return yaml.dump(record)

    assert 0


def test_roundtrip():
    e = sourmash_lib.Estimators(n=1, ksize=20)
    e.add("AT" * 10)
    sig = SourmashSignature('titus@idyll.org', e)
    s = save_signatures([sig])
    siglist = load_signatures(s)
    sig2 = siglist[0]
    e2 = sig2.estimator

    assert sig.similarity(sig2) == 1.0
    assert sig2.similarity(sig) == 1.0


def test_md5():
    e = sourmash_lib.Estimators(n=1, ksize=20)
    e.mh.add_hash(5)
    sig = SourmashSignature('titus@idyll.org', e)
    print(sig.save())
    assert sig.md5sum() == 'e4da3b7fbbce2345d7772b0674a318d5', sig.md5sum()
