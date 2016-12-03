#! /usr/bin/env python
"""
Save and load MinHash sketches in a YAML format, along with some metadata.
"""
import yaml
import hashlib
import sourmash_lib

import io
import gzip
import bz2file

SIGNATURE_VERSION=0.4


class FakeHLL(object):
    def __init__(self, cardinality):
        self.cardinality = int(cardinality)

    def estimate_cardinality(self):
        return self.cardinality

    def consume_string(self):
        raise Exception("cannot add to this HLL")

    def __eq__(self, other):
        return self.cardinality == other.cardinality


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
        m = hashlib.md5()
        m.update(str(self.estimator.ksize).encode('ascii'))
        for k in self.estimator.mh.get_mins():
            m.update(str(k).encode('utf-8'))
        return m.hexdigest()

    def __eq__(self, other):
        for k in self.d:
            if self.d[k] != other.d[k]:
                return False
            
        return self.estimator == other.estimator

    def name(self):
        "Return as nice a name as possible, defaulting to md5 prefix."
        if 'name' in self.d:
            return self.d.get('name')
        elif 'filename' in self.d:
            return self.d.get('filename')
        else:
            return self.md5sum()[:8]

    def _save(self):
        "Return metadata and a dictionary containing the sketch info."
        e = dict(self.d)
        estimator = self.estimator

        sketch = {}
        sketch['ksize'] = int(estimator.ksize)
        sketch['num'] = len(estimator.mh)
        sketch['mins'] = list(map(int, estimator.mh.get_mins()))
        sketch['md5sum'] = self.md5sum()

        if estimator.mh.is_protein():
            sketch['molecule'] = 'protein'
        else:
            sketch['molecule'] = 'dna'

        if estimator.hll is not None:
            sketch['cardinality'] = estimator.hll.estimate_cardinality()

        e['signature'] = sketch

        return self.d.get('email'), self.d.get('name'), \
            self.d.get('filename'), sketch

    def similarity(self, other):
        "Compute similarity with the other MinHash signature."
        return self.estimator.similarity(other.estimator)
    jaccard = similarity


def _guess_open(filename):
    """
    Make a best-effort guess as to how to parse the given sequence file.

    Handles '-' as shortcut for stdin.
    Deals with .gz and .bz2 as well as plain text.
    """
    magic_dict = {
        b"\x1f\x8b\x08": "gz",
        b"\x42\x5a\x68": "bz2",
    }  # Inspired by http://stackoverflow.com/a/13044946/1585509

    if filename == '-':
        filename = '/dev/stdin'

    bufferedfile = io.open(file=filename, mode='rb', buffering=8192)
    num_bytes_to_peek = max(len(x) for x in magic_dict)
    file_start = bufferedfile.peek(num_bytes_to_peek)
    compression = None
    for magic, ftype in magic_dict.items():
        if file_start.startswith(magic):
            compression = ftype
            break
    if compression is 'bz2':
        sigfile = bz2file.BZ2File(filename=bufferedfile)
    elif compression is 'gz':
        if not bufferedfile.seekable():
            bufferedfile.close()
            raise ValueError("gziped data not streamable, pipe through zcat \
                            first")
        sigfile = gzip.GzipFile(filename=filename)
    else:
        sigfile = bufferedfile

    return sigfile


def load_signatures(data, select_ksize=None, select_moltype=None,
                    ignore_md5sum=False):
    """Load a YAML string with signatures into classes.

    Returns list of SourmashSignature objects.

    Note, the order is not necessarily the same as what is in the source file.
    """

    # is it a data string?
    if hasattr(data, 'find') and data.find(b'class: sourmash_signature') == -1:
        try:                                  # is it a file handle?
            data.read
        except AttributeError:                # no - treat it like a filename.
            data = _guess_open(data)

    # at this point, whatever 'data' is, it should be loadable!

    # record header
    x = yaml.load_all(data)
    siglist = []
    for n, d in enumerate(x): # allow empty records & concat of signatures
        if n > 0 and n % 100 == 0:
            print('...sig loading {}'.format(n))
        if not d:
            continue
        if d.get('class') != 'sourmash_signature':
            raise Exception("incorrect class: %s" % d.get('class'))
        email = d['email']

        name = ''
        if 'name' in d:
            name = d['name']

        filename = ''
        if 'filename' in d:
            filename = d['filename']

        if 'signatures' not in d:
            raise Exception("invalid format")

        if d['version'] != SIGNATURE_VERSION:
            raise Exception("cannot load version %s" % (d['version']))

        for sketch in d['signatures']:
            sig = _load_one_signature(sketch, email, name, filename,
                                          ignore_md5sum)
            if not select_ksize or select_ksize == sig.estimator.ksize:
                if not select_moltype or \
                     sig.estimator.is_molecule_type(select_moltype):
                    yield sig


def _load_one_signature(sketch, email, name, filename, ignore_md5sum=False):
    """Helper function to unpack and check one signature block only."""
    ksize = sketch['ksize']
    mins = list(map(int, sketch['mins']))
    n = int(sketch['num'])
    molecule = sketch.get('molecule', 'dna')
    if molecule == 'protein':
        is_protein = True
    elif molecule == 'dna':
        is_protein = False
    else:
        raise Exception("unknown molecule type: {}".format(molecule))

    e = sourmash_lib.Estimators(ksize=ksize, n=n, protein=is_protein)
    for m in mins:
        e.mh.add_hash(m)
    if 'cardinality' in sketch:
        e.hll = FakeHLL(int(sketch['cardinality']))

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


def save_signatures(siglist, fp=None):
    "Save multiple signatures into a YAML string (or into file handle 'fp')"
    top_records = {}
    for sig in siglist:
        email, name, filename, sketch = sig._save()
        k = (email, name, filename)

        x = top_records.get(k, [])
        x.append(sketch)
        top_records[k] = x

    records = []
    for (email, name, filename), sketches in top_records.items():
        record = {}
        record['email'] = email
        if name:
            record['name'] = name
        if filename:
            record['filename'] = filename
        record['signatures'] = sketches

        record['version'] = SIGNATURE_VERSION
        record['class'] = 'sourmash_signature'
        record['type'] = 'mrnaseq'
        record['hash_function'] = '0.murmur64'

        records.append(record)

    s = yaml.dump_all(records)
    if fp:
        fp.write(s.encode('utf-8'))
    if len(records):
        if fp:
            fp.write('---\n'.encode('utf-8'))
        else:
            s += '---\n'

    return s


def test_roundtrip():
    e = sourmash_lib.Estimators(n=1, ksize=20)
    e.add("AT" * 10)
    sig = SourmashSignature('titus@idyll.org', e)
    s = save_signatures([sig])
    siglist = list(load_signatures(s))
    sig2 = siglist[0]
    e2 = sig2.estimator

    assert sig.similarity(sig2) == 1.0
    assert sig2.similarity(sig) == 1.0


def test_roundtrip_empty_email():
    e = sourmash_lib.Estimators(n=1, ksize=20)
    e.add("AT" * 10)
    sig = SourmashSignature('', e)
    s = save_signatures([sig])
    siglist = list(load_signatures(s))
    sig2 = siglist[0]
    e2 = sig2.estimator

    assert sig.similarity(sig2) == 1.0
    assert sig2.similarity(sig) == 1.0


def test_md5():
    e = sourmash_lib.Estimators(n=1, ksize=20)
    e.mh.add_hash(5)
    sig = SourmashSignature('titus@idyll.org', e)
    print(sig._save())
    assert sig.md5sum() == 'eae27d77ca20db309e056e3d2dcd7d69', sig.md5sum()


def test_name():
    e = sourmash_lib.Estimators(n=1, ksize=20)
    sig = SourmashSignature('titus@idyll.org', e, name='foo')
    assert sig.name() == 'foo'


def test_name_2():
    e = sourmash_lib.Estimators(n=1, ksize=20)
    sig = SourmashSignature('titus@idyll.org', e, filename='foo.txt')
    assert sig.name() == 'foo.txt'


def test_name_3():
    e = sourmash_lib.Estimators(n=1, ksize=20)
    sig = SourmashSignature('titus@idyll.org', e, name='foo',
                            filename='foo.txt')
    assert sig.name() == 'foo'


def test_name_4():
    e = sourmash_lib.Estimators(n=1, ksize=20)
    sig = SourmashSignature('titus@idyll.org', e)
    assert sig.name() == sig.md5sum()[:8]


def test_save_load_multisig():
    e1 = sourmash_lib.Estimators(n=1, ksize=20)
    sig1 = SourmashSignature('titus@idyll.org', e1)
    
    e2 = sourmash_lib.Estimators(n=1, ksize=20)
    sig2 = SourmashSignature('titus2@idyll.org', e2)

    x = save_signatures([sig1, sig2])
    y = list(load_signatures(x))

    print(x)

    assert len(y) == 2
    assert sig1 in y                      # order not guaranteed, note.
    assert sig2 in y
    assert sig1 != sig2
