import sys
import io
import json
import ijson
import sourmash
from sourmash.signature import SourmashSignature
from sourmash.signature_json import (_json_next_atomic_array,
                                         _json_next_signature,
                                         load_signature_json,
                                         load_signatures_json,
                                         load_signatureset_json_iter,
                                         save_signatures_json)
from collections import OrderedDict

def test__json_next_atomic_array():
    t = (2,3,4,5,6)
    s = json.dumps(t)
    if sys.version_info[0] < 3:
        s = unicode(s)
    it = ijson.parse(io.BytesIO(s.encode('utf-8')))
    a = _json_next_atomic_array(it)
    assert len(t) == len(a)
    assert all(x == y for x,y in zip(t, a))

# integration test more than a unit test...
def test__json_next_signature():

    name = 'Foo Bar'
    filename = '/tmp/foobar'

    minhash = (2,3,4,5,6)
    t = OrderedDict((('ksize', 21),
                     ('num', len(minhash)),
                     #('md5sum', ),
                     ('cardinality', 123456),
                     ('mins', minhash)))
    s = json.dumps(t)
    if sys.version_info[0] < 3:
        s = unicode(s)
    it = ijson.parse(io.BytesIO(s.encode('utf-8')))
    # no MD5SUM
    sig = _json_next_signature(it, name, filename,
                               ignore_md5sum=True,
                               ijson=ijson)

    ## check MD5SUM
    minhash = (5,)
    t = OrderedDict((('ksize', 20),
                     ('num', len(minhash)),
                     ('md5sum', 'eae27d77ca20db309e056e3d2dcd7d69'),
                     ('cardinality', 123456),
                     ('mins', minhash)))
    s = json.dumps(t)
    if sys.version_info[0] < 3:
        s = unicode(s)
    it = ijson.parse(io.BytesIO(s.encode('utf-8')))
    sig = _json_next_signature(it, name, filename,
                               ignore_md5sum=False,
                               ijson=ijson)

# integration test more than a unit test
def test_load_signature_json():
    name = 'Foo Bar'
    filename = '/tmp/foobar'

    minhash = (2,3,4,5,6)
    t = OrderedDict((('name', name),
                     ('filename', filename),
                     ('signatures',
                      (
                          OrderedDict((('ksize', 21),
                                       ('num', len(minhash)),
                                       #('md5sum', ),
                                       ('cardinality', 123456),
                                       ('mins', minhash))),
                      ))))
    s = json.dumps(t)
    if sys.version_info[0] < 3:
        s = unicode(s)
    it = ijson.parse(io.BytesIO(s.encode('utf-8')))
    # no MD5SUM
    sig_entry = load_signature_json(it, ignore_md5sum=True)

# integration test more than a unit test
def test_load_signaturesset_json_iter():

    t = list()
    for name, filename in (('Foo', '/tmp/foo'),
                           ('Bar', '/tmp/bar')):
        minhash = (2,3,4,5,6)
        t.append(OrderedDict((
            ('class', 'sourmash_signature'),
            ('name', name),
            ('filename', filename),
            ('signatures',
             (
                 OrderedDict((('ksize', 21),
                              ('num', len(minhash)),
                              #('md5sum', ),
                              ('cardinality', 123456),
                              ('mins', minhash))),
             )))))

    s = json.dumps(t)
    if sys.version_info[0] < 3:
        s = unicode(s)
    # no MD5SUM
    sig_entries = tuple(load_signatureset_json_iter(io.BytesIO(s.encode('utf-8')),
                                                    ignore_md5sum=True,
                                                    ijson=ijson))
    assert len(sig_entries) == 2


# integration test more than a unit test
def test_load_signaturesset_json_iter_molecules():

    t = list()
    molecules = 'DNA', 'protein', 'dayhoff', 'hp'
    names = "Foo", 'Bar', "Biz", "Baz"
    filenames = '/tmp/foo', '/tmp/bar', '/tmp/biz', '/tmp/baz'

    for molecule, name, filename in zip(molecules, names, filenames):
        minhash = (2,3,4,5,6)
        t.append(OrderedDict((
            ('class', 'sourmash_signature'),
            ('name', name),
            ('filename', filename),
            ('signatures',
             (
                 OrderedDict((('ksize', 21),
                              ('num', len(minhash)),
                              #('md5sum', ),
                              ('molecule', molecule),
                              ('cardinality', 123456),
                              ('mins', minhash))),
             )))))

    s = json.dumps(t)
    if sys.version_info[0] < 3:
        s = unicode(s)
    # no MD5SUM
    sig_entries = tuple(load_signatureset_json_iter(io.BytesIO(s.encode('utf-8')),
                                                    ignore_md5sum=True,
                                                    ijson=ijson))
    # Ensure all molecule types were read properly
    assert len(sig_entries) == 4

def test_save_load_multisig_json():
    e1 = sourmash.MinHash(n=1, ksize=20)
    sig1 = SourmashSignature(e1)

    e2 = sourmash.MinHash(n=1, ksize=25)
    sig2 = SourmashSignature(e2)

    x = save_signatures_json([sig1, sig2])
    y = list(load_signatures_json(x))

    print(x)

    assert len(y) == 2
    assert sig1 in y                      # order not guaranteed, note.
    assert sig2 in y
    assert sig1 != sig2
