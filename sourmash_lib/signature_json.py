"""
Extension to sourmash.signature using JSON (making load times of collection of signatures
10 to 20 times faster). -- Laurent Gautier
"""

# This was written for Python 3, may be there is a chance it will work with Python 2...
import sys
import warnings
if sys.version_info[0] < 3:
    warnings.warn("The module 'signature_json' was written for Python 3 and you Python version is older.")

import sourmash_lib

import io
import json
import ijson

def _json_next_atomic_array(iterable, prefix_item = 'item', ijson = ijson):
    """
    - iterable: iterator as returned by ijson.parse
    - prefix_item: prefix found for items in the JSON array
    - ijson: ijson backend
    """
    l = list()
    prefix, event, value = next(iterable)
    while event != 'start_array':
        prefix, event, value = next(iterable)
    prefix, event, value = next(iterable)
    while event != 'end_array':
        assert prefix == prefix_item
        l.append(value)
        prefix, event, value = next(iterable)
    return tuple(l)

def test__json_next_atomic_array():
    t = (2,3,4,5,6)
    s = json.dumps(t)
    it = ijson.parse(io.StringIO(s))
    a = _json_next_atomic_array(it)
    assert len(t) == len(a)
    assert all(x == y for x,y in zip(t, a))

def _json_next_signature(iterable,
                         email = None,
                         name = None,
                         filename = None,
                         ignore_md5sum=False,
                         prefix_item='mins.item',
                         ijson = ijson):
    """Helper function to unpack and check one signature block only.
    - iterable: an iterable such the one returned by ijson.parse()
    - email: 
    - name:
    - filename:
    - ignore_md5sum:
    - prefix_item: required when parsing nested JSON structures
    - ijson: ijson backend to use.
    """
    from .signature import FakeHLL, SourmashSignature, SIGNATURE_VERSION

    d = dict()
    prefix, event, value = next(iterable)
    if event == 'start_map':
        prefix, event, value = next(iterable)
    while event != 'end_map':
        assert event == 'map_key'
        key = value
        if key == 'mins':
            value = _json_next_atomic_array(iterable,
                                            prefix_item=prefix_item, ijson=ijson)
        elif key == 'abundances':
            value = _json_next_atomic_array(iterable,
                                            prefix_item=prefix_item, ijson=ijson)
        else:
            prefix, event, value = next(iterable)
        d[key] = value
        prefix, event, value = next(iterable)
        
    ksize = d['ksize']
    mins = d['mins']
    n = d['num']
    e = sourmash_lib.Estimators(ksize=ksize, n=n)
    for m in mins:
        e.mh.add_hash(m)
    if 'cardinality' in d:
        e.hll = FakeHLL(d['cardinality'])

    sig = SourmashSignature(email, e)

    if not ignore_md5sum:
        md5sum = d['md5sum']
        if md5sum != sig.md5sum():
            raise Exception('error loading - md5 of estimator does not match')

    if name:
        sig.d['name'] = name
    if filename:
        sig.d['filename'] = filename

    return sig

# integration test more than a unit test...
def test__json_next_signature():
    from collections import OrderedDict
    email = 'foo@bar.com'
    name = 'Foo Bar'
    filename = '/tmp/foobar'

    minhash = (2,3,4,5,6)
    t = OrderedDict((('ksize', 21),
                     ('num', len(minhash)),
                     #('md5sum', ),
                     ('cardinality', 123456),
                     ('mins', minhash)))
    s = json.dumps(t)
    it = ijson.parse(io.StringIO(s))
    # no MD5SUM
    sig = _json_next_signature(it, email, name, filename,
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
    it = ijson.parse(io.StringIO(s))
    sig = _json_next_signature(it, email, name, filename,
                               ignore_md5sum=False,
                               ijson=ijson)
        
def load_signature_json(iterable,
                        ignore_md5sum=False,
                        prefix_item='signatures.item.mins.item',
                        ijson = ijson):
    """
    - iterable:  an iterable such as the one returned by `ijson.parse()`
    - ignore_md5sum:
    - prefix_item: prefix required to parse nested JSON structures
    - ijson: ijson backend to use
    """
    d = dict()
    prefix, event, value = next(iterable)
    if event != 'start_map':
        raise ValueError('expected "start_map".')
    
    prefix, event, value = next(iterable)
    while event != 'end_map':
        assert event == 'map_key'
        key = value
        if key == 'signatures':
            signatures = list()
            prefix, event, value = next(iterable)
            assert event == 'start_array'
            while event != 'end_array':
                sig = _json_next_signature(iterable,
                                           email = None,
                                           name = None,
                                           filename = None,
                                           ignore_md5sum=ignore_md5sum,
                                           prefix_item=prefix_item,
                                           ijson=ijson)
                signatures.append(sig)
                prefix, event, value = next(iterable)
            value = signatures
        else:
            prefix, event, value = next(iterable)
        d[key] = value
        prefix, event, value = next(iterable)

    # email, name, and filename not assumed to be parsed before the 'signatures'
    for sig in signatures:
        sig.d['email'] = d['email']
        if 'name' in d:
            sig.d['name'] = d['name']
        if 'filename' in d:
            sig.d['filename'] = d['filename']

    return d


# integration test more than a unit test
def test_load_signature_json():
    from collections import OrderedDict
    email = 'foo@bar.com'
    name = 'Foo Bar'
    filename = '/tmp/foobar'

    minhash = (2,3,4,5,6)
    t = OrderedDict((('email', email),
                     ('name', name),
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
    it = ijson.parse(io.StringIO(s))
    # no MD5SUM
    sig_entry = load_signature_json(it, ignore_md5sum=True)
        

def load_signatureset_json_iter(data, select_ksize=None, ignore_md5sum=False, ijson=ijson):
    """
    - data: file handle (or file handle-like) object
    - select_ksize:
    - ignore_md5sum:
    - ijson: ijson backend    
    """
    
    parser = ijson.parse(data)

    prefix, event, value = next(parser)
    assert prefix == '' and event == 'start_array' and value is None

    siglist = []    
    n = 0
    while True:
        try:
            sig = load_signature_json(parser,
                                      prefix_item = 'item.signatures.item.mins.item',
                                      ignore_md5sum=ignore_md5sum,
                                      ijson=ijson)
            if not select_ksize or select_ksize == sig.estimator.ksize:
                yield sig
        except ValueError:
            # possible end of the array of signatures
            prefix, event, value = next(parser)
            assert event == 'end_array'
            break
        n += 1

def load_signatures_json(data, select_ksize=None, ignore_md5sum=True, ijson=ijson):
    """
    - data: file handle (or file handle-like) object
    - select_ksize:
    - ignore_md5sum:
    - ijson: ijson backend
    """
    n = 0

    if isinstance(data, str):
        data = io.StringIO(data)
        
    it = load_signatureset_json_iter(data, select_ksize=select_ksize,
                                     ignore_md5sum=ignore_md5sum,
                                     ijson=ijson)

    for n, sigset in enumerate(it):
        if n > 0 and n % 100 == 0:
            print('\r...sig loading {:,}'.format(n), end='', flush=True)
        for sig in sigset['signatures']:
            yield sig
            
    if n > 1:
        print('\r...sig loading {:,}'.format(n), flush=True)

    
# integration test more than a unit test
def test_load_signaturesset_json_iter():
    from collections import OrderedDict

    t = list()
    for email, name, filename in (('foo@foo.com', 'Foo', '/tmp/foo'),
                                  ('bar@bar.com', 'Bar', '/tmp/bar')):
        minhash = (2,3,4,5,6)
        t.append(OrderedDict((
            ('class', 'sourmash_signature'),
            ('email', email),
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
    # no MD5SUM
    sig_entries = tuple(load_signatureset_json_iter(io.StringIO(s),
                                                    ignore_md5sum=True,
                                                    ijson=ijson))
    assert len(sig_entries) == 2
        
        
def save_signatures_json(siglist, fp=None, indent=4, sort_keys=True):
    """ Save multiple signatures into a JSON string (or into file handle 'fp')
    - siglist: sequence of SourmashSignature objects
    - fp:
    - indent: indentation spaces (an integer) or if None no indentation
    - sort_keys: sort the keys in mappings before writting to JSON
    """
    from .signature import SIGNATURE_VERSION

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

    if fp:
        s = json.dump(records, fp, indent=indent, sort_keys=sort_keys)
    else:
        s = json.dumps(records, indent=indent, sort_keys=sort_keys)

    return s


def test_save_load_multisig_json():
    from .signature import SourmashSignature

    e1 = sourmash_lib.Estimators(n=1, ksize=20)
    sig1 = SourmashSignature('lalala@land.org', e1)

    e2 = sourmash_lib.Estimators(n=1, ksize=20)
    sig2 = SourmashSignature('lalala2@land.org', e2)

    x = save_signatures_json([sig1, sig2])
    y = list(load_signatures_json(x))

    print(x)

    assert len(y) == 2
    assert sig1 in y                      # order not guaranteed, note.
    assert sig2 in y
    assert sig1 != sig2
