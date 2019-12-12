"""
Extension to sourmash.signature using JSON (making load times of collection of signatures
10 to 20 times faster). -- Laurent Gautier
"""

# This was written for Python 3, may be there is a chance it will work with Python 2...
from __future__ import print_function, unicode_literals

import io
import json
import time
try:
    import ijson.backends.yajl2 as ijson
except ImportError:
    import ijson


from . import DEFAULT_SEED, MinHash
from .logging import notify


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
        #assert prefix == prefix_item
        l.append(value)
        prefix, event, value = next(iterable)
    return tuple(l)


def _json_next_signature(iterable,
                         name = None,
                         filename = None,
                         ignore_md5sum=False,
                         prefix_item='abundances.item',
                         ijson = ijson):
    """Helper function to unpack and check one signature block only.
    - iterable: an iterable such the one returned by ijson.parse()
    - name:
    - filename:
    - ignore_md5sum:
    - prefix_item: required when parsing nested JSON structures
    - ijson: ijson backend to use.
    """
    from .signature import SourmashSignature

    d = dict()
    prefix, event, value = next(iterable)
    if event == 'start_map':
        prefix, event, value = next(iterable)
    while event != 'end_map':
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
    if n == 0xffffffff:               # load legacy signatures where n == -1
        n = 0
    max_hash = d.get('max_hash', 0)
    seed = d.get('seed', DEFAULT_SEED)

    molecule = d.get('molecule', 'DNA')
    if molecule == 'protein':
        is_protein = True
        dayhoff = False
        hp = False
    elif molecule == "dayhoff":
        is_protein = True
        dayhoff = True
        hp = False
    elif molecule == "hp":
        is_protein = True
        dayhoff = False
        hp = True
    elif molecule.upper() == 'DNA':
        is_protein = False
        dayhoff = False
        hp = False
    else:
        raise Exception("unknown molecule type: {}".format(molecule))


    track_abundance = False
    if 'abundances' in d:
        track_abundance = True

    e = MinHash(ksize=ksize, n=n, is_protein=is_protein,
                dayhoff=dayhoff, hp=hp,
                track_abundance=track_abundance,
                max_hash=max_hash, seed=seed)

    if not track_abundance:
        for m in mins:
            e.add_hash(m)
    else:
        abundances = list(map(int, d['abundances']))
        e.set_abundances(dict(zip(mins, abundances)))

    sig = SourmashSignature(e)

    if not ignore_md5sum:
        md5sum = d['md5sum']
        if md5sum != sig.md5sum():
            raise Exception('error loading - md5 of minhash does not match')

    if name:
        sig.d['name'] = name
    if filename:
        sig.d['filename'] = filename

    return sig

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

    # name, and filename not assumed to be parsed before the 'signatures'
    for sig in signatures:
        if 'name' in d:
            sig.d['name'] = d['name']
        if 'filename' in d:
            sig.d['filename'] = d['filename']

    # hardcode in support only for CC0 going forward
    if d.get('license', 'CC0') != 'CC0':
        raise Exception("sourmash only supports CC0-licensed signatures.")

    sig.d['license'] = d.get('license', 'CC0')

    return d


def load_signatureset_json_iter(data, ksize=None, ignore_md5sum=False, ijson=ijson):
    """
    - data: file handle (or file handle-like) object
    - ksize:
    - ignore_md5sum:
    - ijson: ijson backend
    """

    parser = ijson.parse(data)

    prefix, event, value = next(parser)
    assert prefix == '' and event == 'start_array' and value is None

    n = 0
    while True:
        try:
            sig = load_signature_json(parser,
                                      prefix_item = 'item.signatures.item.mins.item',
                                      ignore_md5sum=ignore_md5sum,
                                      ijson=ijson)
            if not ksize or ksize == sig.minhash.ksize:
                yield sig
        except ValueError:
            # possible end of the array of signatures
            try:
                prefix, event, value = next(parser)
                assert event == 'end_array'
            except StopIteration:
                pass
            finally:
                break
        n += 1

def load_signatures_json(data, ksize=None, ignore_md5sum=True, ijson=ijson):
    """
    - data: file handle (or file handle-like) object
    - ksize:
    - ignore_md5sum:
    - ijson: ijson backend
    """
    n = 0

    if isinstance(data, str):
        data = io.BytesIO(data.encode('utf-8'))

    it = load_signatureset_json_iter(data, ksize=ksize,
                                     ignore_md5sum=ignore_md5sum,
                                     ijson=ijson)

    for n, sigset in enumerate(it):
        if n > 0 and n % 100 == 0:
            notify('\r...sig loading {:,}', n, end='', flush=True)
        for sig in sigset['signatures']:
            yield sig

    if n > 1:
        notify('\r...sig loading {:,}', n, flush=True)


def add_meta_save(siglist):
    """ Convert one signature into a JSON dict
    - siglist: sequence of SourmashSignature objects
    - index: index of siglist to save
    """
    from .signature import SIGNATURE_VERSION
    records = []
    top_records = {}

    for sig in siglist:
        name, filename, sketch = sig._save()
        k = (name, filename)
        x = top_records.get(k, [])
        x.append(sketch)
        top_records[k] = x

    if not top_records:
        return records

    for (name, filename), sketches in top_records.items():
        record = {}
        if name:
            record['name'] = name
        if filename:
            record['filename'] = filename
        record['signatures'] = sketches

        record['version'] = SIGNATURE_VERSION
        record['class'] = 'sourmash_signature'
        record['hash_function'] = '0.murmur64'
        record['license'] = 'CC0'
        record['email'] = ''
        records.append(record)

    return records


def write_records_to_json(records, fp=None, indent=None, sort_keys=True):
    s = json.dumps(records, indent=indent, sort_keys=sort_keys, separators=(str(','), str(':')))
    if fp:
        try:
            fp.write(s)
        except TypeError:  # fp is opened in binary mode
            try:
                fp.write(s.encode('utf-8'))
            except TypeError:  # Python 2
                fp.write(unicode(s))
        return None
    return s


def save_signatures_json(
        siglist, fp=None, indent=None, sort_keys=True):
    """ Save multiple signatures into a JSON string (or into file handle 'fp')
    - siglist: sequence of SourmashSignature objects
    - fp: file handle to the location of a sig file
    - indent: indentation spaces (an integer) or if None no indentation
    - sort_keys: sort the keys in mappings before writting to JSON
    """
    startt = time.time()
    records = add_meta_save(siglist)
    if records == []:
        return ""
    s = write_records_to_json(records, fp, indent, sort_keys)
    notify("time taken to save signatures is {:.5f} seconds", time.time() - startt, end="\r")
    return s
