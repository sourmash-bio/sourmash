#! /usr/bin/env python
"""
Save and load MinHash sketches in a JSON format, along with some metadata.
"""
from __future__ import print_function
import hashlib

import gzip
import bz2file
import io
import sys

from . import signature_json
from .logging import error


SIGNATURE_VERSION=0.4


class SourmashSignature(object):
    "Main class for signature information."

    def __init__(self, minhash, name='', filename=''):
        self.d = {}
        self.d['class'] = 'sourmash_signature'
        if name:
            self.d['name'] = name
        if filename:
            self.d['filename'] = filename

        self.minhash = minhash
        self.d['license'] = 'CC0'

    def __hash__(self):
        return hash(self.md5sum())

    def __str__(self):
        name = self.name()
        md5pref = self.md5sum()[:8]
        if name != md5pref:
            return "SourmashSignature('{}', {})".format(name, md5pref)
        return "SourmashSignature({})".format(md5pref)
    __repr__ = __str__

    def md5sum(self):
        "Calculate md5 hash of the bottom sketch, specifically."
        m = hashlib.md5()
        m.update(str(self.minhash.ksize).encode('ascii'))
        for k in self.minhash.get_mins():
            m.update(str(k).encode('utf-8'))
        return m.hexdigest()

    def __eq__(self, other):
        allkeys = set(self.d.keys()).union(set(other.d.keys()))
        for k in allkeys:
            if self.d.get(k) != other.d.get(k):
                return False

        return self.minhash == other.minhash

    def __ne__(self, other):
        return not self == other

    def name(self):
        "Return as nice a name as possible, defaulting to md5 prefix."
        if 'name' in self.d:
            return self.d.get('name')
        elif 'filename' in self.d:
            return self.d.get('filename')
        else:
            return self.md5sum()[:8]

    def _display_name(self, max_length):
        if 'name' in self.d:
            name = self.d['name']
            if len(name) > max_length:
                name = name[:max_length - 3] + '...'
        elif 'filename' in self.d:
            name = self.d['filename']
            if len(name) > max_length:
                name = '...' + name[-max_length + 3:]
        else:
            name = self.md5sum()[:8]
        assert len(name) <= max_length
        return name

    def _save(self):
        "Return metadata and a dictionary containing the sketch info."
        e = dict(self.d)
        minhash = self.minhash

        sketch = {}
        sketch['ksize'] = int(minhash.ksize)
        sketch['num'] = minhash.num
        sketch['max_hash'] = minhash.max_hash
        sketch['seed'] = int(minhash.seed)
        if self.minhash.track_abundance:
            values = minhash.get_mins(with_abundance=True)
            sketch['mins'] = list(map(int, values.keys()))
            sketch['abundances'] = list(map(int, values.values()))
        else:
            sketch['mins'] = list(map(int, minhash.get_mins()))
        sketch['md5sum'] = self.md5sum()

        if minhash.is_protein and not minhash.dayhoff and not minhash.hp:
            sketch['molecule'] = 'protein'
        elif minhash.dayhoff:
            sketch['molecule'] = 'dayhoff'
        elif minhash.hp:
            sketch['molecule'] = 'hp'
        else:
            sketch['molecule'] = 'DNA'

        e['signature'] = sketch

        return self.d.get('name'), self.d.get('filename'), sketch

    def similarity(self, other, ignore_abundance=False, downsample=False):
        "Compute similarity with the other MinHash signature."
        try:
            return self.minhash.similarity(other.minhash, ignore_abundance)
        except ValueError as e:
            if 'mismatch in max_hash' in str(e) and downsample:
                xx = self.minhash.downsample_max_hash(other.minhash)
                yy = other.minhash.downsample_max_hash(self.minhash)
                return xx.similarity(yy, ignore_abundance)
            else:
                raise

    def jaccard(self, other):
        "Compute Jaccard similarity with the other MinHash signature."
        return self.minhash.similarity(other.minhash, True)

    def contained_by(self, other, downsample=False):
        "Compute containment by the other signature. Note: ignores abundance."
        try:
            return self.minhash.contained_by(other.minhash)
        except ValueError as e:
            if 'mismatch in max_hash' in str(e) and downsample:
                xx = self.minhash.downsample_max_hash(other.minhash)
                yy = other.minhash.downsample_max_hash(self.minhash)
                return xx.contained_by(yy)
            else:
                raise


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


def load_signatures(data, ksize=None, select_moltype=None,
                    ignore_md5sum=False, do_raise=False, quiet=False):
    """Load a JSON string with signatures into classes.

    Returns list of SourmashSignature objects.

    Note, the order is not necessarily the same as what is in the source file.
    """
    if ksize:
        ksize = int(ksize)

    if not data:
        return

    is_fp = False
    if hasattr(data, 'find') and data.find('sourmash_signature') == -1:   # filename
        done = False
        try:                                  # is it a file handle?
            data.read
            is_fp = True
            done = True
        except AttributeError:
            pass

        # not a file handle - treat it like a filename.
        if not done:
            try:
                data = _guess_open(data)
                is_fp = True
            except OSError as excinfo:
                if not quiet: error(str(excinfo))
                if do_raise:
                    raise
                return
    else:  # file-like
        if hasattr(data, 'mode'):  # file handler
            if 't' in data.mode:  # need to reopen handler as binary
                if sys.version_info >= (3, ):
                    data = data.buffer

    try:
        # JSON format
        for sig in signature_json.load_signatures_json(data,
                                                     ignore_md5sum=ignore_md5sum):
            if not ksize or ksize == sig.minhash.ksize:
                if not select_moltype or \
                     sig.minhash.is_molecule_type(select_moltype):
                    if select_moltype == 'protein':
                        if any(sig.minhash.is_molecule_type(t) for t in ('dayhoff', 'hp')):
                            # dayhoff and hp are also protein MHs. only yield
                            # sig if it is exactly one of (protein, hp, dayhoff)
                            continue
                    yield sig
    except Exception as e:
        if not quiet:
            error("Error in parsing signature; quitting.")
            error("Exception: {}", str(e))
        if do_raise:
            raise
    finally:
        if is_fp:
            data.close()


def load_one_signature(data, ksize=None, select_moltype=None,
                       ignore_md5sum=False):
    sigiter = load_signatures(data, ksize=ksize,
                              select_moltype=select_moltype,
                              ignore_md5sum=ignore_md5sum)

    try:
        first_sig = next(sigiter)
    except StopIteration:
        raise ValueError("no signatures to load")

    try:
        next(sigiter)
    except StopIteration:
        return first_sig

    raise ValueError("expected to load exactly one signature")


def save_signatures(siglist, fp=None):
    "Save multiple signatures into a JSON string (or into file handle 'fp')"
    return signature_json.save_signatures_json(siglist, fp)
