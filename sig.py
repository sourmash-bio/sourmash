#! /usr/bin/env python
import yaml
import hashlib
import sys
import sourmash
import hashlib

class SourmashSignature(object):
    def __init__(self, email, estimator, name='', filename=''):
        self.d = {}
        self.d['class'] = 'sourmash_signature'
        self.d['type'] = 'mrnaseq'
        self.d['email'] = email
        self.d['version'] = '0.1'
        if name:
            self.d['name'] = name
        if filename:
            self.d['filename'] = filename
        
        self.estimator = estimator

    def md5sum(self):
        m = hashlib.md5()
        for k in self.estimator._mins:
            m.update(str(k))
        return m.hexdigest()

    def name(self):
        if 'name' in self.d:
            return self.d.get('name')
        elif 'filename' in self.d:
            return self.d.get('filename')
        else:
            return self.md5sum()[:8]

    def save(self):
        e = dict(self.d)
        estimator = self.estimator
        
        sketch = {}
        sketch['ksize'] = int(estimator._kh.ksize())
        sketch['prime'] = estimator.p
        sketch['mins'] = estimator._mins
        sketch['md5sum'] = self.md5sum()
        e['signature'] = sketch
                           
        return yaml.dump(e)

    def jaccard(self, other):
        return self.estimator.jaccard(other.estimator)


def load_signature(data, ignore_md5sum=False):
    d = yaml.safe_load(data)
    if d.get('class') != 'sourmash_signature':
        raise Exception("incorrect class: %s" % d.get('class'))
    email = d['email']

    sketch = d['signature']
    ksize = sketch['ksize']
    prime = sketch['prime']
    mins = list(map(int, sketch['mins']))
    e = sourmash.Estimators(ksize=ksize, max_prime=prime, n=len(mins))
    e._mins = mins

    sig = SourmashSignature(email, e)
    if not ignore_md5sum:
        md5sum = sketch['md5sum']
        if md5sum != sig.md5sum():
            raise Exception('error loading - md5 of estimator does not match')

    keys = list(d.keys())
    keys.remove('email')
    keys.remove('signature')
    for k in keys:
        sig.d[k] = d[k]

    return sig

def test_roundtrip():
    e = sourmash.Estimators()
    sig = SourmashSignature('titus@idyll.org', e)
    s = sig.save()
    sig2 = load_signature(s)
    e2 = sig2.estimator
    
    assert e.jaccard(e2) == 1.0

def test_md5():
    e = sourmash.Estimators()
    sig = SourmashSignature('titus@idyll.org', e)
    print(sig.save())
    assert sig.md5sum() == '9e1f3d700a7344f61db1d37fa11f60f2', sig.md5sum()
