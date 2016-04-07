#! /usr/bin/env python
import yaml
import hashlib
import sys
import sourmash

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
        
        self.sketch = {}
        self.sketch['ksize'] = int(estimator._kh.ksize())
        self.sketch['prime'] = estimator.p
        self.sketch['mins'] = estimator._mins
        self.d['signature'] = self.sketch
        self.estimator = estimator

    def name(self):
        if 'name' in self.d:
           return self.d.get('name')
        elif 'filename' in self.d:
           return self.d.get('filename')
        else:
           return ""

    def save(self):
        return yaml.dump(self.d)

    def jaccard(self, other):
        return self.estimator.jaccard(other.estimator)


def load_signature(data):
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

    keys = list(d.keys())
    keys.remove('email')
    keys.remove('signature')
    for k in keys:
        sig.d[k] = d[k]

    return sig

def test_roundtrip():
    e = sourmash.Estimators()
    sig = SourmashSignature('titus@idyll.org', e)
    print(sig.save())
    s = sig.save()
    sig2 = load_signature(s)
    e2 = sig2.estimator
    
    assert e.jaccard(e2) == 1.0
