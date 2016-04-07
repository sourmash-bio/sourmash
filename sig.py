#! /usr/bin/env python
import yaml
import hashlib
import sys
from cStringIO import StringIO
import sourmash

class SourmashSignature(object):
    def __init__(self, email, estimator):
        self.d = {}
        self.d['class'] = 'sourmash_signature'
        self.d['type'] = 'mrnaseq'
        self.d['email'] = email
        self.d['version'] = '0.1'
        
        self.sketch = {}
        self.sketch['ksize'] = int(estimator._kh.ksize())
        self.sketch['prime'] = estimator.p
        self.sketch['mins'] = estimator._mins
        self.d['signature'] = self.sketch
        self.estimator = estimator

    def save(self):
        return yaml.dump(self.d)


def load_signature(data):
    d = yaml.safe_load(data)
    if d.get('class') != 'sourmash_signature':
        raise Exception("incorrect class: %s" % d.get('class'))
    email = d['email']

    sketch = d['signature']
    ksize = sketch['ksize']
    prime = sketch['prime']
    mins = map(long, sketch['mins'])
    e = sourmash.Estimators(ksize=ksize, max_prime=prime, n=len(mins))
    e._mins = mins

    return SourmashSignature(email, e)

def test_roundtrip():
    e = sourmash.Estimators()
    sig = SourmashSignature('titus@idyll.org', e)
    print(sig.save())
    s = sig.save()
    sig2 = load_signature(s)
    e2 = sig2.estimator
    
    assert e.jaccard(e2) == 1.0
