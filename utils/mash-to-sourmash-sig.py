#! /usr/bin/env python
import sys
import json
import sourmash_lib, sourmash_lib.signature

x = json.loads(open(sys.argv[1]).read())
ksize = x['kmer']
num = x['sketchSize']

assert x['hashType'] == "MurmurHash3_x64_128"
assert x['hashBits'] == 64
assert x['hashSeed'] == 42

xx = x['sketches'][0]
hashes = xx['hashes']

# assert stuff HERE


E = sourmash_lib.Estimators(ksize=ksize, n=num, is_protein=False)
mh = E.mh
for h in hashes:
    mh.add_hash(h)

s = sourmash_lib.signature.SourmashSignature('', E, name=sys.argv[1])
print(sourmash_lib.signature.save_signatures([s]))

