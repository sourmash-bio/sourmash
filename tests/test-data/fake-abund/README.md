signature:
```
sourmash compute -k 31 --track-abundance ../podar-ref/1.fa ../podar-ref/1.fa ../podar-ref/1.fa ../podar-ref/1.fa ../podar-ref/1.fa ../podar-ref/63.fa --merge "query" -o query.sig --scaled=1000
```

database:

```
sourmash lca index ../podar-lineage.csv matches.lca.json ../podar-ref/{1,63}.fa.sig -k 31 -C 3 --split-identifiers
```
