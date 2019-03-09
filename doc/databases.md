# Prepared search databases

## RefSeq microbial genomes - SBT

These database are formatted for use with `sourmash search` and
`sourmash gather`. They are calculated with a scaled value of 2000.

Approximately 60,000 microbial genomes (including viral and fungal)
from NCBI RefSeq.

* [RefSeq k=21, 2018.03.29][0] - 7 GB
* [RefSeq k=31, 2018.03.29][1] - 7 GB
* [RefSeq k=51, 2018.03.29][2] - 7.1 GB

## Genbank microbial genomes - SBT

These database are formatted for use with `sourmash search` and
`sourmash gather`.

Approximately 100,000 microbial genomes (including viral and fungal)
from NCBI Genbank.

* [Genbank k=21, 2018.03.29][3]- 8.4 GB
* [Genbank k=31, 2018.03.29][4] - 8.4 GB
* [Genbank k=51, 2018.03.29][5] - 8.4 GB

### Details

The individual signatures for the above SBTs were calculated as follows:

```
sourmash compute -k 4,5 \
                 -n 2000 \
                 --track-abundance \
                 --name-from-first \
                 -o {output} \
                 {input}

sourmash compute -k 21,31,51 \
                 --scaled 2000 \
                 --track-abundance \
                 --name-from-first \
                 -o {output} \
                 {input}
```

See [github.com/dib-lab/sourmash_databases](https://github.com/dib-lab/sourmash_databases) for a Snakemake workflow
to build the databases.

[0]: https://s3-us-west-2.amazonaws.com/sourmash-databases/2018-03-29/refseq-d2-k21.tar.gz
[1]: https://s3-us-west-2.amazonaws.com/sourmash-databases/2018-03-29/refseq-d2-k31.tar.gz
[2]: https://s3-us-west-2.amazonaws.com/sourmash-databases/2018-03-29/refseq-d2-k51.tar.gz
[3]: https://s3-us-west-2.amazonaws.com/sourmash-databases/2018-03-29/genbank-d2-k21.tar.gz
[4]: https://s3-us-west-2.amazonaws.com/sourmash-databases/2018-03-29/genbank-d2-k31.tar.gz
[5]: https://s3-us-west-2.amazonaws.com/sourmash-databases/2018-03-29/genbank-d2-k51.tar.gz

## Genbank LCA Database

These databases are formatted for use with `sourmash lca`; they are
v2 LCA databases and will work with sourmash v2.0a11 and later.
They are calculated with a scaled value of 10000 (1e5).

Approximately 87,000 microbial genomes (including viral and fungal)
from NCBI Genbank.

* [Genbank k=21, 2017.11.07](https://osf.io/d7rv8/download), 109 MB
* [Genbank k=31, 2017.11.07](https://osf.io/4f8n3/download), 120 MB
* [Genbank k=51, 2017.11.07](https://osf.io/nemkw/download), 125 MB

### Details

The above LCA databases were calculated as follows:

```
sourmash lca index genbank-genomes-taxonomy.2017.05.29.csv \
    genbank-k21.lca.json.gz -k 21 --scaled=10000 \
    -f --traverse-directory .sbt.genbank-k21 --split-identifiers
```

See
[github.com/dib-lab/2018-ncbi-lineages](https://github.com/dib-lab/2018-ncbi-lineages)
for information on preparing the genbank-genomes-taxonomy file.
