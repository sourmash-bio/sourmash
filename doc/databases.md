# Prepared search databases

We provide several databases for download. Note that these databases can
be used with both sourmash v3.5 and sourmash v4.0.

## RefSeq microbial genomes - SBT

These database are formatted for use with `sourmash search` and
`sourmash gather`. They are calculated with a scaled value of 2000.

Approximately 91,000 microbial genomes (including viral and fungal)
from NCBI RefSeq.

* [RefSeq k=21, 2018.03.29][0] - 3.3 GB - [manifest](https://osf.io/wamfk/download)
* [RefSeq k=31, 2018.03.29][1] - 3.3 GB - [manifest](https://osf.io/x3aut/download)
* [RefSeq k=51, 2018.03.29][2] - 3.4 GB - [manifest](https://osf.io/zpkau/download)

## Genbank microbial genomes - SBT

These database are formatted for use with `sourmash search` and
`sourmash gather`.

Approximately 98,000 microbial genomes (including viral and fungal)
from NCBI Genbank.

* [Genbank k=21, 2018.03.29][3] - 3.9 GB - [manifest](https://osf.io/vm5kb/download)
* [Genbank k=31, 2018.03.29][4] - 3.9 GB - [manifest](https://osf.io/p87ec/download)
* [Genbank k=51, 2018.03.29][5] - 3.9 GB - [manifest](https://osf.io/cbxg9/download)

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

[0]: https://sourmash-databases.s3-us-west-2.amazonaws.com/zip/refseq-k21.sbt.zip
[1]: https://sourmash-databases.s3-us-west-2.amazonaws.com/zip/refseq-k31.sbt.zip
[2]: https://sourmash-databases.s3-us-west-2.amazonaws.com/zip/refseq-k51.sbt.zip

[3]: https://sourmash-databases.s3-us-west-2.amazonaws.com/zip/genbank-k21.sbt.zip
[4]: https://sourmash-databases.s3-us-west-2.amazonaws.com/zip/genbank-k31.sbt.zip
[5]: https://sourmash-databases.s3-us-west-2.amazonaws.com/zip/genbank-k51.sbt.zip

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
