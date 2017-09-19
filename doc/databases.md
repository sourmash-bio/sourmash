# Prepared search databases

## RefSeq microbial genomes

These database are formatted for use with `sourmash sbt_search` and
`sourmash sbt_gather`.

Approximately 60,000 microbial genomes (including viral and fungal)
from NCBI RefSeq.

* [RefSeq k=21, 2017.05.09](https://s3-us-west-1.amazonaws.com/spacegraphcats.ucdavis.edu/microbe-refseq-sbt-k21-2017.05.09.tar.gz) - 3.5 GB
* [RefSeq k=31, 2017.05.09](https://s3-us-west-1.amazonaws.com/spacegraphcats.ucdavis.edu/microbe-refseq-sbt-k31-2017.05.09.tar.gz) - 3.5 GB
* [RefSeq k=51, 2017.05.09](https://s3-us-west-1.amazonaws.com/spacegraphcats.ucdavis.edu/microbe-refseq-sbt-k51-2017.05.09.tar.gz) - 3.5 GB

## Genbank microbial genomes

These database are formatted for use with `sourmash sbt_search` and
`sourmash sbt_gather`.

Approximately 100,000 microbial genomes (including viral and fungal)
from NCBI Genbank.

* [Genbank k=21, 2017.05.09](https://s3-us-west-1.amazonaws.com/spacegraphcats.ucdavis.edu/microbe-genbank-sbt-k21-2017.05.09.tar.gz) - 4.2 GB
* [Genbank k=31, 2017.05.09](https://s3-us-west-1.amazonaws.com/spacegraphcats.ucdavis.edu/microbe-genbank-sbt-k31-2017.05.09.tar.gz) - 4.2 GB
* [Genbank k=51, 2017.05.09](https://s3-us-west-1.amazonaws.com/spacegraphcats.ucdavis.edu/microbe-genbank-sbt-k51-2017.05.09.tar.gz) - 4.2 GB

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
