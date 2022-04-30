# Prepared databases

```{contents}
```

We provide a number of pre-built collections and indexed databases
that you can use with sourmash.

## Types of databases

For each k-mer size, three types of databases may be available: Zipfile (`.zip`), SBT (`.sbt.zip`), and LCA (`.lca.jzon.gz`).
We recommend using the Zipfile databases for `sourmash gather` and the SBT databases for `sourmash search`. You must use the LCA databases for `sourmash lca` operations.

You can read more about the different database and index types [here](https://sourmash.readthedocs.io/en/latest/command-line.html#indexed-databases).

Note that the SBT and LCA databases can be used with sourmash v3.5 and later, while Zipfile collections can only be used with sourmash v4.1 and up.

## Downloading and using the databases

All databases below can be downloaded via the command line with `curl -L <url> -o <output>`, where `<url>` is the URL below, and `<output>` is the filename you want to use locally.

The databases do not need to be unpacked or prepared in any way after download.

You can verify that they've been successfully downloaded with `sourmash sig summarize <output>`.

## GTDB R07-RS207 - DNA databases

[GTDB R07-RS207](https://forum.gtdb.ecogenomic.org/t/announcing-gtdb-r07-rs207/264) consists of 317,542 genomes organized into 65,703 species clusters.

The lineage spreadsheet (for `sourmash tax` commands) is available [at the species level](https://osf.io/v3zmg/download) and [at the strain level](https://osf.io/r87td/download).

### GTDB R07-RS207 genomic representatives (66k)

The GTDB genomic representatives are a low-redundancy subset of Genbank genomes, with 65,703 species-level genomes.

| K-mer size | Zipfile collection | SBT | LCA |
| -------- | -------- | -------- | ---- |
| 21 | [download (1.7 GB)](https://osf.io/f2wzc/download) | [download (3.5 GB)](https://osf.io/zsypg/download) | [download (181 MB)](https://osf.io/pm35d/download) |
| 31 | [download (1.7 GB)](https://osf.io/3a6gn/download) | [download (3.5 GB)](https://osf.io/ernct/download) | [download (181 MB)](https://osf.io/p9ezm/download) |
| 51 | [download (1.7 GB)](https://osf.io/f23qn/download) | [download (3.5 GB)](https://osf.io/yq7dc/download) | [download (181 MB)](https://osf.io/8qhgy/download) |

### GTDB R07-RS207 all genomes (318k)

These are databases for the full GTDB release, each containing 317,542 genomes.

| K-mer size | Zipfile collection | SBT | LCA |
| -------- | -------- | -------- | ---- |
| 21 | [download (9.4 GB)](https://osf.io/9gpck/download) | [download (19 GB)](https://osf.io/wr8pk/download) | [download (351 MB)](https://osf.io/su9za/download) |
| 31 | [download (9.4 GB)](https://osf.io/k2u8s/download) | [download (19 GB)](https://osf.io/748ew/download) | [download (351 MB)](https://osf.io/tf3ah/download) |
| 51 | [download (9.4 GB)](https://osf.io/ubt7p/download) | [download (19 GB)](https://osf.io/78hdr/download) | [download (351 MB)](https://osf.io/vc8ua/download) |

## Genbank genomes from March 2022

The below zip files contain signatures for all microbial Genbank genomes as of March 2022, based on the assembly_summary files provided [here](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/).

Since some of the files are extremely large, we only provide them in Zip format.

Taxonomic spreadsheets for each domain are provided below as well.

### Genbank viral

47,952 genomes:

[genbank-2022.03-viral-k21.zip](https://dweb.link/ipfs/bafybeicjyx6qkhdtw6q4cxs6fyl46gqfhd4q5eqje5lkswf2npljnyytzi)

[genbank-2022.03-viral-k31.zip](https://dweb.link/ipfs/bafybeibqsldwsztjf66rwvwnb6hamjtsfkmdk5bmfqbzwrod6wwwkqz2ya)

[genbank-2022.03-viral-k51.zip](https://dweb.link/ipfs/bafybeibgifuv4q3mihfubnhjhwm2esjnoseudpnzwahlkp3hvlbmtd4s2q)

[genbank-2022.03-viral.lineages.csv.gz](https://osf.io/j4tsu/download)

### Genbank archaeal

8,750 genomes:

[genbank-2022.03-archaea-k21.zip](https://dweb.link/ipfs/bafybeiepywe7c6zjzgh3rksqiwpo5zyb7uuefbgvbn5nkgiq77iaavpzl4)

[genbank-2022.03-archaea-k31.zip](https://dweb.link/ipfs/bafybeidn6epju7yrdxrktq5wjko2yiwp6nrx3mq37htiuwecm7lffrbcdi)

[genbank-2022.03-archaea-k51.zip](https://dweb.link/ipfs/bafybeifyrwbx5dnay4mflboc5zai2de3xrvcxtgiu4j7adzj6qrxhb3zva)

[genbank-2022.03-archaea.lineages.csv.gz](https://osf.io/kcbpn/download)


### Genbank protozoa

1193 genomes:

[genbank-2022.03-protozoa-k21.zip](https://dweb.link/ipfs/bafybeicfh4xl4wuxd4xy2tf73hfamxlqa3s2higghnsjay4t5wtlmsdo5y)

[genbank-2022.03-protozoa-k31.zip](https://dweb.link/ipfs/bafybeicpxjhfrzem7f34eghbbwm3vglz2njxo72vpqcw7foilfomexsghi)

[genbank-2022.03-protozoa-k51.zip](https://dweb.link/ipfs/bafybeigfpxkmzyq6sdkob53l6ztiy5ro44dzkad7dxakhuaao6cw4gp4eu)

[genbank-2022.03-protozoa.lineages.csv.gz](https://osf.io/2x8u4/download)


### Genbank fungi

10,286 genomes:

[genbank-2022.03-fungi-k21.zip](https://dweb.link/ipfs/bafybeibrirvek4lxn6hh3wgsmtsd5vz5gtmewpzeg364bix3hojghwmygq)

[genbank-2022.03-fungi-k31.zip](https://dweb.link/ipfs/bafybeidhhwvwujkteno5ugwgjy4brhrv5dff2aumifcuew73qolfktdndq)

[genbank-2022.03-fungi-k51.zip](https://dweb.link/ipfs/bafybeibnrtt45f7wez2xb3fy5rxhatpeevc3rilm2gs65u5h6gc4u72fam)

[genbank-2022.03-fungi.lineages.csv.gz](https://osf.io/s4b85/download)

### Genbank bacterial:

1,148,011 genomes:

[genbank-2022.03-bacteria-k21.zip](https://dweb.link/ipfs/bafybeif2hdztfrevkngnfqk3bsoyajxxf67o57u4dezbz647jwcf6gnwoy)

[genbank-2022.03-bacteria-k31.zip](https://dweb.link/ipfs/bafybeigkcvizvhe3xzxsuzv3ryf3ogvgvcmms2e5nfk7epl5egts22jyue)

[genbank-2022.03-bacteria-k51.zip](https://dweb.link/ipfs/bafybeie3eyyectnh5xqxz44oa3qj5vura3bffqdwfqk6jjuzzadkh7e2sq)

[genbank-2022.03-bacteria.lineages.csv.gz](https://osf.io/4agsp/download)

## GTDB R06-RS202 - DNA databases

All files below are available under https://osf.io/wxf9z/. The GTDB taxonomy spreadsheet (in a format suitable for `sourmash lca index`) is available [here](https://osf.io/p6z3w/).

### GTDB R06-RS202 genomic representatives (47.8k)

The GTDB genomic representatives are a low-redundancy subset of Genbank genomes.

| K-mer size | Zipfile collection | SBT | LCA |
| -------- | -------- | -------- | ---- |
| 21 | [download (1.3 GB)](https://osf.io/jp5zh/download) | [download (2.6 GB)](https://osf.io/py92w/download) | [download (114 MB)](https://osf.io/gk2za/download) | 
| 31 | [download (1.3 GB)](https://osf.io/nqmau/download) | [download (2.6 GB)](https://osf.io/w4bcm/download) | [download (131 MB)](https://osf.io/ypsjq/download) | 
| 51 | [download (1.3 GB)](https://osf.io/px6qd/download) | [download (2.6 GB)](https://osf.io/rv9zp/download) | [download (137 MB)](https://osf.io/297dp/download) | 

### GTDB all genomes (258k)

These databases contain the complete GTDB collection of 258,406 genomes.

| K-mer size | Zipfile collection | SBT | LCA |
| -------- | -------- | -------- | ---- |
| 21 | [download (7.8 GB)](https://osf.io/vgex4/download) | [download (15 GB)](https://osf.io/ar67j/download) | [download (266 MB)](https://osf.io/hm3c4/download) | 
| 31 | [download (7.8 GB)](https://osf.io/94mzh/download) | [download (15 GB)](https://osf.io/dmsz8/download) | [download (286 MB)](https://osf.io/9xdg2/download) | 
| 51 | [download (7.8 GB)](https://osf.io/x9cdp/download) | [download (15 GB)](https://osf.io/8fc3t/download) | [download (299 MB)](https://osf.io/3cdp6/download)  | 

## Appendix: database use and construction details

Database release workflows are being archived at [sourmash-bio/database-releases](https://github.com/sourmash-bio/database-releases).

Some more details on database use and construction:

* Zipfile collections can be used for a linear search. The signatures were calculated with a scaled of 1000, which robustly supports searches for ~10kb or larger matches.
* SBT databases are indexed versions of the Zipfile collections that support faster search. They are also indexed with scaled=1000.
* LCA databases are indexed versions of the Zipfile collections that also contain taxonomy information and can be used with regular search as well as with [the `lca` subcommands for taxonomic analysis](https://sourmash.readthedocs.io/en/latest/command-line.html#sourmash-lca-subcommands-for-taxonomic-classification). They are indexed with scaled=10,000, which robustly supports searches for 100kb or larger matches.

## Appendix: Memory and time requirements

The detailed memory usage of sourmash depends on the type of search, the query, and the database you're searching, but to help guide you here is a range of numbers:

| Search type | Query | Database | Max RAM | Time |
| -------- | -------- | -------- | -------- | -------- |
| gather | Bacterial genome | GTDB complete (280k) | 1 GB | 6 minutes |
| gather | Simple metagenome | GTDB reps .zip (65k)     |   2 GB   | 6 minutes |
| gather | Real metagenome | All Genbank (1.2m) | 100 GB | 3 hours
| lca summarize     |Simple metagenome | GTDB reps .sql (65k)     |   400 MB   | 20 seconds |
| lca summarize | Simple metagenome | GTDB reps .json (65k) | 6.2 GB | 1m 20 seconds |


Please see [sourmash#1958](https://github.com/sourmash-bio/sourmash/issues/1958) for detailed GTDB numbers and [gather paper#47](https://github.com/dib-lab/2020-paper-sourmash-gather/issues/47) for detailed Genbank numbers.

## Appendix: legacy databases

Legacy databases are available [here](legacy-databases.md).
