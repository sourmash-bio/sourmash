# Prepared databases

## GTDB R06-rs202 - DNA databases

All files below are available under https://osf.io/wxf9z/. The GTDB taxonomy spreadsheet (in a format suitable for `sourmash lca index`) is available [here](https://osf.io/p6z3w/).

For each k-mer size, three databases are available.

* Zipfile collections can be used for a linear search. The signatures were calculated with a scaled of 1000, which robustly supports searches for ~10kb or larger matches.
* SBT databases are indexed versions of the Zipfile collections that support faster search. They are also indexed with scaled=1000.
* LCA databases are indexed versions of the Zipfile collections that also contain taxonomy information and can be used with regular search as well as with [the `lca` subcommands for taxonomic analysis](https://sourmash.readthedocs.io/en/latest/command-line.html#sourmash-lca-subcommands-for-taxonomic-classification). They are indexed with scaled=10,000, which robustly supports searches for 100kb or larger matches.

You can read more about the different database and index types [here](https://sourmash.readthedocs.io/en/latest/command-line.html#indexed-databases).

Legacy databases are available [here](legacy-databases.md)

Note that the SBT and LCA databases can be used with sourmash v3.5 and later, while Zipfile collections can only be used with sourmash v4.1.0 and up.

### GTDB genomic representatives (47.8k genomes)

The GTDB genomic representatives are a low-redundancy subset of Genbank genomes.

| K-mer size | Zipfile collection | SBT | LCA |
| -------- | -------- | -------- | ---- |
| 21 | [download (1.3 GB)](https://osf.io/jp5zh/download) | [download (2.6 GB)](https://osf.io/py92w/download) | [download (114 MB)](https://osf.io/gk2za/download) | 
| 31 | [download (1.3 GB)](https://osf.io/nqmau/download) | [download (2.6 GB)](https://osf.io/w4bcm/download) | [download (131 MB)](https://osf.io/ypsjq/download) | 
| 51 | [download (1.3 GB)](https://osf.io/px6qd/download) | [download (2.6 GB)](https://osf.io/rv9zp/download) | [download (137 MB)](https://osf.io/297dp/download) | 

### GTDB all genomes (258k genomes)

These databases contain the complete GTDB collection of 258,406 genomes.

| K-mer size | Zipfile collection | SBT | LCA |
| -------- | -------- | -------- | ---- |
| 21 | [download (7.8 GB)](https://osf.io/vgex4/download) | [download (15 GB)](https://osf.io/ar67j/download) | [download (266 MB)](https://osf.io/hm3c4/download) | 
| 31 | [download (7.8 GB)](https://osf.io/94mzh/download) | [download (15 GB)](https://osf.io/dmsz8/download) | [download (286 MB)](https://osf.io/9xdg2/download) | 
| 51 | [download (7.8 GB)](https://osf.io/x9cdp/download) | [download (15 GB)](https://osf.io/8fc3t/download) | [download (299 MB)](https://osf.io/3cdp6/download)  | 
