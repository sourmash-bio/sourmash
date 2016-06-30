# Misc utility scripts

## Misc scripts

* trim-noV.sh - a script to do trimming of short reads. requires khmer >= 2.0.
* setname.py - a script to set the 'name' in .sig files.

# Bulk download SRA scripts

Files for bulk downloading of echinoderm (sea urchin & friends) RNA
sequences from the Sequence Read Archive/ENA:

```
name-urchin.py
select-urchin.py
slurp_sra.py
```

## Instructions

The script `slurp_sra.py` will take a file like this:

```
"Experiment Accession","Experiment Title","Organism Name","Instrument","Submitter","Study Accession","Study Title","Sample Accession","Sample Title","Total Size, Mb","Total RUNs","Total Spots","Total Bases","Library Name","Library Strategy","Library Source","Library Selection"
"SRX1625120","RNA-Seq of  Ophiolimna perfida: field-collected adult body","Ophiolimna perfida","Illumina HiSeq 2000","Museum Victoria","SRP071599","Transcriptome-based phylogeny of the echinoderm class Ophiuroidea","SRS1334413","","1778.76","1","14928719","2985743800","MVF188866","RNA-Seq","TRANSCRIPTOMIC","RANDOM"
"SRX1625119","RNA-Seq of  Ophiocoma wendtii: field-collected adult body","Ophiocoma wendtii","Illumina HiSeq 2000","Museum Victoria","SRP071599","Transcriptome-based phylogeny of the echinoderm class Ophiuroidea","SRS1334414","","1940.88","1","16000000","3200000000","MVF193471","RNA-Seq","TRANSCRIPTOMIC","RANDOM"
"SRX1625118","RNA-Seq of  Ophioleuce brevispinum: field-collected adult body","Ophioleuce brevispinum","Illumina HiSeq 2000","Museum Victoria","SRP071599","Transcriptome-based phylogeny of the echinoderm class Ophiuroidea","SRS1334415","","1706.99","1","14372240","2874448000","MVF188879","RNA-Seq","TRANSCRIPTOMIC","RANDOM"
```

that contains a list of SRA records, and produce a file `ftp_list.csv` that looks like this:

```
SRX1625117,SRR3217922,ftp.sra.ebi.ac.uk/vol1/fastq/SRR321/002/SRR3217922/SRR3217922_1.fastq.gz,d9375ad599dbcc24dc29570ace7c328a,1167260213
SRX1625117,SRR3217922,ftp.sra.ebi.ac.uk/vol1/fastq/SRR321/002/SRR3217922/SRR3217922_2.fastq.gz,0c41ce2f0d7e80257ed45a91bc0c5a69,1172062623
SRX1625116,SRR3217921,ftp.sra.ebi.ac.uk/vol1/fastq/SRR321/001/SRR3217921/SRR3217921_1.fastq.gz,afa3f0c4763dfbd43fc6137c691fa927,1672839396
```

These URLs (third column) can be grabbed directly with curl or
wget. You generally want to take only URLs that have _1.fastq.gz in
them - _2 is the other end of fragments in _1 and hence correlated,
and no _1 or _2 is older-style sequences that are shorter and probably
less useful.

The way you get the first sra_result.csv file is by searching the SRA like so,

```
https://www.ncbi.nlm.nih.gov/sra/?term=txid7586%5BOrganism%3Aexp%5D+illumina
```

and then doing 'send to' (upper right) 'File'. There's probably a way
to do this programmatically but this works.

CTB 6/2016

