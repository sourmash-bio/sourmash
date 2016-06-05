====================================
Using sourmash from the command line
====================================

An example
==========

Grab three bacterial genomes from NCBI::

   curl -L -O ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Escherichia_coli/reference/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
   curl -L -O ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Salmonella_enterica/reference/GCF_000006945.1_ASM694v1/GCF_000006945.1_ASM694v1_genomic.fna.gz
   curl -L -O ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Sphingobacteriaceae_bacterium_DW12/latest_assembly_versions/GCF_000783305.1_ASM78330v1/GCF_000783305.1_ASM78330v1_genomic.fna.gz

Compute signatures for each::

   sourmash compute -f *.fna.gz

This will produce three `.sig` files containing MinHash signatures at k=31;
the `-f` bypasses an error where the last of the genomes has some non-ATCGN
characters in it.

Next, compare all the signatures to each other::

   sourmash compare *.sig -o cmp

Finally, plot a dendrogram::

   ./plot-comparison.py cmp --pdf

Mention:

* reads fa, fq, gz, bz2
* show png
  
