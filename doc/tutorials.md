# A sourmash tutorial

This tutorial should run without modification on Ubuntu 15.10 ("wily");
it was developed on AWS ami-05384865 (us-west region).

You'll need about 30 GB of free disk space to download the database,
and about 1-2 GB of RAM to search it.  The tutorial should take about
20 minutes total to run.

## Installing sourmash

To install sourmash, run:

```
sudo apt-get -y update && \
sudo apt-get install -y python3.5-dev python3.5-venv make \
    libc6-dev g++ zlib1g-dev
```

this installs Python 3.5.

Now, create a local software install and populate it with Jupyter and
other dependencies:

```

python3.5 -m venv ~/py3
. ~/py3/bin/activate
pip install -U pip
pip install -U Cython
pip install -U jupyter jupyter_client ipython pandas matplotlib scipy scikit-learn khmer

pip install -U https://github.com/dib-lab/sourmash/archive/master.zip

```

## Generate a signature for Illumina reads

Download some reads and a reference genome:

```
mkdir ~/data
cd ~/data
wget https://s3.amazonaws.com/public.ged.msu.edu/ecoli_ref-5m.fastq.gz
wget https://s3.amazonaws.com/public.ged.msu.edu/ecoliMG1655.fa.gz
```

Compute a scaled MinHash signature from our reads:

```
mkdir ~/sourmash
cd ~/sourmash

sourmash compute --scaled 10000 ~/data/ecoli_ref*pe*.fq.gz -o ecoli-reads.sig -k 31
```

## Compare reads to assemblies

Use case: how much of the read content is contained in the reference genome?

Build a signature for an E. coli genome:

```
sourmash compute --scaled 10000 -k 31 ~/data/ecoliMG1655.fa.gz -o ecoli-genome.sig
```

and now evaluate *containment*, that is, what fraction of the read content is
contained in the genome:

```
sourmash search -k 31 ecoli-reads.sig ecoli-genome.sig --containment
```

and you should see:

```
# running sourmash subcommand: search
loaded query: /home/ubuntu/data/ecoli_ref-5m... (k=31, DNA)
loaded 1 signatures from ecoli-genome.sig
1 matches:
similarity   match
----------   -----
 46.6%       /home/ubuntu/data/ecoliMG1655.fa.gz
```


Try the reverse - why is it bigger?

```
sourmash search -k 31 ecoli-genome.sig ecoli-reads.sig --containment
```

## Make and search a database quickly.

Suppose that we have a collection of signatures (made with `sourmash
compute` as above) and we want to search it with our newly assembled
genome (or the reads, even!). How would we do that?

Let's grab a sample collection of 50 E. coli genomes and unpack it --

```
mkdir ecoli_many_sigs
cd ecoli_many_sigs

curl -O -L https://github.com/dib-lab/sourmash/raw/master/data/eschericia-sigs.tar.gz

tar xzf eschericia-sigs.tar.gz
rm eschericia-sigs.tar.gz

cd ../

```

This will produce 50 files named `ecoli-N.sig` in the `ecoli_many_sigs` --

```
ls ecoli_many_sigs
```

Let's turn this into an easily-searchable database with `sourmash index` --

```
sourmash index -k 31 ecolidb ecoli_many_sigs/*.sig
```

and now we can search!

```
sourmash search ecoli-genome.sig ecolidb.sbt.json -n 20
```

You should see output like this:

```
select query k=31 automatically.
loaded query: /home/ubuntu/data/ecoliMG1655.... (k=31, DNA)
loaded 0 signatures and 1 databases total.                                     

49 matches; showing first 20:
similarity   match
----------   -----
 75.4%       NZ_JMGW01000001.1 Escherichia coli 1-176-05_S4_C2 e117605...
 72.2%       NZ_GG774190.1 Escherichia coli MS 196-1 Scfld2538, whole ...
 71.4%       NZ_JMGU01000001.1 Escherichia coli 2-011-08_S3_C2 e201108...
 70.1%       NZ_JHRU01000001.1 Escherichia coli strain 100854 100854_1...
 69.0%       NZ_JH659569.1 Escherichia coli M919 supercont2.1, whole g...
 64.9%       NZ_JNLZ01000001.1 Escherichia coli 3-105-05_S1_C1 e310505...
 63.0%       NZ_MOJK01000001.1 Escherichia coli strain 469 Cleandata-B...
 62.9%       NZ_MOGK01000001.1 Escherichia coli strain 676 BN4_676_1_(...
 62.0%       NZ_JHDG01000001.1 Escherichia coli 1-176-05_S3_C1 e117605...
 59.9%       NZ_MIWF01000001.1 Escherichia coli strain AF7759-1 contig...
 52.7%       NZ_KE700241.1 Escherichia coli HVH 147 (4-5893887) acYxy-...
 51.7%       NZ_APWY01000001.1 Escherichia coli 178200 gec178200.conti...
 49.3%       NZ_LVOV01000001.1 Escherichia coli strain swine72 swine72...
 49.3%       NZ_MIWP01000001.1 Escherichia coli strain K6412 contig_00...
 49.0%       NZ_LQWB01000001.1 Escherichia coli strain GN03624 GCID_EC...
 48.9%       NZ_JHGJ01000001.1 Escherichia coli O45:H2 str. 2009C-4780...
 48.1%       NZ_CP011331.1 Escherichia coli O104:H4 str. C227-11, comp...
 47.7%       NZ_JHNB01000001.1 Escherichia coli O103:H25 str. 2010C-45...
 47.7%       NZ_JHRE01000001.1 Escherichia coli strain 302014 302014_1...
 47.6%       NZ_JHHE01000001.1 Escherichia coli O103:H2 str. 2009C-327...
```

## Compare many signatures and build a tree.

Compare all the things:

```
sourmash compare ecoli_many_sigs/* -o ecoli_cmp
```

and then plot:

```
sourmash plot --pdf --labels ecoli_cmp
```

which will produce a file `ecoli_cmp.matrix.pdf` and `ecoli_cmp.dendro.pdf`
which you can then download via your file browser and view on your local
computer.

Here's a PNG version:

![E. coli comparison plot](_static/ecoli_cmp.matrix.png)

## What's in my metagenome?

Download and unpack a newer version of the k=31 RefSeq index described
in
[CTB's blog post][0]
-- this one contains sketches of all 100k Genbank microbes. (See
[available databases][1] for more information.)

```
curl -O https://s3-us-west-1.amazonaws.com/spacegraphcats.ucdavis.edu/microbe-genbank-sbt-k31-2017.05.09.tar.gz
tar xzf microbe-genbank-sbt-k31-2017.05.09.tar.gz
```

This produces a file `genbank-k31.sbt.json` and a whole bunch of hidden
files in the directory `.sbt.genbank-k31`.

Next, run the 'gather' command to see what's in your ecoli genome --
```
sourmash gather -k 31 ecoli-genome.sig genbank-k31.sbt.json
```

and you should get:

```
# running sourmash subcommand: gather
loaded query: /home/ubuntu/data/ecoliMG1655.... (k=31, DNA)

overlap     p_query p_match
---------   ------- --------
4.9 Mbp     100.0%   99.8%      CP011320.1 Escherichia coli strain SQ...

found 1 matches total;
the recovered matches hit 100.0% of the query
```

In this case, the output is kind of boring because this is a single
genome.  But! You can use this on metagenomes (assembled and
unassembled) as well; you've just got to make the signature files.

To see this in action, here is gather running on a signature generated
from some sequences that assemble (but don't align to known genomes)
from the
[Shakya et al. 2013 mock metagenome paper.][2]

```
wget https://github.com/dib-lab/sourmash/raw/master/doc/_static/shakya-unaligned-contigs.sig
sourmash gather -k 31 shakya-unaligned-contigs.sig genbank-k31.sbt.json
```

This should yield:
```
loaded query: mqc500.QC.AMBIGUOUS.99.unalign... (k=31, DNA)
loaded 0 signatures and 1 databases total.


overlap     p_query p_match 
---------   ------- --------
1.4 Mbp      11.0%   58.0%      JANA01000001.1 Fusobacterium sp. OBRC...
1.0 Mbp       7.7%   25.9%      CP001957.1 Haloferax volcanii DS2 pla...
0.9 Mbp       7.4%   11.8%      BA000019.2 Nostoc sp. PCC 7120 DNA, c...
0.7 Mbp       5.9%   23.0%      FOVK01000036.1 Proteiniclasticum rumi...
0.7 Mbp       5.3%   17.6%      AE017285.1 Desulfovibrio vulgaris sub...
0.6 Mbp       4.9%   11.1%      CP001252.1 Shewanella baltica OS223, ...
0.6 Mbp       4.8%   27.3%      AP008226.1 Thermus thermophilus HB8 g...
0.6 Mbp       4.4%   11.2%      CP000031.2 Ruegeria pomeroyi DSS-3, c...
480.0 kbp     3.8%    7.6%      CP000875.1 Herpetosiphon aurantiacus ...
410.0 kbp     3.3%   10.5%      CH959317.1 Sulfitobacter sp. NAS-14.1...
1.4 Mbp       2.2%   11.8%      LN831027.1 Fusobacterium nucleatum su...
0.5 Mbp       2.1%    5.3%      CP000753.1 Shewanella baltica OS185, ...
420.0 kbp     1.9%    7.7%      FNDZ01000023.1 Proteiniclasticum rumi...
150.0 kbp     1.2%    4.5%      CP015081.1 Deinococcus radiodurans R1...
150.0 kbp     1.2%    8.2%      CP000969.1 Thermotoga sp. RQ2, comple...
290.0 kbp     1.1%    4.1%      CH959311.1 Sulfitobacter sp. EE-36 sc...
1.2 Mbp       1.0%    5.0%      CP013328.1 Fusobacterium nucleatum su...
110.0 kbp     0.9%    3.5%      FREL01000833.1 Enterococcus faecalis ...
0.6 Mbp       0.8%    2.8%      CP000527.1 Desulfovibrio vulgaris DP4...
340.0 kbp     0.6%    3.3%      KQ235732.1 Fusobacterium nucleatum su...
70.0 kbp      0.6%    1.2%      CP000850.1 Salinispora arenicola CNS-...
60.0 kbp      0.5%    0.7%      CP000270.1 Burkholderia xenovorans LB...
50.0 kbp      0.4%    2.6%      CP001080.1 Sulfurihydrogenibium sp. Y...
50.0 kbp      0.4%    3.2%      L77117.1 Methanocaldococcus jannaschi...
found less than 40.0 kbp in common. => exiting

found 24 matches total;
the recovered matches hit 73.1% of the query


```

It is straightforward to build your own databases for use with `search`
and `gather`; ping us if you want us to write that up.

[Return to index][3]

[0]:http://ivory.idyll.org/blog/2016-sourmash-sbt-more.html
[1]:databases.html
[2]:https://www.ncbi.nlm.nih.gov/pubmed/233877
[3]:index.html
