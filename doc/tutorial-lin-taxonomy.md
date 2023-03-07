# Analyzing Metagenome Composition using the LIN taxonomic framework  

Tessa Pierce Ward

March 2023

requires sourmash v4.8+

---

```{contents}
   :depth: 2
```

This tutorial uses the `sourmash taxonomy` module, which was introduced via [blog post](https://bluegenes.github.io/sourmash-tax/)
and was recently shown to perfom well for taxonomic profiling of long (and short) reads in [Evaluation of taxonomic classification and profiling methods for long-read shotgun metagenomic sequencing datasets](https://link.springer.com/article/10.1186/s12859-022-05103-0), Portik et al., 2022.


In this tutorial, we'll use sourmash gather to analyze metagenomes using the [LIN taxonomic framework](https://dl.acm.org/doi/pdf/10.1145/3535508.3545546).
Specifically, we will analyze plant metagenomes with a low-level pathogen spike-in.
The goal is to see if we can correctly assign the pathogen sequence to its LINgroup, which includes
all known pathogenic strains.

- `barcode1` - highest spike-in (75 picogram/microliter pathogen DNA)
- `barcode3` - lower spike-in (7.5 picogram/microliter pathogen DNA)
- `barcode5` - no spike-in

The pathogen is `Ralstonia solanacearum` in the `Phylum IIB sequevar 1` group. 

This data is courtesy of [The Laboratory of Plant & Atmospheric Microbiology & (Meta)Genomics](https://sites.google.com/vt.edu/lab-vinatzer/home) in collaboration with USDA APHIS.

## Install sourmash

First, we need to install the software! We'll use conda/mamba to do this.

The below command installs [sourmash](http://sourmash.readthedocs.io/).

Install the software:
```
# create a new environment
mamba create -n smash -y -c conda-forge -c bioconda sourmash
```

then activate the conda environment:
```
conda activate smash
```

> Victory conditions: your prompt should start with
> `(smash) `
> and you should now be able to run `sourmash` and have it output usage information!!

## Create a working subdirectory

Make a directory named `smash_lin`, change into it:
```
mkdir -p ~/smash_lin
cd ~/smash_lin
```

Now make a couple useful folders:
```
mkdir -p inputs
mkdir -p databases
```

## Download relevant data

### First, download a database and taxonomic information

Here, we know the spike-in is a pathogenic seqevar of Ralstonia. We will download a database
containing signatures of 27 Ralstonia genomes (pathogenic and not) and the corresponding taxonomic and lingroup information.

```
# database
curl -JLO https://osf.io/vxsta/download
mv ralstonia*.zip ./databases/ralstonia.zip

# taxonomy csv
curl -JLO https://raw.githubusercontent.com/bluegenes/2023-demo-sourmash-LIN/main/databases/ralstonia-lin.taxonomy.GCA-GCF.csv
mv ralstonia-lin.taxonomy.GCA-GCF.csv ./databases

# lingroup csv
curl -JLO https://raw.githubusercontent.com/bluegenes/2023-demo-sourmash-LIN/main/inputs/ralstonia.lingroups.csv
mv ralstonia.lingroups.csv ./databases

ls databases # look at the database files
```

### Next, download pre-made sourmash signatures made from the input metagenomes

```
# download barcode 1 sig
curl -JLO https://osf.io/ujntr/download
mv barcode1_22142.sig.zip ./inputs/

# download barcode 3 signature
curl -JLO https://osf.io/2h9wx/download
mv barcode3_31543.sig.zip ./inputs

# download barcode 5 signature
curl -JLO https://osf.io/k8nw5/download
mv barcode5_36481.sig.zip ./inputs

# look at available input files
ls inputs
```

## Look at the signatures

Let's start with the `barcode1` (highest spike-in) sample

### First, let's look at the metagenome signature.

By running `sourmash sig fileinfo`, we can see information on the signatures available within the zip file.

Here, you can see I've generated the metagenome signature with `scaled=1000` and built two ksizes, `k=31` and `k=51`

Run:
```
sourmash sig fileinfo ./inputs/barcode1_22142.sig.zip
```

In the output, you should see:
```
** loading from './inputs/barcode1_22142.sig.zip'
path filetype: ZipFileLinearIndex
location: /home/jovyan/smash_lin/inputs/barcode1_22142.sig.zip
is database? yes
has manifest? yes
num signatures: 2
total hashes: 914328
summary of sketches:
   1 sketches with DNA, k=31, scaled=1000, abund      426673 total hashes
   1 sketches with DNA, k=51, scaled=1000, abund      487655 total hashes
```

### We can also look at the database

Here, you can see I've generated the database with `scaled=1000` and built three ksizes, `k=21`, `k=31` and `k=51`

Run:
```
sourmash sig fileinfo ./databases/ralstonia.zip
```

In the output, you should see:

```
** loading from './databases/ralstonia.zip'
path filetype: ZipFileLinearIndex
location: /home/jovyan/databases/ralstonia.zip
is database? yes
has manifest? yes
num signatures: 81
** examining manifest...
total hashes: 445041
summary of sketches:
   27 sketches with DNA, k=21, scaled=1000, abund     148324 total hashes
   27 sketches with DNA, k=31, scaled=1000, abund     148111 total hashes
   27 sketches with DNA, k=51, scaled=1000, abund     148606 total hashes
```
There's a lot of things to digest in this output but the two main ones are:
* there are 27 genomes represented in this database, each of which are sketched at k=21,k=31,k=51
* this database represents ~445 *million* k-mers (multiply number of hashes by the scaled number)


## Run sourmash gather using ksize 51

Now let's run `sourmash gather` to find the closest reference genome(s) in the database.
If you want to read more about what sourmash is doing, please see [Lightweight compositional analysis of metagenomes with FracMinHash and minimum metagenome covers](https://www.biorxiv.org/content/10.1101/2022.01.11.475838v2), Irber et al., 2022.

Run:
```
query="inputs/barcode1_22142.sig.zip"
database="databases/ralstonia.zip"

gather_csv_output="barcode1_22141.k51.gather.csv"

sourmash gather $query $database -k 51 -o $gather_csv_output
```

You should see the following output:
```
selecting specified query k=51
loaded query: barcode1_22142... (k=51, DNA)
--ading from 'databases/ralstonia.zip'...
loaded 81 total signatures from 1 locations.
after selecting signatures compatible with search, 27 remain.
Starting prefetch sweep across databases.

Found 7 signatures via prefetch; now doing gather.

overlap     p_query p_match avg_abund
---------   ------- ------- ---------
105.0 kbp      0.0%    2.0%       1.0    GCA_002251655.1 Ralstonia solanacear...
found less than 50.0 kbp in common. => exiting

found 1 matches total;
the recovered matches hit 0.0% of the abundance-weighted query.
the recovered matches hit 0.0% of the query k-mers (unweighted).
```

The first step of gather found all potential matches (7), and the greedy algorithm narrowed this to a single best match, `GCA_002251655.1` which shared an estimated 105 kbp with the metagenome (a very small percentage of the total dataset.) This is expected, though, since the dataset is a plant metagenome with a small `Ralstonia` spike-in.

## Add taxonomic information and summarize up lingroups

`sourmash gather` finds the smallest set of reference genomes that contains all the known information (k-mers) in the metagenome. In most cases, `gather` will find many metagenome matches. Here, we're only looking for `Ralstonia` matches and we only have a single gather result. Regardless, let's use `sourmash tax metagenome` to add taxonomic information and see if we've correctly assigned the pathogenic sequence.

### First, let's look at the relevant taxonomy files.

These commands will show the first few lines of each file. If you prefer, you can look at a more human-friendly view by opening the files in a spreadsheet program.

- **taxonomy_csv:** `databases/ralstonia-lin.taxonomy.GCA-GCF.csv`
  - the essential columns are `lin` (`14;1;0;...`) and `ident` (`GCF_00`...)
- **lingroups information:** `databases/ralstonia.lingroups.csv`
  - both columns are essential (`name`, `lin`)


Look at the taxonomy file:
```
head -n 5 databases/ralstonia-lin.taxonomy.GCA-GCF.csv
```

You should see:
```
lin,species,strain,filename,accession,ident
14;1;0;0;0;0;0;0;0;0;6;0;1;0;1;0;0;0;0;0,Ralstonia solanacearum,OE1_1,GCF_001879565.1_ASM187956v1_genomic.fna,GCF_001879565.1,GCF_001879565.1
14;1;0;0;0;0;0;0;0;0;6;0;1;0;0;0;0;0;0;0,Ralstonia solanacearum,PSS1308,GCF_001870805.1_ASM187080v1_genomic.fna,GCF_001870805.1,GCF_001870805.1
14;1;0;0;0;0;0;0;0;0;2;1;0;0;0;0;0;0;0;0,Ralstonia solanacearum,FJAT_1458,GCF_001887535.1_ASM188753v1_genomic.fna,GCF_001887535.1,GCF_001887535.1
14;1;0;0;0;0;0;0;0;0;2;0;0;4;4;0;0;0;0;0,Ralstonia solanacearum,Pe_13,GCF_012062595.1_ASM1206259v1_genomic.fna,GCF_012062595.1,GCF_012062595.1
```
> The key columns are:
> - `ident`, containing identifiers matching the database sketches
> - `lin`, containing the species information.

Now, let's look at the lingroups file
```
head -n5 databases/ralstonia.lingroups.csv
```

You should see:
```
name,lin
Phyl II,14;1;0;0;0;3;0
Phyl IIA,14;1;0;0;0;3;0;1;0;0
Phyl IIB,14;1;0;0;0;3;0;0
Phyl IIB seq1 and seq2,14;1;0;0;0;3;0;0;0;0;1;0;0;0;0
```
> Here, we have two columns:
> - `name` - the name for each lingroup. 
> - `lin` - the LIN prefix corresponding to each group.


### Now, run `sourmash tax metagenome` to integrate taxonomic information into `gather` results

Using the `gather` output we generated above, we can integrate taxonomic information and summarize up "ranks" (lin positions). We can produce several different types of outputs, including a `lingroup` report.

`lingroup` format summarizes the taxonomic information at each `lingroup`, and produces a report with 4 columns: 
- `name` (from lingroups file)
- `lin` (from lingroups file)
- `percent_containment` - total % of the file matched to this lingroup
- `num_bp_contained` - estimated number of bp matched to this lingroup

> Since sourmash assigns all k-mers to individual genomes, no reads/base pairs are "assigned" to higher taxonomic ranks or lingroups (as with Kraken-style LCA). Here, "percent_containment" and "num_bp_contained" is calculated by summarizing the assignments made to all genomes in a lingroup. This is akin to the "contained" information in Kraken-style reports.

Run `tax metagenome`:
```
gather_csv_output="barcode1_22141.k51.gather.csv"
taxonomy_csv="databases/ralstonia-lin.taxonomy.GCA-GCF.csv"
lingroups_csv="databases/ralstonia.lingroups.csv"

sourmash tax metagenome -g $gather_csv_output -t $taxonomy_csv \
                        --lins --lingroup $lingroups_csv
```

You should see:
```
loaded 1 gather results from 'barcode1_22141.k51.gather.csv'.
loaded results for 1 queries from 1 gather CSVs
Starting summarization up rank(s): 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0 
Read 11 lingroup rows and found 11 distinct lingroup prefixes.
```

and the results:
```
name	lin	percent_containment	num_bp_contained
Phyl II	14;1;0;0;0;3;0	0.02	108000
Phyl IIB	14;1;0;0;0;3;0;0	0.02	108000
Phyl IIB seq1 and seq2	14;1;0;0;0;3;0;0;0;0;1;0;0;0;0	0.02	108000
IIB seq1	14;1;0;0;0;3;0;0;0;0;1;0;0;0;0;0;0	0.02	108000
```

Here, the most specific lingroup we assign to is `Phyl IIB seq1`, which is the pathogenic lingroup that was spiked in, yay! Note that the other groups in the output all contain this group.


#### Now output the lingroup report to a file (instead of to the terminal)

use `-o` to provide an output basename for taxonomic output.

```
gather_csv_output="barcode1_22141.k51.gather.csv"
taxonomy_csv="databases/ralstonia-lin.taxonomy.GCA-GCF.csv"
lingroups_csv="databases/ralstonia.lingroups.csv"

sourmash tax metagenome -g $gather_csv_output -t $taxonomy_csv \
                        --lins --lingroup $lingroups_csv \
                        -o "barcode1"
```

> You should see `saving 'lingroup' output to 'barcode1.lingroup.tsv'` in the output.

#### Optionally, write multiple output formats

You can use `-F` to specify additional output formats. Here, I've added `csv_summary`. Note that while the `lingroup` format will be generated automatically if you specify the `--lingroup` file, you can also specify it with `-F lingroup` if you want, as I've done here.

Run:
```
gather_csv_output="barcode1_22141.k51.gather.csv"
taxonomy_csv="databases/ralstonia-lin.taxonomy.GCA-GCF.csv"
lingroups_csv="databases/ralstonia.lingroup.csv"

sourmash tax metagenome -g $gather_csv_output -t $taxonomy_csv \
                        --lins --lingroup $lingroups_csv \
                        -F lingroup csv_summary -o "barcode1"
```


You should see the following in the output:

```
saving 'csv_summary' output to 'barcode1.summarized.csv'.
saving 'lingroup' output to 'barcode1.lingroup.txt'.
```

The `csv_summary` format is the **full** summary of this sample, e.g. the summary at each taxonomic rank (LIN position). It also includes an entry with the `unclassified` portion at each rank.

> Note: Multiple output formats require the `-o` `--output-base` to be specified, as each must be written to a file.

Here's an abbreviated version of the `gather` results for  `barcode1`, with lingroup information added:

|         | **ksize** | **scaled** | **best overlap** | **gather match(es)** | **lingroup** | **lin**                            |
| ------- | --------- | ---------- | ---------------- | -------------------- | ------------ | ---------------------------------- |
| **bc1** | 51        | 1000       | 105 kb           | GCA_002251655.1      | IIB seq1     | 14;1;0;0;0;3;0;0;0;0;1;0;0;0;0;0;0 |
| **bc1** | 31        | 1000       | 173 kb           | GCA_002251655.1      | IIB seq1     | 14;1;0;0;0;3;0;0;0;0;1;0;0;0;0;0;0 |


### Now run with `barcode3` sample

#### sourmash gather
Run:
```
query="inputs/barcode3_31543.sig.zip"
database="databases/ralstonia.zip"

gather_csv_output="barcode3_31543.dna.k51.gather.csv"

sourmash gather $query $database -k 51 -o $gather_csv_output
```

You should see:
```
selecting specified query k=51
loaded query: barcode3_31543... (k=51, DNA)
loading from 'databases/ralstonia.zip'...
loaded 81 total signatures from 1 locations.
after selecting signatures compatible with search, 27 remain.
Starting prefetch sweep across databases.
Found 0 signatures via prefetch; now doing gather.
found less than 50.0 kbp in common. => exiting


found 0 matches total;
the recovered matches hit 0.0% of the query k-mers (unweighted).
```

#### gather found no sequence matches! But, we can lower the detection threshold:

```
query="inputs/barcode3_31543.sig.zip"
database="databases/ralstonia.zip"
gather_csv_output="barcode3_31543.k51.gather.csv"

# use a 10kb detection threshold
sourmash gather $query $database -k 51 --threshold-bp 10000 -o $gather_csv_output
```

This time, you should see:
```
selecting specified query k=51
loaded query: barcode3_31543... (k=51, DNA)
loading from 'databases/ralstonia.zip'...
loaded 81 total signatures from 1 locations.
after selecting signatures compatible with search, 27 remain.
Starting prefetch sweep across databases.
Found 6 signatures via prefetch; now doing gather.


overlap     p_query p_match avg_abund
---------   ------- ------- ---------
12.0 kbp       0.0%    0.2%       1.0    GCA_000750575.1 Ralstonia solanacear...

found 1 matches total;
the recovered matches hit 0.0% of the abundance-weighted query.
the recovered matches hit 0.0% of the query k-mers (unweighted).

```

You'll notice that while we have an estimated ~12kbp overlap, the matched genome (`GCA_000750575.1`) is different from the one matched above for `barcode5`. If you run `sourmash tax metagenome` on this output, you'll see that this genome belongs to `Phyl IIB seq 2` group, which is a sister group to the correct `Phyl IIB seq 1` group that we expected. So we have a match but it's not the right one -- why not?

### What happened? Use `prefetch` to investigate

`sourmash gather` has two steps: first, it runs a `prefetch` to find ALL genome matches, and then uses a greedy approach to select the smallest set of genomes that contain ('cover') all known sequence content. Let's run `prefetch` independently so we can look at the results of the first step. Here, let's use `--threshold-bp 0` to get all possible matches.

Run:
```
query="inputs/barcode3_31543.sig.zip"
prefetch_csv_output="barcode3_31543.k51.prefetch.csv"
database="databases/ralstonia.zip"

sourmash prefetch $query $database -k 51 --threshold-bp 0 -o $prefetch_csv_output 
```

You should see:
```
selecting specified query k=51
loaded query: barcode3_31543... (k=51, DNA)
query sketch has scaled=1000; will be dynamically downsampled as needed.
--tal of 10 matching signatures so far.tonia.zip'
loaded 81 total signatures from 1 locations.
after selecting signatures compatible with search, 27 remain.
--
total of 15 matching signatures.
saved 15 matches to CSV file 'barcode3_31543.k51.prefetch.csv'
of 487043 distinct query hashes, 12 were found in matches above threshold.
a total of 487031 query hashes remain unmatched.
final scaled value (max across query and all matches) is 1000
```

Here, the output is telling us we found matches to 15 of the 27 Ralstonia genomes. But only **12 k-mers** were shared between the metagenome sample and the genomes. Remember that sourmash uses a representative subsample of all k-mers, so here these 12 k-mers represent ~ 12kb of sequence (12 * scaled). We've found that this is sufficient to detect presence of an organism, but at this low level, it can be hard to distinguish between closely-related genomes. Let's open the prefetch output to see how those 12 k-mers matched between different genomes.

#### Look at the `barcode3_31543.k51.prefetch.csv` file

> Use a spreadsheet program on your computer or use `less -S barcode3_31543.k51.prefetch.csv` to see the file on the terminal. If using `less`, hit `q` when you want to exit and return to your terminal prompt.

The first column contains the estimated number of base pairs matched between our query and each matching reference genome. You'll notice there are four genomes that match 12kb of sequence, one of which is the "correct" genome (`GCA_002251655.1`, which is in the `IIB seq1` lingroup).

**What is happening here?**

When faced with equally good matches, `sourmash gather` makes a random choice about which genome to assign these k-mers to. This happens primarily with highly similar genomes and/or very small sequence matches. If this happens and you need to distinguish between these genomes, we recommend trying a lower scaled value.

To see if we could robustly assign the correct sequevar for `barcode3` using a higher resolution sketch, I also ran `gather` using scaled=100.

Here's an abbreviated version of the `gather` results for  `barcode3`, with lingroup information added:



|         | **ksize** | **scaled** | **best overlap** | **gather match(es)** | **lingroup** | **lin**                            |
| ------- | --------- | ---------- | ---------------- | -------------------- | ------------ | ---------------------------------- |
| **bc3** | 51        | 1000       | 12kb             | GCA_000750575.1      | IIB seq2     | 14;1;0;0;0;3;0;0;0;0;1;0;0;0;0;1;0 |
| **bc3** | 31        | 1000       | 28 kb            | GCA_002251655.1      | IIB seq1     | 14;1;0;0;0;3;0;0;0;0;1;0;0;0;0;0;0 |
| **bc3** | 51        | 100        | 14.8 kb          | GCA_002251655.1      | IIB seq1     | 14;1;0;0;0;3;0;0;0;0;1;0;0;0;0;0;0 |
| **bc3** | 31        | 100        | 21.1 kb          | GCA_002251655.1      | IIB seq1     | 14;1;0;0;0;3;0;0;0;0;1;0;0;0;0;0;0 |

We typically use k=51 for strain-level matching and k=31 for species-level matching. Notice that running at k=31 with scaled 1000 found the right match. However, if you run prefetch for this sample, you see there are three matches with `28kb` overlap, so we just got lucky that `gather` selected the "right" one for this case.

In contrast, running at `scaled=100` had sufficient information to correctly assign the sequence to the `IIB seq1` lingroup.


### Now try the `barcode5` sample

You can also run the `barcode5` file using the same commands as above and see that no matches are found. If you drop the threshold-bp  to 0 (`--threshold-bp 0`), you can find ~1kbp overlap (a single k-mer match!). **Note, we do not recommend trusting/using results with fewer than 3 k-mer matches (3kbp at scaled=1000)**.

I then ran this file at higher resolution to see how the results changed. In each case, very few k-mers matched and we could not robustly identify a specific `Ralstonia` genome or lingroup. As it turns out, `barcode5` does not have a `Ralstonia` spike-in, so this is a good thing!

Here's an abbreviated version of the `gather` results for  `barcode5`, with lingroup information added in cases with a single gather match:

|         | **ksize** | **scaled** | **best overlap** | **gather match(es)** | **lingroup** | **lin**                            |
| ------- | --------- | ---------- | ---------------- | -------------------- | ------------ | ---------------------------------- |
| **bc5** | 51        | 1000       | 1 kbp            | GCA_000750575.1      | IIB seq2     | 14;1;0;0;0;3;0;0;0;0;1;0;0;0;0;1;0 |
| **bc5** | 31        | 1000       | 0                | N/A                  |              |                                    |
| **bc5** | 51        | 100        | 300bp            | all                  |              |                                    |
| **bc5** | 31        | 100        | 1.2 kb           | all                  |              |                                    |
| **bc5** | 51        | 10         | 120 bp           | all                  |              |                                    |
| **bc5** | 31        | 10         | 670 bp           | all                  |              |                                    |
| **bc5** | 51        | 5          | 150 bp           | all                  |              |                                    |
| **bc5** | 31        | 5          | 500 bp           | all                  |              |                                    |


**Again, while I've used a threshold-bp of 0 to get the gather match at scaled=1000, we do not typically trust gather matches with less than `3*scaled` overlap (< 3 k-mers matched).**

## Summary and concluding thoughts

The LIN taxonomic framework may be useful distinguishing groups below the species level.
We can now use LINs and lingroups with `sourmash tax metagenome`. For low level matches, the gather greedy
approach can struggle. We are working on ways to better warn users about this behavior and welcome
feedback and suggestions on our [issue tracker](https://github.com/sourmash-bio/sourmash/issues/new).