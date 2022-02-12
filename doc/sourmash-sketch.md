# `sourmash sketch` documentation

Most of the commands in sourmash work with **signatures**, which contain information about genomic or proteomic sequences. Each signature contains one or more **sketches**, which are compressed versions of these sequences. Using sourmash, you can search, compare, and analyze these sequences in various ways.

To create a signature with one or more sketches, you use the `sourmash sketch` command. There are three main commands:

```
sourmash sketch dna
sourmash sketch protein
sourmash sketch translate
```

The `sketch dna` command reads in **DNA sequences** and outputs **DNA sketches**.

The `sketch protein` command reads in **protein sequences** and outputs **protein sketches**.

The `sketch translate` command reads in **DNA sequences**, translates them in all six frames, and outputs **protein sketches**.

All `sourmash sketch` commands take FASTA or FASTQ sequences as input;
input data can be uncompressed, compressed with gzip, or compressed
with bzip2. The output will be one or more JSON signature files that
can be used with the other sourmash commands.

## Quickstart

### DNA sketches for genomes and reads

To create a DNA sketch for a genome, run:
```
sourmash sketch dna genome.fna
```
This will create an output file `genome.fna.sig` in the current directory, containing a single DNA signature for the entire genome, calculated using the default parameters.


Sourmash can work with unassembled reads; run
```
sourmash sketch dna -p k=21,k=31,k=51,abund metagenome.fq.gz
```
to create three abundance-weighted sketches at k=21, 31, and 51, for the given FASTQ file.

By default, `sketch dna` ignores bad k-mers (e.g. non-ACGT characters
in DNA). If `--check-sequence` is provided, `sketch dna` will error
exit on the first bad k-mer.

### Protein sketches for genomes and proteomes

Likewise,
```
sourmash sketch translate genome.fna
```
will output a protein sketch in `./genome.fna.sig`, calculated by translating the genome sequence in all six frames and then using the default protein sketch parameters.  K-mers may include stop codons and stop codons are considered valid protein-coding sequence.

And
```
sourmash sketch protein -p k=25,scaled=500 -p k=27,scaled=250 genome.faa
```
outputs two protein sketches to `./genome.faa.sig`, one calculated with k=25 and scaled=500, the other calculated with k=27 and scaled=250.

If you want to use different encodings, you can specify them in a few ways; here is a parameter string that specifies a dayhoff encoding for the k-mers:
```
sourmash sketch protein -p k=25,scaled=500,dayhoff genome.faa
```

## More detailed documentation

### Input formats

`sourmash sketch` auto-detects and reads FASTQ or FASTA files, either uncompressed or compressed with gzip or bzip2. The filename doesn't matter; `sourmash sketch` will figure out the format from the file contents.

You can also stream any of these formats into `sourmash sketch` via stdin by using `-` as the input filename. For example,
```
gunzip -c data/GCF*.fna.gz | sourmash sketch dna - -o out.sig
```
will make a single DNA signature from all of the FASTA sequences in
`data/GCF*.fna.gz`.

Note, for signatures calculated from stdin, the signature filename attribute
will be left empty, and `sourmash sig describe` will output `** no name **`.

### Input contents and output signatures

By default, `sourmash sketch` will produce signatures for each input
*file*. If the file contains multiple FASTA/FASTQ records, these
records will be merged into the output signature.  You can provide a
*list of FASTA files* in a text file to `sourmash sketch` by passing
the text file path in via `--from-file`.

If you specify `--singleton`, `sourmash sketch` will produce signatures for each *record*.

If you specify `--merge <name>`, sourmash sketch will produce signatures for all input files and combine them into one signature.

The output signature(s) will be saved in locations that depend on your input parameters. By default, `sourmash sketch` will put the signatures in the current directory, in a file named for the input file with a `.sig` suffix. If you specify `-o`, all of the signatures will be placed in that file.

### Protein encodings

`sourmash sketch protein` and `sourmash sketch translate` output protein sketches by default, but can also use the `dayhoff` and `hp` encodings.  The [Dayhoff encoding](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-367/tables/1) collapses multiple amino acids into a smaller alphabet so that amino acids that share biochemical properties map to the same character. The hp encoding divides amino acids into hydrophobic and polar (hydrophilic) amino acids, collapsing amino acids with hydrophobic side chains together and doing the same for polar amino acids.

We are still in the process of benchmarking these encodings; ask [on the issue tracker](https://github.com/sourmash-bio/sourmash/issues) if you are interested in updates.

Note that stop characters (`*`) are considered valid in all three
encodings, and are not truncated. For example, amino acid sequences
that contain stop characters at the end will produce a k-mer containing
the stop character, and that k-mer will be hashed and potentially included
in the sketch.

### Parameter strings

The `-p` argument to `sourmash sketch` provides parameter strings to sourmash, and these control what signatures and sketches are calculated and output. Zero or more parameter strings can be given to sourmash. Each parameter string produces at least one sketch.

A parameter string is a space-delimited collection that can contain one or more fields, comma-separated.
* `k=<ksize>` - create a sketch at this k-mer size; can provide more than one time in a parameter string. Typically `ksize` is between 4 and 100.
* `scaled=<int>` - create a scaled MinHash with k-mers sampled deterministically at 1 per `<scaled>` value. This controls sketch compression rates and resolution; for example, a 5 Mbp genome sketched with a scaled of 1000 would yield approximately 5,000 k-mers. `scaled` is incompatible with `num`. See [our guide to signature resolution](using-sourmash-a-guide.md#what-resolution-should-my-signatures-be-how-should-i-create-them) for more information.
* `num=<int>` - create a standard MinHash with no more than `<num>` k-mers kept. This will produce sketches identical to [mash sketches](https://mash.readthedocs.io/en/latest/). `num` is incompatible with `scaled`. See [our guide to signature resolution](using-sourmash-a-guide.md#what-resolution-should-my-signatures-be-how-should-i-create-them) for more information.
* `abund` / `noabund` - create abundance-weighted (or not) sketches. See [Classify signatures: Abundance Weighting](classifying-signatures.md#abundance-weighting) for details of how this works.
* `dna`, `protein`, `dayhoff`, `hp` - create this kind of sketch. Note that `sourmash sketch dna -p protein` and `sourmash sketch protein -p dna` are invalid; please use `sourmash sketch translate` for the former.

For all field names but `k`, if multiple fields in a parameter string are provided, the last one encountered overrides the previous values. For `k`, if multiple ksizes are specified in a single parameter string, sketches for all ksizes specified are created.

If a field isn't specified, then the default value for that sketch type is used; so, for example, `sourmash sketch dna -p abund` would calculate a sketch with `k=31,scaled=1000,abund`. See below for the defaults.

### Default parameters

The default parameters for sketches are as follows:

* dna: `k=31,scaled=1000,noabund`
* protein: `k=10,scaled=200,noabund`
* dayhoff: `k=16,scaled=200,noabund`
* hp: `k=42,scaled=200,noabund`

These were chosen by a committee of PhDs as being good defaults for an initial analysis, so, beware :).

More seriously, the DNA parameters were chosen based on the analyses done by Koslicki and Falush in [MetaPalette: a k-mer Painting Approach for Metagenomic Taxonomic Profiling and Quantification of Novel Strain Variation](https://msystems.asm.org/content/1/3/e00020-16).

The protein, dayhoff, and hp parameters were selected based on unpublished research results and/or magic formulas. We are working on publishing the results! Please ask on the [issue tracker](https://github.com/sourmash-bio/sourmash/issues) if you are curious.

### More complex parameter string examples

Below are some more complicated `sourmash sketch` command lines:

* `sourmash sketch dna -p k=51` - default to a scaled=1000 and noabund for a k-mer size of 51 (based on moltype/command)
* `sourmash sketch dna -p k=31,k=51,k=21` - create one signature with multiple ksizes, using the defaults otherwise
* `sourmash sketch translate -p k=20,num=500,protein -p k=19,num=400,dayhoff,abund -p k=30,scaled=200,hp` - create three signatures with different ksizes, moltypes, and scaled/num.

### Signature naming

Signature names are displayed in the output for search, gather, and
compare, and can be specified in a few different ways.

With default arguments, `sourmash sketch` does not set a name, and the
filename is used in display output.

You can set a name using `--name`, but this has the side effect of
merging the sequence records before signature creation. So, for example,
`sourmash sketch dna genome1.fa genome2.fa --name genome1 -o
genome.sig` would produce one signature after combining `genome1.fa`
and `genome2.fa`.

The option `--name-from-first` will set the signature name from the
first record header encountered in each file.  When used with `--singleton`,
this will name each signature based on the record that it is created from.

You can examine the signature name using `sourmash sig describe`.

Individual signature renaming can be done from the command line using
`sourmash sig split` to create individual files for each signature,
and then `sourmash sig rename`.

### Locations for output files

Signature files can contain multiple signatures and sketches. Use `sourmash sig describe` to get details on the contents of a file.

You can use `-o <filename>` to specify a file output location for all the output signatures; `-o -` means stdout. This does not merge signatures unless `--merge` is provided.

Specify `--outdir` to put all the signatures in a specific directory.

### Downsampling and flattening signatures

Creating signatures is probably the most time consuming part of using sourmash, and it is the only part that requires access to the raw data. Moreover, the output signatures are generally much smaller than the input data. So, we generally suggest creating a large set of signatures once.

To support this, sourmash can do two kinds of signature conversion without going back to the raw data.

First, you can downsample `num` and `scaled` signatures using `sourmash sig downsample`.  For any sketch created with `num` parameter, you can decrease that `num`. And, for any `scaled` parameter, you can increase the `scaled`. This will decrease the size of the sketch accordingly; for example, going from a `num` of 5000 to a `num` of 1000 will decrease the sketch size by a factor of 5, and going from a `scaled` of 1000 to a `scaled` of 10000 will decrease the sketch size by a factor of 10.

(Note that decreasing `num` or increasing `scaled` will increase calculation speed and lower the accuracy of your results.)

Second, you can flatten abundances using `sourmash sig flatten`. For any sketch created with `abund`, you can convert it to a `noabund` sketch.  This will decrease the sketch size, although not necessarily by a lot.

Unfortunately, changing the k-mer size or using different DNA/protein encodings cannot be done on a sketch, and you need to create new signatures from the raw data for that.

### Examining the output of `sourmash sketch`

You can use `sourmash sig describe` to get detailed information about the contents of a signature file. This can help if you want to see exactly what a particular `sourmash sketch` command does!

### Filing issues and asking for help

We try to provide good documentation and error messages, but may not succeed in answer all your questions! So we're happy to help out!

Please post questions [on the sourmash issue tracker](https://github.com/sourmash-bio/sourmash/issues). If you find something confusing or buggy about the documentation or about sourmash, we'd love to fix it -- for you *and* for everyone else!
