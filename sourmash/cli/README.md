# Proposed CLI design

This code demonstrates the mechanics of implementing a CLI with multiple levels of subcommand nesting.

Note that at the moment the actual arguments are incomplete or incorrect for most of the sourmash commands/subcommands.
Fixing all of that will be a tedious copy/paste job.

## Examples

Invoke `./cli-sandbox` or `./cli-sandbox --help`

```
usage: sourmash [-h] [-v] cmd ...

Compute, compare, manipulate, and analyze MinHash sketches of DNA sequences.

Commands:
  Invoke "sourmash <cmd> --help" for more details on executing each command.

  cmd            compare -- compute -- gather -- info -- lca -- plot -- sbt --
                 search -- signature

Options:
  -h, --help     show this help message and exit
  -v, --version  show program's version number and exit
```

Invoke `./cli-sandbox compare --help`

```
usage: sourmash compare [-h] [-q] [-k K] [--protein] [--dayhoff] [--hp] [-o F]
                        [--ignore-abundance] [--traverse-directory] [--csv F]
                        [-p N]
                        signatures [signatures ...]

positional arguments:
  signatures            list of signatures to compare

optional arguments:
  -h, --help            show this help message and exit
  -q, --quiet           suppress non-error output
  -k K, --ksize K       k-mer size; default=31
  --protein             choose a protein signature; by default, a nucleotide
                        signature is used
  --dayhoff             build Dayhoff-encoded amino acid signatures
  --hp, --hydrophobic-polar
                        build hydrophobic-polar-encoded amino acid signatures
  -o F, --output F      file to which output will be written; default is
                        terminal (standard output)
  --ignore-abundance    do NOT use k-mer abundances even if present
  --traverse-directory  compare all signatures underneath directories
  --csv F               write matrix to specified file in CSV format (with
                        column headers)
  -p N, --processes N   Number of processes to use to calculate similarity

```

Invoke `./cli-sandbox compute --help`

```
usage: sourmash compute [-h] [--protein] [--dayhoff] [--hp] [-q]
                        [--input-is-protein] [-k KSIZES] [-n NUM_HASHES]
                        [--check-sequence] [-f] [-o OUTPUT] [--singleton]
                        [--merge FILE] [--name-from-first] [--input-is-10x]
                        [--count-valid-reads COUNT_VALID_READS]
                        [--write-barcode-meta-csv WRITE_BARCODE_META_CSV]
                        [-p PROCESSES] [--save-fastas SAVE_FASTAS]
                        [--line-count LINE_COUNT] [--track-abundance]
                        [--scaled SCALED] [--seed SEED] [--randomize]
                        [--license LICENSE] [--rename-10x-barcodes FILE]
                        [--barcodes-file FILE]
                        filenames [filenames ...]

positional arguments:
  filenames             file(s) of sequences

optional arguments:
  -h, --help            show this help message and exit
  --protein             choose a protein signature; by default, a nucleotide
                        signature is used
  --dayhoff             build Dayhoff-encoded amino acid signatures
  --hp, --hydrophobic-polar
                        build hydrophobic-polar-encoded amino acid signatures
  -q, --quiet           suppress non-error output
  --input-is-protein    Consume protein sequences - no translation needed.
  -k KSIZES, --ksizes KSIZES
                        comma-separated list of k-mer sizes; default=21,31,51
  -n NUM_HASHES, --num-hashes NUM_HASHES
                        number of hashes to use in each sketch; default=1500
  --check-sequence      complain if input sequence is invalid
  -f, --force           recompute signatures even if the file exists
  -o OUTPUT, --output OUTPUT
                        output computed signatures to this file
  --singleton           compute a signature for each sequence record
                        individually
  --merge FILE, --name FILE
                        merge all input files into one signature file with the
                        specified name
  --name-from-first     name the signature generated from each file after the
                        first record in the file
  --input-is-10x        input is 10x single cell output folder
  --count-valid-reads COUNT_VALID_READS
                        For 10x input only (i.e input-is-10x flag is True), A
                        barcode is only considered a valid barcode read and
                        its signature is written if number of umis are greater
                        than count-valid-reads. It is used to weed out cell
                        barcodes with few umis that might have been due to
                        false rna enzyme reactions
  --write-barcode-meta-csv WRITE_BARCODE_META_CSV
                        For 10x input only (i.e input-is-10x flag is True),
                        for each of the unique barcodes, Write to a given
                        path, number of reads and number of umis per barcode.
  -p PROCESSES, --processes PROCESSES
                        For 10x input only (i.e input-is-10x flag is True,
                        Number of processes to use for reading 10x bam file
  --save-fastas SAVE_FASTAS
                        For 10x input only (i.e input-is-10x flag is True),
                        save merged fastas for all the unique barcodes to
                        {CELL_BARCODE}.fasta in the absolute path given by
                        this flag, By default, fastas are not saved
  --line-count LINE_COUNT
                        For 10x input only (i.e input-is-10x flag is True),
                        line count for each bam shard
  --track-abundance     track k-mer abundances in the generated signature
  --scaled SCALED       choose number of hashes as 1 in FRACTION of input
                        k-mers
  --seed SEED           seed used by MurmurHash; default=42
  --randomize           shuffle the list of input filenames randomly
  --license LICENSE     signature license. Currently only CC0 is supported.
  --rename-10x-barcodes FILE
                        Tab-separated file mapping 10x barcode name to new
                        name, e.g. with channel or cell annotation label
  --barcodes-file FILE  Barcodes file if the input is unfiltered 10x bam file
```

Invoke `./cli-sandbox lca --help`

```
usage: sourmash lca [-h] subcmd ...

Subcommands:
  Invoke "sourmash lca <subcmd> --help" for more details on executing each
  subcommand.

  subcmd      classify -- compare -- gather -- index -- rankinfo -- summarize

Options:
  -h, --help  show this help message and exit
```

Invoke `./cli-sandbox lca index --help`

```
usage: sourmash lca index [-h] csv lca_db_out signatures [signatures ...]

positional arguments:
  csv         taxonomy spreadsheet
  lca_db_out  output database name
  signatures  one or more sourmash signatures

optional arguments:
  -h, --help  show this help message and exit
```
