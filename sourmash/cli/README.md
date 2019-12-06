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
