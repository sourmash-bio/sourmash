# Proposed CLI design

This code demonstrates the mechanics of implementing a CLI with multiple levels of subcommand nesting.
The `cli/` directory here would be placed in the sourmash root directory.

Note that the actual arguments are incomplete for the examples below, and incorrect for all the commands/subcommands not shown below.
Fixing all of that will be a tedious copy/paste job.

## Examples

Invoke `./sourmash` or `./sourmash --help`

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

Invoke `./sourmash compare --help`

```
usage: sourmash compare [-h] [-o OUT] [--ignore-abundance]
                        signatures [signatures ...]

positional arguments:
  signatures            list of signatures

optional arguments:
  -h, --help            show this help message and exit
  -o OUT, --output OUT
  --ignore-abundance    do NOT use k-mer abundances if present
```

Invoke `./sourmash lca --help`

```
usage: sourmash lca [-h] subcmd ...

Subcommands:
  Invoke "sourmash lca <subcmd> --help" for more details on executing each
  subcommand.

  subcmd      classify -- compare -- gather -- index -- rankinfo -- summarize

Options:
  -h, --help  show this help message and exit
```

Invoke `./sourmash lca index --help`

```
usage: sourmash lca index [-h] csv lca_db_out signatures [signatures ...]

positional arguments:
  csv         taxonomy spreadsheet
  lca_db_out  output database name
  signatures  one or more sourmash signatures

optional arguments:
  -h, --help  show this help message and exit
```
