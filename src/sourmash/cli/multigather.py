"'sourmash multigather' - gather many signatures against multiple databases."

usage="""

The `multigather` subcommand runs 'gather' for multiple query sequences
against the same collection of sequences.  The main use for multigather
is to amortize the cost of loading databases over many gather queries,
so it is most useful when searching against databases that are slow to load.

Usage:
```
sourmash multigather --query <query1.sig> [<query2.sig> ...] --db <db1> <db2>
```

For each query signature, the following output files are created in the
current working directory:

* <base>.csv - 'gather' CSV output, same as 'gather -o'
* <base>.matches.sig - 'gather' matching sigs, same as 'gather --save-matches'
* <base>.unassigned.sig - 'gather' unassigned hashes, same as
       'gather --output-unassigned'

where 'base' is the basename of the 'source file' from the query, or,
if empty, the md5sum from the signature - use `sourmash sig describe` to
retrieve these.

The following commands:
```
sourmash gather query1.sig db1
sourmash gather query2.sig db1
```
can be turned into a multigather command like so:
```
sourmash multigather --query query1.sig query2.sig --db db1
```

"""

from sourmash.cli.utils import add_ksize_arg, add_moltype_args, add_scaled_arg


def subparser(subparsers):
    subparser = subparsers.add_parser('multigather')
    subparser.add_argument(
        '--query', nargs='*', default=[], action='append',
        help='query signature'
    )
    subparser.add_argument(
        '--query-from-file',
        help='file containing list of signature files to query'
    )
    subparser.add_argument(
        '--db', nargs='+', action='append',
        help='signatures/SBTs to search',
    )
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-d', '--debug', action='store_true'
    )
    subparser.add_argument(
        '--threshold-bp', metavar='REAL', type=float, default=5e4,
        help='threshold (in bp) for reporting results (default=50,000)'
    )
    subparser.add_argument(
        '--ignore-abundance',  action='store_true',
        help='do NOT use k-mer abundances if present'
    )
    subparser.add_argument(
        '--estimate-ani-ci', action='store_true',
        help='also output confidence intervals for ANI estimates'
    )
    subparser.add_argument(
        '--fail-on-empty-database', action='store_true',
        help='stop at databases that contain no compatible signatures'
    )
    subparser.add_argument(
        '--no-fail-on-empty-database', action='store_false',
        dest='fail_on_empty_database',
        help='continue past databases that contain no compatible signatures'
    )
    subparser.set_defaults(fail_on_empty_database=True)

    subparser.add_argument(
        '--output-dir', '--outdir',
        help='output CSV results to this directory',
    )

    add_ksize_arg(subparser)
    add_moltype_args(subparser)
    add_scaled_arg(subparser, 0)


def main(args):
    import sourmash
    return sourmash.commands.multigather(args)
