"""search a metagenome signature against dbs"""

usage="""

The `gather` subcommand selects the best reference genomes to use for
a metagenome analysis, by finding the smallest set of non-overlapping
matches to the query in a database.  This is specifically meant for
metagenome and genome bin analysis.  (See "Classifying Signatures" [1]
in the command line documentation for more information on the
different approaches that can be used here.)

If the input signature was created with `-p abund`, output
will be abundance weighted (unless `--ignore-abundances` is
specified).  `-o/--output` will create a CSV file containing the
matches.

`gather`, like `search`, will load all of provided signatures into
memory.  You can use `sourmash index` to create a Sequence Bloom Tree
(SBT) that can be quickly searched on disk; this is the same format in
which we provide GenBank and other databases.

Command line usage:
```
sourmash gather query.sig [ list of signatures or SBTs ]
```

Example output:
```
overlap     p_query p_match
---------   ------- --------
1.4 Mbp      11.0%%  58.0%%     JANA01000001.1 Fusobacterium sp. OBRC...
1.0 Mbp       7.7%%  25.9%%     CP001957.1 Haloferax volcanii DS2 pla...
0.9 Mbp       7.4%%  11.8%%     BA000019.2 Nostoc sp. PCC 7120 DNA, c...
0.7 Mbp       5.9%%  23.0%%     FOVK01000036.1 Proteiniclasticum rumi...
0.7 Mbp       5.3%%  17.6%%     AE017285.1 Desulfovibrio vulgaris sub...
```

The command line option `--threshold-bp` sets the threshold below
which matches are no longer reported; by default, this is set to
50kb. see the Appendix in Classifying Signatures [1] for details.

Note:

Use `sourmash gather` to classify a metagenome against a collection of
genomes with no (or incomplete) taxonomic information.  Use `sourmash
lca summarize` to classify a metagenome using a collection of genomes
with taxonomic information.

[1] https://sourmash.readthedocs.io/en/latest/classifying-signatures.html

---
"""

from sourmash.cli.utils import (add_ksize_arg, add_moltype_args,
                                add_picklist_args, add_scaled_arg,
                                add_pattern_args)


def subparser(subparsers):
    subparser = subparsers.add_parser('gather', description=__doc__, usage=usage)
    subparser.add_argument('query', help='query signature')
    subparser.add_argument(
        'databases', nargs='+',
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
        '-n', '--num-results', default=None, type=int, metavar='N',
        help='number of results to report (default: terminate at --threshold-bp)'
    )
    subparser.add_argument(
        '-o', '--output', metavar='FILE',
        help='output CSV containing matches to this file'
    )
    subparser.add_argument(
        '--save-matches', metavar='FILE',
        help='save gather matched signatures from the database to the '
        'specified file'
    )
    subparser.add_argument(
        '--save-prefetch', metavar='FILE',
        help='save all prefetch-matched signatures from the databases to the '
        'specified file or directory'
    )
    subparser.add_argument(
        '--save-prefetch-csv', metavar='FILE',
        help='save a csv with information from all prefetch-matched signatures '
        'to the specified file'
    )
    subparser.add_argument(
        '--threshold-bp', metavar='REAL', type=float, default=5e4,
        help='reporting threshold (in bp) for estimated overlap with remaining query (default=50kb)'
    )
    subparser.add_argument(
        '--output-unassigned', metavar='FILE',
        help='output unassigned portions of the query as a signature to the '
        'specified file'
    )
    subparser.add_argument(
        '--ignore-abundance',  action='store_true',
        help='do NOT use k-mer abundances if present'
    )
    subparser.add_argument(
        '--md5', default=None,
        help='select the signature with this md5 as query'
    )
    subparser.add_argument(
        '--cache-size', default=0, type=int, metavar='N',
        help='number of internal SBT nodes to cache in memory (default: 0, cache all nodes)'
    )

    # advanced parameters
    subparser.add_argument(
        '--linear', dest="linear", action='store_true',
        help="force a low-memory but maybe slower database search",
    )
    subparser.add_argument(
        '--no-linear', dest="linear", action='store_false',
    )
    subparser.add_argument(
        '--no-prefetch', dest="prefetch", action='store_false',
        help="do not use prefetch before gather; see documentation",
    )
    subparser.add_argument(
        '--prefetch', dest="prefetch", action='store_true',
        help="use prefetch before gather; see documentation",
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

    add_ksize_arg(subparser)
    add_moltype_args(subparser)
    add_picklist_args(subparser)
    add_pattern_args(subparser)
    add_scaled_arg(subparser, 0)


def main(args):
    import sourmash
    return sourmash.commands.gather(args)
