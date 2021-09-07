"""search a metagenome signature against dbs"""

usage="""
The `gather` subcommand selects the best reference genomes to use for
a metagenome analysis, by finding the smallest set of non-overlapping
matches to the query in a database.  This is specifically meant for
metagenome and genome bin analysis.  (See
[Classifying Signatures](classifying-signatures.md) for more
information on the different approaches that can be used here.)

If the input signature was created with `-p abund`, output
will be abundance weighted (unless `--ignore-abundances` is
specified).  `-o/--output` will create a CSV file containing the
matches.

`gather`, like `search`, will load all of provided signatures into
memory.  You can use `sourmash index` to create a Sequence Bloom Tree
(SBT) that can be quickly searched on disk; this is
[the same format in which we provide GenBank and other databases](databases.md).

Usage:
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
50kb. see the Appendix in
[Classifying Signatures](classifying-signatures.md) for details.

As of sourmash 4.2.0, `gather` supports `--picklist`, to
[select a subset of signatures based on a CSV file](#using-picklists-to-subset-large-collections-of-signatures). This
can be used to search only a small subset of a large collection, or to
exclude a few signatures from a collection, without modifying the
collection itself.

Note:

Use `sourmash gather` to classify a metagenome against a collection of
genomes with no (or incomplete) taxonomic information.  Use `sourmash
lca summarize` to classify a metagenome using a collection of genomes
with taxonomic information.

### Alternative search mode for low-memory (but slow) search: `--linear`

By default, `sourmash gather` uses all information available for
faster search. In particular, for SBTs, `prefetch` will prune the search
tree.  This can be slow and/or memory intensive for very large databases,
and `--linear` asks `sourmash prefetch` to instead use a linear search
across all leaf nodes in the tree.

The results are the same whether `--no-linear` or `--linear` is
used.

### Alternative search mode: `--no-prefetch`

By default, `sourmash gather` does a "prefetch" to find *all* candidate
signatures across all databases, before removing overlaps between the
candidates. In rare circumstances, depending on the databases and parameters
used, this may be slower or more memory intensive than doing iterative
overlap removal. Prefetch behavior can be turned off with `--no-prefetch`.

The results are the same whether `--prefetch` or `--no-prefetch` is
used.  This option can be used with or without `--linear` (although
`--no-prefetch --linear` will generally be MUCH slower).

---
"""

from sourmash.cli.utils import (add_ksize_arg, add_moltype_args,
                                add_picklist_args, add_scaled_arg)


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

    add_ksize_arg(subparser, 31)
    add_moltype_args(subparser)
    add_picklist_args(subparser)
    add_scaled_arg(subparser, 0)


def main(args):
    import sourmash
    return sourmash.commands.gather(args)
