"""search a signature against other signatures"""

usage="""

The `search` subcommand searches a collection of signatures or SBTs for
matches to the query signature.  It can search for matches with either
high [Jaccard similarity](https://en.wikipedia.org/wiki/Jaccard_index)
or containment; the default is to use Jaccard similarity, unless
`--containment` is specified.  `-o/--output` will create a CSV file
containing the matches.

`search` will load all of provided signatures into memory, which can
be slow and somewhat memory intensive for large collections.  You can
use `sourmash index` to create a Sequence Bloom Tree (SBT) that can
be quickly searched on disk; this is [the same format in which we provide
GenBank and other databases](databases.md).

Usage:
```
sourmash search query.sig [ list of signatures or SBTs ]
```

Example output:

```
49 matches; showing first 20:
similarity   match
----------   -----
 75.4 percent       NZ_JMGW01000001.1 Escherichia coli
 1-176-05_S4_C2 e117605...
 72.2 percent       NZ_GG774190.1 Escherichia coli MS 196-1 Scfld2538, whole ...
 71.4 percent       NZ_JMGU01000001.1 Escherichia coli 2-011-08_S3_C2 e201108...
 70.1 percent       NZ_JHRU01000001.1 Escherichia coli strain 100854 100854_1...
 69.0 percent       NZ_JH659569.1 Escherichia coli M919 supercont2.1, whole g...
...  
```

Note, as of sourmash 4.2.0, `search` supports `--picklist`, to
[select a subset of signatures to search, based on a CSV file](#using-picklists-to-subset-large-collections-of-signatures). This
can be used to search only a small subset of a large collection, or to
exclude a few signatures from a collection, without modifying the
collection itself.

---
"""

from sourmash.cli.utils import (add_ksize_arg, add_moltype_args,
                                add_picklist_args, add_scaled_arg)


def subparser(subparsers):
    subparser = subparsers.add_parser('search', description=__doc__, usage=usage)
    subparser.add_argument(
        'query', help='query signature'
    )
    subparser.add_argument(
        'databases', nargs='+',
        help='signatures/SBTs to search',
    )
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '--threshold', metavar='T', default=0.08, type=float,
        help='minimum threshold for reporting matches; default=0.08'
    )
    subparser.add_argument(
        '--save-matches', metavar='FILE',
        help='output matching signatures to the specified file'
    )
    subparser.add_argument(
        '--best-only', action='store_true',
        help='report only the best match (with greater speed)'
    )
    subparser.add_argument(
        '-n', '--num-results', default=3, type=int, metavar='N',
        help='number of results to report'
    )
    subparser.add_argument(
        '--containment', action='store_true',
        help='score based on containment rather than similarity'
    )
    subparser.add_argument(
        '--max-containment', action='store_true',
        help='score based on max containment rather than similarity'
    )
    subparser.add_argument(
        '--ignore-abundance', action='store_true',
        help='do NOT use k-mer abundances if present; note: has no effect if '
        '--containment or --max-containment is specified'
    )
    subparser.add_argument(
        '-o', '--output', metavar='FILE',
        help='output CSV containing matches to this file'
    )
    subparser.add_argument(
        '--md5', default=None,
        help='select the signature with this md5 as query'
    )
    add_ksize_arg(subparser, 31)
    add_moltype_args(subparser)
    add_picklist_args(subparser)
    add_scaled_arg(subparser, 0)


def main(args):
    import sourmash
    return sourmash.commands.search(args)
