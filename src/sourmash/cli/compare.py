"""compare sequence signatures made by compute"""

usage="""

The `compare` subcommand compares one or more signatures (created with
`sketch`) using estimated Jaccard index [1] or (if signatures are
created with `-p abund`) the angular similarity [2]).

The default output is a text display of a similarity matrix where each
entry `[i, j]` contains the estimated Jaccard index between input
signature `i` and input signature `j`.  The output matrix can be saved
to a file with `--output` and used with the `sourmash plot` subcommand
(or loaded with `numpy.load(...)`.  Using `--csv` will output a CSV
file that can be loaded into other languages than Python, such as R.

Command line usage:
```
sourmash compare file1.sig [ file2.sig ... ]
```

**Note:** compare by default produces a symmetric similarity matrix that can be used as an input to clustering. With `--containment`, however, this matrix is no longer symmetric and cannot formally be used for clustering.

[1] https://en.wikipedia.org/wiki/Jaccard_index
[2] https://en.wikipedia.org/wiki/Cosine_similarity#Angular_distance_and_similarity

---
"""

from sourmash.cli.utils import (add_ksize_arg, add_moltype_args,
                                add_picklist_args)


def subparser(subparsers):
    subparser = subparsers.add_parser('compare', description=__doc__, usage=usage)
    subparser.add_argument(
        'signatures', nargs='*', help='list of signatures to compare',
        default=[]
    )
    subparser.add_argument(
        '-q', '--quiet', action='store_true', help='suppress non-error output'
    )
    add_ksize_arg(subparser)
    add_moltype_args(subparser)
    subparser.add_argument(
        '-o', '--output', metavar='F',
        help='file to which output will be written; default is terminal '
        '(standard output)'
    )
    subparser.add_argument(
        '--ignore-abundance', action='store_true',
        help='do NOT use k-mer abundances even if present'
    )
    subparser.add_argument(
        '--containment', action='store_true',
        help='calculate containment instead of similarity'
    )
    subparser.add_argument(
        '--max-containment', action='store_true',
        help='calculate max containment instead of similarity'
    )
    subparser.add_argument(
        '--from-file',
        help='a text file containing a list of files to load signatures from'
    )
    subparser.add_argument(
        '-f', '--force', action='store_true',
        help='continue past errors in file loading'
    )
    subparser.add_argument(
        '--csv', metavar='F',
        help='write matrix to specified file in CSV format (with column '
        'headers)'
    )
    subparser.add_argument(
        '-p', '--processes', metavar='N', type=int, default=None,
        help='Number of processes to use to calculate similarity')
    add_picklist_args(subparser)


def main(args):
    import sourmash
    return sourmash.commands.compare(args)
