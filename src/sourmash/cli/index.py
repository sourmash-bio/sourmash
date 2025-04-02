"""index signatures for rapid search"""

usage = """

   sourmash index -k 31 dbname *.sig -F dbtype

Create an on-disk database of signatures that can be searched quickly
& in low memory. All signatures must be scaled, and must be the same
k-mer size and molecule type; the standard signature selectors
(-k/--ksize, --scaled, --dna/--protein) choose which signatures to be
added.

The key options for index are:

 * `-k/--ksize <int>`: k-mer size to select
 * `--dna` or --protein`: nucleotide or protein signatures (default `--dna`)
 * `-F <dbtype>`: 'SBT' (default), 'rocksdb', or 'zip'. 'rocksdb' is recommended and will be come the default in sourmash v5.

---
"""

from sourmash.cli.utils import (
    add_ksize_arg,
    add_moltype_args,
    add_picklist_args,
    add_scaled_arg,
)


def subparser(subparsers):
    subparser = subparsers.add_parser("index", description=__doc__, usage=usage)
    subparser.add_argument(
        "-F",
        "--index-type",
        help="type of index to build (default: SBT)",
        choices=["SBT", "rocksdb", "zip"],
        default="SBT",
    )

    subparser.add_argument(
        "name",
        help="name to save index under; defaults to {name}.sbt.zip",
    )
    subparser.add_argument("signatures", nargs="*", help="signatures to load into SBT")
    subparser.add_argument(
        "--from-file",
        help="a text file containing a list of files to load signatures from",
    )
    subparser.add_argument(
        "-q", "--quiet", action="store_true", help="suppress non-error output"
    )
    subparser.add_argument(
        "-d",
        "--n_children",
        metavar="D",
        type=int,
        default=2,
        help="number of children for internal nodes; default=2",
    )
    subparser.add_argument(
        "--append",
        action="store_true",
        default=False,
        help="add signatures to an existing SBT",
    )
    subparser.add_argument(
        "-x",
        "--bf-size",
        metavar="S",
        type=float,
        default=1e5,
        help="Bloom filter size used for internal nodes",
    )
    subparser.add_argument(
        "-f",
        "--force",
        action="store_true",
        help='try loading *all* files in provided subdirectories, not just .sig files"',
    )
    subparser.add_argument(
        "-s",
        "--sparseness",
        metavar="FLOAT",
        type=float,
        default=0.0,
        help="What percentage of internal nodes will not be saved; ranges "
        "from 0.0 (save all nodes) to 1.0 (no nodes saved)",
    )
    add_ksize_arg(subparser)
    add_moltype_args(subparser)
    add_picklist_args(subparser)
    add_scaled_arg(subparser, 0)


def main(args):
    import sourmash

    return sourmash.commands.index(args)
