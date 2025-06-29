"""collect manifest information across many files"""

usage = """

    sourmash sig collect <filenames> -o all.sqlmf

This will collect manifests from across many files and save the information
into a standalone manifest database.

By default, 'sig collect' requires a pre-existing manifest for collections;
this prevents potentially slow manifest rebuilding. You
can turn this check off with '--no-require-manifest'.

"""

from sourmash.cli.utils import (
    add_moltype_args,
    add_ksize_arg,
    add_picklist_args,
    add_pattern_args,
    add_v5_args,
)


def subparser(subparsers):
    subparser = subparsers.add_parser("collect", usage=usage)
    subparser.add_argument("locations", nargs="*", help="locations of input signatures")
    subparser.add_argument("-o", "--output", help="manifest output file", required=True)
    subparser.add_argument(
        "-q", "--quiet", action="store_true", help="suppress non-error output"
    )
    subparser.add_argument(
        "-d", "--debug", action="store_true", help="provide debugging output"
    )
    subparser.add_argument(
        "--from-file",
        help="a text file containing a list of files to load signatures from",
    )
    subparser.add_argument(
        "--no-require-manifest",
        help="do not require a manifest; generate dynamically if needed",
        action="store_true",
    )
    subparser.add_argument(
        "-F",
        "--manifest-format",
        help="format of manifest output file; default is 'csv')",
        default="sql",
        choices=["csv", "sql"],
    )

    subparser.add_argument(
        "--merge-previous",
        action="store_true",
        help="merge new manifests into existing",
    )
    subparser.add_argument(
        "--abspath",
        "--use-absolute-paths",
        help="convert all locations to absolute paths",
        action="store_true",
        default=None,  # None: not user specified, so can switch v4/v5
    )
    subparser.add_argument(
        "--no-abspath",
        help="do not convert all locations to absolute paths",
        action="store_false",
        dest="abspath",
    )
    subparser.add_argument(
        "--relpath",
        "--use-relative-paths",
        help="convert all locations to paths relative to the output manifest",
        action="store_true",
    )
    subparser.add_argument(
        "--no-relpath",
        help="do not convert all locations to paths relative to the output manifest",
        action="store_false",
        dest="relpath",
    )

    add_ksize_arg(subparser)
    add_moltype_args(subparser)
    add_v5_args(subparser)


def main(args):
    import sourmash

    return sourmash.sig.__main__.collect(args)
