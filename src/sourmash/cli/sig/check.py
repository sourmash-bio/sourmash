"""check signature collections against a picklist"""

usage="""

    sourmash sig check <filenames> --picklist ... -o miss.csv -m manifest.csv

This will check the signature contents of <filenames> against the given
picklist, optionally outputting the unmatched picklist rows to 'miss.csv'
and optionally outputting a manifest of the matched signatures to
'manifest.csv'.

By default, 'sig check' requires a pre-existing manifest for collections;
this prevents potentially slow manifest rebuilding. You
can turn this check off with '--no-require-manifest'.

"""

from sourmash.cli.utils import (add_moltype_args, add_ksize_arg,
                                add_picklist_args, add_pattern_args)


def subparser(subparsers):
    subparser = subparsers.add_parser('check', usage=usage)
    subparser.add_argument('signatures', nargs='*')
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-d', '--debug', action='store_true',
        help='provide debugging output'
    )
    subparser.add_argument(
        '-o', '--output-missing', metavar='FILE',
        help='output picklist with remaining unmatched entries to this file',
    )
    subparser.add_argument(
        '-f', '--force', action='store_true',
        help='try to load all files as signatures'
    )
    subparser.add_argument(
        '--from-file',
        help='a text file containing a list of files to load signatures from'
    )
    subparser.add_argument(
        '-m', '--save-manifest-matching',
        help='save a manifest of the matching entries to this file.'
    )
    subparser.add_argument(
        '--fail-if-missing', action='store_true',
        help='exit with an error code (-1) if there are any missing picklist values.'
    )
    subparser.add_argument(
        '--no-require-manifest',
        help='do not require a manifest; generate dynamically if needed',
        action='store_true'
    )
    subparser.add_argument(
        '-F', '--manifest-format',
        help="format of manifest output file; default is 'csv')",
        default='csv',
        choices=['csv', 'sql'],
    )

    add_ksize_arg(subparser)
    add_moltype_args(subparser)
    add_pattern_args(subparser)
    add_picklist_args(subparser)


def main(args):
    import sourmash
    return sourmash.sig.__main__.check(args)
