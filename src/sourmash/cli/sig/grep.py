"""extract one or more signatures by substr/regex match"""

usage="""
    sourmash sig grep <pattern> <filename> [... <filenames>]

This will search for the provided pattern in the files or databases,
using the signature metadata, and output matching signatures.
Currently 'grep' searches the 'name', 'filename', and 'md5' fields as
displayed by `sig describe`.

'pattern' can be a string or a regular expression.

'sig grep' uses the built-in Python regexp module, 're', to implement
regexp searching. See https://docs.python.org/3/howto/regex.html and
https://docs.python.org/3/library/re.html for details.

The '-v' (exclude), '-i' (case-insensitive), and `-c` (count) options
of 'grep' are supported.

'-o/--output' can be used to output matching signatures to a specific
location.

By default, 'sig grep' requires a pre-existing manifest for collections;
this prevents potentially slow manifest rebuilding. You
can turn this check off with '--no-require-manifest'.

"""

from sourmash.cli.utils import (add_moltype_args, add_ksize_arg,
                                add_picklist_args)


def subparser(subparsers):
    subparser = subparsers.add_parser('grep', usage=usage)
    subparser.add_argument('pattern', help='search pattern (string/regex)')
    subparser.add_argument('signatures', nargs='*')
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-d', '--debug', action='store_true',
        help='output debug information'
    )
    subparser.add_argument(
        '-o', '--output', metavar='FILE',
        help='output matching signatures to this file (default stdout)',
        default='-',
    )
    subparser.add_argument(
        '-f', '--force', action='store_true',
        help='try to load all files as signatures, independent of filename'
    )
    subparser.add_argument(
        '--from-file',
        help='a text file containing a list of files to load signatures from'
    )
    subparser.add_argument(
        '-v', '--invert-match',
        help="select non-matching signatures",
        action="store_true"
    )
    subparser.add_argument(
        '-i', '--ignore-case',
        help="ignore case distinctions (search lower and upper case both)",
        action="store_true"
    )
    subparser.add_argument(
        '--no-require-manifest',
        help='do not require a manifest; generate dynamically if needed',
        action='store_true'
    )
    subparser.add_argument(
        '--csv',
        help='save CSV file containing signature data in manifest format'
    )
    subparser.add_argument(
        '--silent', '--no-signatures-output',
        help="do not output signatures",
        action='store_true',
    )
    subparser.add_argument(
        '-c', '--count',
        help="only output a count of discovered signatures; implies --silent",
        action='store_true'
    )
    add_ksize_arg(subparser)
    add_moltype_args(subparser)
    add_picklist_args(subparser)


def main(args):
    import sourmash.sig.grep
    return sourmash.sig.grep.main(args)
