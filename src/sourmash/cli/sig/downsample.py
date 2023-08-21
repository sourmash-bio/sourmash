"""downsample one or more signatures"""

usage="""

### `sourmash signature downsample` - decrease the size of a signature

Downsample one or more signatures.

With `downsample`, you can --

* increase the `scaled` value for a signature created with `-p scaled=SCALED`, shrinking it in size;
* decrease the `num` value for a traditional num MinHash, shrinking it in size;
* try to convert a `scaled` signature to a `num` signature;
* try to convert a `num` signature to a `scaled` signature.

For example,

sourmash signature downsample file1.sig file2.sig --scaled 100000 -o downsampled.sig

will output each signature, downsampled to a scaled value of 100000, to
`downsampled.sig`; and

sourmash signature downsample --num 500 scaled_file.sig -o downsampled.sig

will try to convert a scaled MinHash to a num MinHash.

"""

from sourmash.cli.utils import (add_moltype_args, add_ksize_arg,
                                add_picklist_args, add_num_arg)


def subparser(subparsers):
    subparser = subparsers.add_parser('downsample', description=__doc__, usage=usage)
    subparser.add_argument('signatures', nargs="*")
    subparser.add_argument(
        '--scaled', type=int, default=0,
        help='scaled value to downsample to'
    )
    subparser.add_argument(
        '--from-file',
        help='a text file containing a list of files to load signatures from'
    )
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-o', '--output', metavar='FILE',
        help='output signature to this file (default stdout)',
        default='-',
    )
    subparser.add_argument(
        '-f', '--force', action='store_true',
        help='try to load all files as signatures'
    )
    add_ksize_arg(subparser)
    add_moltype_args(subparser)
    add_picklist_args(subparser)
    add_num_arg(subparser, 0)


def main(args):
    import sourmash
    return sourmash.sig.__main__.downsample(args)
