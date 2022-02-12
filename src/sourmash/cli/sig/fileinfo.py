"""provide summary information on the given file"""


def subparser(subparsers):
    subparser = subparsers.add_parser('fileinfo')
    subparser.add_argument('path')
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-f', '--force', action='store_true',
        help='try to load all files as signatures'
    )


def main(args):
    import sourmash
    return sourmash.sig.__main__.fileinfo(args)
