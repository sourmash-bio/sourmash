"""ingest/import a mash or other signature"""


def subparser(subparsers):
    # Dirty hack to simultaneously support new and previous interface
    # If desired, this function can be removed with a major version bump.
    for cmd in ('ingest', 'import'):
        subparser = subparsers.add_parser(cmd)
        subparser.add_argument('--csv', action='store_true',
                               help='import in Mash CSV format')
        subparser.add_argument('filenames', nargs='+')
        subparser.add_argument(
            '-q', '--quiet', action='store_true',
            help='suppress non-error output'
        )
        subparser.add_argument(
            '-o', '--output', metavar='FILE',
            help='output signature to this file (default stdout)'
        )


def main(args):
    import sourmash
    return sourmash.sig.__main__.sig_import(args)
