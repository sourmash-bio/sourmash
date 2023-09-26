"""ingest/import a mash or other signature"""

usage="""

   sourmash sig ingest --csv <input filename> [ <more inputs> ] -o <output>

Ingest num sketches from a simple CSV format, or alternatively a JSON
formatproduced by 'mash info -d'.  The CSV file should contain one
line per sketch, with the first column containing 'murmur64', the
second being '42', the third and fourth being the k-mer size and the
name, and the remaining columns being the hashes.

"""


def subparser(subparsers):
    # Dirty hack to simultaneously support new and previous interface
    # If desired, this function can be removed with a major version bump.
    for cmd in ('ingest', 'import'):
        subparser = subparsers.add_parser(cmd, usage=usage)
        subparser.add_argument('--csv', action='store_true',
                               help='import in Mash CSV format')
        subparser.add_argument('filenames', nargs='+')
        subparser.add_argument(
            '-q', '--quiet', action='store_true',
            help='suppress non-error output'
        )
        subparser.add_argument(
            '-o', '--output', metavar='FILE', default='-',
            help='output signature to this file (default stdout)'
        )


def main(args):
    import sourmash
    return sourmash.sig.__main__.ingest(args)
