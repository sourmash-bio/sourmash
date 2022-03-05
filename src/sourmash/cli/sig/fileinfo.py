"""provide summary information on the given file"""

usage="""

    sourmash sig fileinfo <filename>

This will provide a summary of the sketch contents in the given file.

JSON output can be generated in place of the normal human-readable output
with '--json-out'.

'sig summarize' and 'sig fileinfo' are aliases for the same command.

"""



def subparser(subparsers):
    subparser = subparsers.add_parser('fileinfo', aliases=['summarize'],
                                      usage=usage)
    subparser.add_argument('path')
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-d', '--debug', action='store_true',
        help='output debug information'
    )
    subparser.add_argument(
        '-f', '--force', action='store_true',
        help='try to load all files as signatures'
    )
    subparser.add_argument(
        '--rebuild-manifest', help='forcibly rebuild the manifest',
        action='store_true'
    )
    subparser.add_argument(
        '--json-out', help='output information in JSON format only',
        action='store_true'
    )


def main(args):
    import sourmash
    return sourmash.sig.__main__.fileinfo(args)
