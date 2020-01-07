"""display sourmash version and other information"""

import os
import screed
import sourmash
from sourmash.logging import notify

def subparser(subparsers):
    subparser = subparsers.add_parser('info')
    subparser.add_argument(
        '-v', '--verbose', action='store_true',
        help='report versions of khmer and screed'
    )


def info(verbose=False):
    notify('sourmash version {}', sourmash.VERSION)
    notify('- loaded from path: {}', os.path.dirname(__file__))
    notify('')

    if verbose:
        import khmer
        notify('khmer version {}', khmer.__version__)
        notify('- loaded from path: {}', os.path.dirname(khmer.__file__))
        notify('')

        notify('screed version {}', screed.__version__)
        notify('- loaded from path: {}', os.path.dirname(screed.__file__))


def main(args):
    info(verbose=args.verbose)
