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
    "Report sourmash version + version of installed dependencies."
    notify(f'sourmash version {sourmash.VERSION}')
    notify(f'- loaded from path: {os.path.dirname(__file__)}')
    notify('')

    if verbose:
        notify('khmer version: None (internal Nodegraph)')
        notify('')

        notify(f'screed version {screed.__version__}')
        notify(f'- loaded from path: {os.path.dirname(screed.__file__)}')


def main(args):
    info(verbose=args.verbose)
