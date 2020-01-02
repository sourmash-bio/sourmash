from . import describe
from . import downsample
from . import extract
from . import filter
from . import flatten
from . import intersect
from . import merge
from . import rename
from . import subtract
from . import ingest
from . import export
from . import overlap
from ..utils import command_list
from argparse import SUPPRESS, RawDescriptionHelpFormatter
import os
import sys


def subparser(subparsers):
    subparser = subparsers.add_parser('sig', formatter_class=RawDescriptionHelpFormatter, usage=SUPPRESS)
    desc = 'Operations\n'
    clidir = os.path.dirname(__file__)
    ops = command_list(clidir)
    for subcmd in ops:
        desc += '    sourmash sig {sc:s} --help\n'.format(sc=subcmd)
    s = subparser.add_subparsers(
        title='Utilities for handling signatures', dest='subcmd', metavar='subcmd', help=SUPPRESS,
        description=desc
    )
    for subcmd in ops:
        getattr(sys.modules[__name__], subcmd).subparser(s)
    subparser._action_groups.reverse()
    subparser._optionals.title = 'Options'
