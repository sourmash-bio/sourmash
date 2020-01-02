from . import classify
from . import compare
from . import gather
from . import index
from . import rankinfo
from . import summarize
from ..utils import command_list
from argparse import SUPPRESS, RawDescriptionHelpFormatter
import os
import sys


def subparser(subparsers):
    subparser = subparsers.add_parser('lca', formatter_class=RawDescriptionHelpFormatter, usage=SUPPRESS)
    desc = 'Operations\n'
    clidir = os.path.dirname(__file__)
    ops = command_list(clidir)
    for subcmd in ops:
        desc += '    sourmash lca {sc:s} --help\n'.format(sc=subcmd)
    s = subparser.add_subparsers(
        title='Lowest common ancestor (LCA) utilities', dest='subcmd', metavar='subcmd', help=SUPPRESS,
        description=desc
    )
    for subcmd in ops:
        getattr(sys.modules[__name__], subcmd).subparser(s)
    subparser._action_groups.reverse()
    subparser._optionals.title = 'Options'
