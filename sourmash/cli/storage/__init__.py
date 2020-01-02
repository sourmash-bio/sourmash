from . import convert
from ..utils import command_list
from argparse import SUPPRESS, RawDescriptionHelpFormatter
import os
import sys


def subparser(subparsers):
    subparser = subparsers.add_parser('storage', formatter_class=RawDescriptionHelpFormatter)
    desc = 'Invoke "sourmash storage <subcmd> --help" for more details on executing each subcommand.\n\n'
    desc += '    Operations\n'
    clidir = os.path.dirname(__file__)
    ops = command_list(clidir)
    for subcmd in ops:
        desc += '        sourmash storage {sc:s} --help\n'.format(sc=subcmd)
    s = subparser.add_subparsers(
        title='Subcommands', dest='subcmd', metavar='subcmd', help=SUPPRESS,
        description=desc
    )
    for subcmd in ops:
        getattr(sys.modules[__name__], subcmd).subparser(s)
    subparser._action_groups.reverse()
    subparser._optionals.title = 'Options'
