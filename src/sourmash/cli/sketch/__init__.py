"""Define the command line interface for sourmash sketch

The top level CLI is defined in ../__init__.py. This module defines the CLI for
`sourmash sketch` operations.
"""

from . import dna
from . import dna as rna
from . import protein
from . import protein as aa
from . import protein as prot
from . import translate
from . import fromfile
from ..utils import command_list
from argparse import SUPPRESS, RawDescriptionHelpFormatter
import os
import sys


def subparser(subparsers):
    subparser = subparsers.add_parser('sketch', formatter_class=RawDescriptionHelpFormatter, usage=SUPPRESS)
    desc = 'Operations\n'
    clidir = os.path.dirname(__file__)
    ops = command_list(clidir)
    for subcmd in ops:
        docstring = getattr(sys.modules[__name__], subcmd).__doc__
        helpstring = 'sourmash sketch {op:s} --help'.format(op=subcmd)
        desc += '        {hs:33s} {ds:s}\n'.format(hs=helpstring, ds=docstring)
    s = subparser.add_subparsers(
        title='Create signatures', dest='subcmd', metavar='subcmd', help=SUPPRESS,
        description=desc
    )
    for subcmd in ops:
        getattr(sys.modules[__name__], subcmd).subparser(s)
    subparser._action_groups.reverse()
    subparser._optionals.title = 'Options'
