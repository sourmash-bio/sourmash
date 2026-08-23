"""Define the command line interface for sourmash sketch

The top level CLI is defined in ../__init__.py. This module defines the CLI for
`sourmash sketch` operations.
"""

import os
import sys
from argparse import SUPPRESS, RawDescriptionHelpFormatter

from ..utils import command_list
from . import dna, fromfile, protein, translate
from . import dna as rna
from . import protein as aa
from . import protein as prot


def subparser(subparsers):
    subparser = subparsers.add_parser(
        "sketch", formatter_class=RawDescriptionHelpFormatter, usage=SUPPRESS
    )
    desc = "Operations\n"
    clidir = os.path.dirname(__file__)
    ops = command_list(clidir)
    for subcmd in ops:
        docstring = getattr(sys.modules[__name__], subcmd).__doc__
        helpstring = f"sourmash sketch {subcmd:s} --help"
        desc += f"        {helpstring:33s} {docstring:s}\n"
    s = subparser.add_subparsers(
        title="Create signatures",
        dest="subcmd",
        metavar="subcmd",
        help=SUPPRESS,
        description=desc,
    )
    for subcmd in ops:
        getattr(sys.modules[__name__], subcmd).subparser(s)
    subparser._action_groups.reverse()
    subparser._optionals.title = "Options"
