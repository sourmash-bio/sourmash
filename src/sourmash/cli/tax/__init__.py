"""Define the command line interface for sourmash tax

The top level CLI is defined in ../__init__.py. This module defines the CLI for
`sourmash tax` operations.
"""

import os
import sys
from argparse import SUPPRESS, RawDescriptionHelpFormatter

from ..utils import command_list
from . import annotate, genome, grep, metagenome, prepare, summarize


def subparser(subparsers):
    subparser = subparsers.add_parser(
        "tax",
        formatter_class=RawDescriptionHelpFormatter,
        usage=SUPPRESS,
        aliases=["taxonomy"],
    )
    desc = "Operations\n"
    clidir = os.path.dirname(__file__)
    ops = command_list(clidir)
    for subcmd in ops:
        docstring = getattr(sys.modules[__name__], subcmd).__doc__
        helpstring = f"sourmash tax {subcmd:s} --help"
        desc += f"        {helpstring:33s} {docstring:s}\n"
    s = subparser.add_subparsers(
        title="Integrate taxonomy information based on 'gather' results",
        dest="subcmd",
        metavar="subcmd",
        help=SUPPRESS,
        description=desc,
    )
    for subcmd in ops:
        getattr(sys.modules[__name__], subcmd).subparser(s)
    subparser._action_groups.reverse()
    subparser._optionals.title = "Options"
