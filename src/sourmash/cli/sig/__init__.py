"""Define the command line interface for sourmash sig

The top level CLI is defined in ../__init__.py. This module defines the CLI for
`sourmash sig` operations.
"""

import os
import sys
from argparse import SUPPRESS, RawDescriptionHelpFormatter

from ..utils import command_list
from . import (
    cat,
    check,
    collect,
    describe,
    downsample,
    export,
    extract,
    fileinfo,
    filter,
    flatten,
    grep,
    inflate,
    ingest,
    intersect,
    kmers,
    manifest,
    merge,
    overlap,
    rename,
    split,
    subtract,
)
from . import fileinfo as summarize


def subparser(subparsers):
    subparser = subparsers.add_parser(
        "sig",
        formatter_class=RawDescriptionHelpFormatter,
        usage=SUPPRESS,
        aliases=["signature"],
    )
    desc = "Operations\n"
    clidir = os.path.dirname(__file__)
    ops = command_list(clidir)
    for subcmd in ops:
        docstring = getattr(sys.modules[__name__], subcmd).__doc__
        helpstring = f"sourmash sig {subcmd:s} --help"
        desc += f"        {helpstring:33s} {docstring:s}\n"
    s = subparser.add_subparsers(
        title="Manipulate signature files",
        dest="subcmd",
        metavar="subcmd",
        help=SUPPRESS,
        description=desc,
    )
    for subcmd in ops:
        getattr(sys.modules[__name__], subcmd).subparser(s)
    subparser._action_groups.reverse()
    subparser._optionals.title = "Options"
