"""Define the command line interface for sourmash signature.

Copy commands over from 'sourmash sig'.

This can be removed once Python 2.7 is no longer supported, in favor of an
'aliases' argument to add_subparser in ../sig/__init__.py.
"""

from ..sig import describe
from ..sig import downsample
from ..sig import extract
from ..sig import filter
from ..sig import flatten
from ..sig import intersect
from ..sig import merge
from ..sig import rename
from ..sig import subtract
from ..sig import ingest
from ..sig import export
from ..sig import overlap
from ..utils import command_list
from argparse import SUPPRESS, RawDescriptionHelpFormatter
import os
import sys


def subparser(subparsers):
    subparser = subparsers.add_parser('signature', formatter_class=RawDescriptionHelpFormatter, usage=SUPPRESS)
    desc = 'Operations\n'
    clidir = os.path.join(os.path.dirname(__file__), '../sig/')
    ops = command_list(clidir)
    for subcmd in ops:
        docstring = getattr(sys.modules[__name__], subcmd).__doc__
        helpstring = 'sourmash signature {op:s} --help'.format(op=subcmd)
        desc += '        {hs:33s} {ds:s}\n'.format(hs=helpstring, ds=docstring)
    s = subparser.add_subparsers(
        title='Manipulate signature files', dest='subcmd', metavar='subcmd', help=SUPPRESS,
        description=desc
    )
    for subcmd in ops:
        getattr(sys.modules[__name__], subcmd).subparser(s)
    subparser._action_groups.reverse()
    subparser._optionals.title = 'Options'
