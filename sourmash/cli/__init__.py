from argparse import ArgumentParser
import sys

import sourmash

from . import utils

# Commands
from . import categorize
from . import compare
from . import compute
from . import gather
from . import info
from . import index
from . import plot
from . import search
from . import watch

# Subcommand groups
from . import lca
from . import sbt
from . import sig


class SourmashParser(ArgumentParser):
    def __init__(self, citation=True, **kwargs):
        super(SourmashParser, self).__init__(**kwargs)
        self.citation = citation
        self._citation_printed = False

    def print_citation(self):
        if self._citation_printed:
            return
        import sourmash
        from sourmash.logging import notify
        notify("== This is sourmash version {version}. ==", version=sourmash.VERSION)
        notify("== Please cite Brown and Irber (2016), doi:10.21105/joss.00027. ==\n")
        self._citation_printed = True

    def parse_args(self, args=None, namespace=None):
        if (args is None and len(sys.argv) == 1) or (args is not None and len(args) == 0):
            self.print_help()
            raise SystemExit(0)
        args = super(SourmashParser, self).parse_args(args=args, namespace=namespace)
        if ('quiet' not in args or not args.quiet) and self.citation:
            self.print_citation()

        # BEGIN: dirty hacks to simultaneously support new and previous interface
        if hasattr(args, 'subcmd') and args.subcmd == 'import':
            args.subcmd = 'ingest'
        if hasattr(args, 'cmd') and args.cmd == 'sbt_combine':
            args.cmd = 'sbt'
            args.subcmd = 'combine'
        if hasattr(args, 'subcmd') and args.subcmd == 'compare_csv':
            args.subcmd = 'compare'
        # END: dirty hacks to simultaneously support new and previous interface
        return args


def get_parser():
    commands = ['compute', 'compare', 'search', 'plot', 'gather', 'index',
                'lca', 'sbt', 'info', 'sig', 'categorize', 'watch']
    commandstr = ' -- '.join(sorted(commands))

    desc = 'Compute, compare, manipulate, and analyze MinHash sketches of DNA sequences.'
    parser = SourmashParser(prog='sourmash', description=desc)
    parser._optionals.title = 'Options'
    parser.add_argument('-v', '--version', action='version', version='sourmash '+ sourmash.VERSION)
    parser.add_argument('-q', '--quiet', action='store_true', help='don\'t print citation information')
    sub = parser.add_subparsers(
        title='Commands', dest='cmd', metavar='cmd', help=commandstr,
        description='Invoke "sourmash <cmd> --help" for more details on executing each command.'
    )
    for cmd in commands:
        getattr(sys.modules[__name__], cmd).subparser(sub)
    # BEGIN: dirty hacks to simultaneously support new and previous interface
    sbt.combine.alt_subparser(sub)
    # END: dirty hacks to simultaneously support new and previous interface
    parser._action_groups.reverse()
    return parser
