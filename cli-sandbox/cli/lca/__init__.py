from . import classify
from . import compare
from . import gather
from . import index
from . import rankinfo
from . import summarize
import sys

subcommands = ['classify', 'compare', 'gather', 'index', 'rankinfo', 'summarize']
subcommandstr = ' -- '.join(sorted(subcommands))

def subparser(subparsers):
    subparser = subparsers.add_parser('lca')
    s = subparser.add_subparsers(
        title='Subcommands', dest='subcmd', metavar='subcmd', help=subcommandstr,
        description='Invoke "sourmash lca <subcmd> --help" for more details on executing each subcommand.'
    )
    for subcmd in subcommands:
        getattr(sys.modules[__name__], subcmd).subparser(s)
    subparser._action_groups.reverse()
    subparser._optionals.title = 'Options'
