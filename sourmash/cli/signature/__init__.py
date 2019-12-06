from . import describe
from . import downsample
from . import extract
from . import flatten
from . import intersect
from . import merge
from . import rename
from . import subtract
from . import ingest
from . import export
from . import overlap
import sys

subcommands = [
    'describe', 'downsample', 'extract', 'flatten', 'intersect', 'merge',
    'rename', 'subtract', 'ingest', 'export', 'overlap'
]
subcommandstr = ' -- '.join(sorted(subcommands))

def subparser(subparsers):
    subparser = subparsers.add_parser('signature')
    s = subparser.add_subparsers(
        title='Subcommands', dest='subcmd', metavar='subcmd', help=subcommandstr,
        description='Invoke "sourmash signature <subcmd> --help" for more details on executing each subcommand.'
    )
    for subcmd in subcommands:
        getattr(sys.modules[__name__], subcmd).subparser(s)
    subparser._action_groups.reverse()
    subparser._optionals.title = 'Options'
