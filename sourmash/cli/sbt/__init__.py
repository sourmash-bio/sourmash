from . import combine
import sys

subcommands = ['combine']
subcommandstr = ' -- '.join(sorted(subcommands))

def subparser(subparsers):
    subparser = subparsers.add_parser('sbt')
    s = subparser.add_subparsers(
        title='Subcommands', dest='subcmd', metavar='subcmd', help=subcommandstr,
        description='Invoke "sourmash sbt <subcmd> --help" for more details on executing each subcommand.'
    )
    for subcmd in subcommands:
        getattr(sys.modules[__name__], subcmd).subparser(s)
    subparser._action_groups.reverse()
    subparser._optionals.title = 'Options'
