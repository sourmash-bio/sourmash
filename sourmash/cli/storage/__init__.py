import sys

from . import convert

subcommands = ['convert']
subcommandstr = ' -- '.join(sorted(subcommands))


def subparser(subparsers):
    subparser = subparsers.add_parser('storage')
    s = subparser.add_subparsers(
        title='Subcommands', dest='subcmd', metavar='subcmd', help=subcommandstr,
        description='Invoke "sourmash storage <subcmd> --help" for more details on executing each subcommand.'
    )
    for subcmd in subcommands:
        getattr(sys.modules[__name__], subcmd).subparser(s)
    subparser._action_groups.reverse()
    subparser._optionals.title = 'Options'
