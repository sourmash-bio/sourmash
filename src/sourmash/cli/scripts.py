"""CLI plugin @CTB"""

# CTB TODO:
# * evaluate how multiple commands in a single package work
# * evaluate/check debugging capability
# * provide suggestions for URLs etc.
# * provide a _sourmash_ mechanism to list _all_ plugins (info?
#   or scripts list?)
# * provide guidance on how to test your CLI plugin at the CLI
#   (sourmash scripts, look for description etc.)

usage="""

CLI plugins @CTB

"""
import argparse
import sourmash

# decorate this module with the various extension objects
# e.g. 'sourmash scripts foo' will look up attribute 'scripts.foo'
# and we will return the extension class object. This is loaded below
# by sourmash.plugins.add_cli_scripts.
_extension_dict = {}
def __getattr__(name):
    if name in _extension_dict:
        return _extension_dict[name]
    raise AttributeError(name)

def subparser(subparsers):
    subparser = subparsers.add_parser('scripts', description=__doc__, usage=usage)

    # get individual help strings:
    descrs = list(sourmash.plugins.get_cli_scripts_descriptions())
    s = subparser.add_subparsers(title="extension commands",
                                 dest='subcmd',
                                 metavar='subcmd',
                                 help=argparse.SUPPRESS,
                                 description="\n".join(descrs))

    _extension_dict.update(sourmash.plugins.add_cli_scripts(s))

    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    # @CTB hmm, this doesn't work...
    subparser.add_argument(
        '-d', '--debug', action='store_true',
        help='provide debugging output'
    )


def main(args):
    raise Exception
    from sourmash.logging import set_quiet
    from sourmash.logging import debug_literal

    args.debug = True

    set_quiet(args.quiet, args.debug)

    # this is what 'sourmash scripts' runs if no subcommand is specified.
    # TODO -
    #   list plugins and descriptions?
    #   provide URLs? or suggest that on -h?
    debug_literal(f'cli.scripts: running {args}')

    if getattr(args, 'func', None):
        args.func(args)
    else:
        print('(default scripts messages goes here.)')
