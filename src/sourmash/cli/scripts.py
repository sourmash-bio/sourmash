"""CLI plugin @CTB"""

usage="""

CLI plugins @CTB

"""
import argparse

def subparser(subparsers):
    subparser = subparsers.add_parser('scripts', description=__doc__, usage=usage)


def main(args):
    from sourmash.logging import debug_literal

    # this is what 'sourmash scripts' runs w/o any subcmd.
    debug_literal(f'cli.scripts: running {args}')
    if getattr(args, 'func', None):
        args.func(args)
    else:
        print('(default scripts messages goes here.)')
