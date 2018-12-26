"""
Command-line entry point for 'python -m sourmash.sig'
"""
import sys
import argparse

from ..logging import set_quiet, error

usage='''
sourmash signature <command> [<args>] - manipulate/work with signature files.

** Commands can be:

merge <signature> [<signature> ...]

** Use '-h' to get subcommand-specific help, e.g.

sourmash signature merge -h
'''

def main(sysv_args):
    set_quiet(False)

    merge_cmd = lambda x: None
    commands = {'merge': merge_cmd}

    parser = argparse.ArgumentParser(
        description='signature file manipulation utilities', usage=usage)
    parser.add_argument('sig_command', nargs='?')
    args = parser.parse_args(sysv_args[0:1])

    if not args.sig_command:
        print(usage)
        sys.exit(1)

    if args.sig_command not in commands:
        error('Unrecognized command: {}', args.sig_command)
        parser.print_help()
        sys.exit(1)

    cmd = commands.get(args.sig_command)
    cmd(sysv_args[1:])

if __name__ == '__main__':
    main(sys.argv[1:])
