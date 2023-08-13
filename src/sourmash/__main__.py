"""
The main entry point for sourmash, defined in pyproject.toml.

Can also be executed as 'python -m sourmash'.
"""


def main(arglist=None):
    import sourmash
    args = sourmash.cli.parse_args(arglist)
    if hasattr(args, 'subcmd'):
        mod = getattr(sourmash.cli, args.cmd)
        submod = getattr(mod, args.subcmd)
        mainmethod = getattr(submod, 'main')
    else:
        mod = getattr(sourmash.cli, args.cmd)
        mainmethod = getattr(mod, 'main')

    retval = mainmethod(args)
    raise SystemExit(retval)


if __name__ == '__main__':
    main()
