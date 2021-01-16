import sourmash


def main(arglist=None):
    args = sourmash.cli.get_parser().parse_args(arglist)
    if hasattr(args, 'subcmd'):
        mod = getattr(sourmash.cli, args.cmd)
        submod = getattr(mod, args.subcmd)
        mainmethod = getattr(submod, 'main')
    else:
        mod = getattr(sourmash.cli, args.cmd)
        mainmethod = getattr(mod, 'main')
    return mainmethod(args)


if __name__ == '__main__':
    main()
