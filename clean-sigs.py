#! /usr/bin/env python
import sig
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('sigfiles', nargs='+')
    parser.add_argument('-f', '--force', action='store_true',
                        help='ignore MD5 and other optional checks')
    parser.add_argument('-n', '--no-save', action='store_true')
    parser.add_argument('-s', '--show', action='store_true',
                        help='print to stdout instead of overwriting')
    parser.add_argument('--name-from', help='transfer name from here')
    args = parser.parse_args()

    ignore_md5sum = False
    if args.force:
        print('ignoring md5sum because of --force')
        ignore_md5sum = True

    for filename in args.sigfiles:
        print('loading', filename)
        siglist = sig.load_signatures(open(filename, 'r'), ignore_md5sum)

        for s in siglist:
            if not s.d.get('name'):
                print('warning, no name for:', filename)

        if args.name_from:              # hackity hack hack
            assert len(siglist) == 1
            other = sig.load_signatures(open(args.name_from, 'r'))
            s = siglist[0]
            s.d['name'] = other[0].d['name']

        data = sig.save_signatures(siglist)

        if not args.no_save:
            print('loaded and saving', filename)
            if args.show:
                print(data)
            else:
                fp = open(filename, 'w')
                fp.write(data)
                fp.close()
        
main()
