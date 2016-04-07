#! /usr/bin/env python
import sig
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('sigfiles', nargs='+')
    parser.add_argument('-f', '--force', action='store_true')
    parser.add_argument('-n', '--no-save', action='store_true')
    args = parser.parse_args()

    ignore_md5sum = False
    if args.force:
        print('ignoring md5sum because of --force')
        ignore_md5sum = True

    for filename in args.sigfiles:
        print('loading', filename)
        s = sig.load_signature(open(filename, 'r'), ignore_md5sum)

        data = s.save()

        if not args.no_save:
            print('loaded and saving', filename)
            fp = open(filename, 'w')
            fp.write(data)
            fp.close()
        
main()
