"""show details of signature"""

from argparse import FileType
import csv

import sourmash
from sourmash.logging import notify, print_results, error


def subparser(subparsers):
    subparser = subparsers.add_parser('describe')
    subparser.add_argument('signatures', nargs='+')
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '--csv', metavar='FILE', type=FileType('wt'),
        help='output information to a CSV file'
    )


def describe(signatures, quiet=True, csvout=None):
    siglist = []
    for sigfile in signatures:
        this_siglist = []
        try:
            this_siglist = sourmash.load_signatures(sigfile, quiet=quiet, do_raise=True)
            for k in this_siglist:
                siglist.append((k, sigfile))
        except Exception as exc:
            error('\nError while reading signatures from {}:'.format(sigfile))
            error(str(exc))
            error('(continuing)')

        notify('loaded {} signatures from {}...', len(siglist), sigfile,
               end='\r')

    notify('loaded {} signatures total.', len(siglist))

    w = None
    if csvout:
        w = csv.DictWriter(
            csvout, [
                'signature_file', 'md5', 'ksize', 'moltype', 'num', 'scaled',
                'n_hashes', 'seed', 'with_abundance','name', 'filename',
                'license'
            ],
            extrasaction='ignore'
        )
        w.writeheader()

    # extract info, write as appropriate.
    for (sig, signature_file) in siglist:
        mh = sig.minhash
        ksize = mh.ksize
        moltype = 'DNA'
        if mh.is_protein:
            if mh.dayhoff:
                moltype = 'dayhoff'
            elif mh.hp:
                moltype = 'hp'
            else:
                moltype = 'protein'
        scaled = mh.scaled
        num = mh.num
        seed = mh.seed
        n_hashes = len(mh)
        with_abundance = 0
        if mh.track_abundance:
            with_abundance = 1
        md5 = sig.md5sum()
        name = sig.name()
        filename = sig.filename
        license = sig.license

        if w:
            w.writerow(locals())

        print_results('''\
---
signature filename: {signature_file}
signature: {name}
source file: {filename}
md5: {md5}
k={ksize} molecule={moltype} num={num} scaled={scaled} seed={seed} track_abundance={with_abundance}
size: {n_hashes}
signature license: {license}
''', **locals())


def main(args):
    describe(args.signatures, quiet=args.quiet, csvout=args.csv)
