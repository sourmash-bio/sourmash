#! /usr/bin/env python
import argparse
import os
import screed
import sourmash
import sourmash_signature as sig

LEN_CUTOFF=80

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('filenames', nargs='+')
    parser.add_argument('-k', '--ksizes',
                        default='31', help='comma-separated list')
    parser.add_argument('-f', '--force', action='store_true')
    args = parser.parse_args()

    ksizes = args.ksizes
    if ',' in ksizes:
        ksizes = ksizes.split(',')
        ksizes = list(map(int, ksizes))
    else:
        ksizes = [int(ksizes)]

    print('Computing signature for ksizes: %s' % str(ksizes))

    for filename in args.filenames:
        sigfile = os.path.basename(filename) + '.sig'
        if os.path.exists(sigfile) and not args.force:
            print('skipping', filename, '- already done')
            continue

        Elist = []
        for k in ksizes:
            E = sourmash.Estimators(ksize=k, n=500)
            Elist.append(E)
            
        for n, record in enumerate(screed.open(filename)):
            if n % 10000 == 0:
                print('...', filename, n)
            if 0 and n % 100000 == 0 and n:
                siglist = []
                for E in Elist:
                    signature = sig.SourmashSignature('titus@idyll.org',
                                              E,
                                              filename=filename)
                    siglist.append(signature)
            
                data = sig.save_signatures(siglist)
                fp = open(sigfile + '.%d' % n, 'w')
                fp.write(data)
                fp.close()

            s = record.sequence
            for i in 'R', 'W', 'Y':
                s = s.replace(i, 'N')
            for E in Elist:
                E.add_sequence(s)

        siglist = []
        for E in Elist:
            signature = sig.SourmashSignature('titus@idyll.org',
                                              E,
                                              filename=filename)
            siglist.append(signature)
            
        data = sig.save_signatures(siglist)
        fp = open(sigfile, 'w')
        fp.write(data)
        fp.close()
    

if __name__ == '__main__':
    main()
    
