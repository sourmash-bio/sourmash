#! /usr/bin/env python3
import urllib.request
import csv


def main():
    acc_dict = {}
    r = csv.reader(open('sra_result.csv', 'r', encoding='utf-8'))
    for row in r:
        if row[0] == 'Experiment Accession':
            continue

        acc = row[0]
        organism = row[2]
#        if organism.startswith('Strongyl'):
#            continue
        acc_dict[acc] = organism.split(' ')[1]

    r = csv.reader(open('ftp_list.csv', 'r', encoding='utf-8'))
    for row in r:
        if row[0] in acc_dict:
            srr = row[1]
            print('../../sourmash clean --set-name %s *%s*' % (acc_dict[row[0]],
                                                           srr))
            print('mv *%s*.sig %s-%s.sig' % (srr, acc_dict[row[0]], srr))

if __name__ == '__main__':
    main()
