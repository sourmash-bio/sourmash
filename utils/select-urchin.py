#! /usr/bin/env python3
import urllib.request
import csv


def main():
    acc_list = []
    r = csv.reader(open('sra_result.csv', 'r', encoding='utf-8'))
    for row in r:
        if row[0] == 'Experiment Accession':
            continue

        organism = row[2]
        if not organism.startswith('Strongyl'):
            continue
        acc_list.append(row[0])

    r = csv.reader(open('ftp_list.csv', 'r', encoding='utf-8'))
    for row in r:
        if row[0] in acc_list:
            url = row[2]
            if not url.endswith('_2.fastq.gz'):
                print("ftp://" + url)

if __name__ == '__main__':
    main()
