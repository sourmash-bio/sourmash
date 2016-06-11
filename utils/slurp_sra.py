#! /usr/bin/env python3
import urllib.request
import csv

filereport_url="https://www.ebi.ac.uk/ena/data/warehouse/filereport"\
               "?accession=%s&result=read_run&"\
               "fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes"


def retrieve_fastq(accession):
    url = filereport_url % accession

    fp = urllib.request.urlopen(url)
    if fp.getcode() != 200:
        print('ERROR retrieving %s' % accession)
        return

    lines = fp.readlines()

    header = lines.pop(0)
    assert header == b'run_accession\tfastq_ftp\tfastq_md5\tfastq_bytes\n'

    for line in lines:
        if line.strip():
            row = line.strip().split(b'\t')
            if len(row) == 1:
                continue
            
            run_acc, ftp_list, md5_list, bytes_list = row
            if b';' in ftp_list:
                ftp_list = ftp_list.split(b';')
                md5_list = md5_list.split(b';')
                bytes_list = bytes_list.split(b';')

                for f,m,b in zip(ftp_list, md5_list, bytes_list):
                    yield run_acc, f, m, b
            else:
                yield run_acc, ftp_list, md5_list, bytes_list

def main():
    fp = open('ftp_list.csv', 'w', encoding='utf-8')
    w = csv.writer(fp, dialect=csv.excel)
    
    r = csv.reader(open('sra.csv', 'r', encoding='utf-8'))
    for row in r:
        if row[0] == 'Experiment Accession':
            continue
        for acc, url, md5, n_bytes in retrieve_fastq(row[0]):
            acc0 = bytes(row[0], 'utf-8')
            x = [acc0, acc, url, md5, n_bytes]
            x = [ i.decode('utf-8') for i in x ]
            print(x[0], x[2])
            w.writerow(x)

if __name__ == '__main__':
    main()
