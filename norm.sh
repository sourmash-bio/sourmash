for i in *.fastq.gz; do normalize-by-median.py -C 20 -M 1e9 -k 21 $i; done
