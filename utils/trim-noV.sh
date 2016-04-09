for i in *.fastq.gz; do
trim-low-abund.py -C 4 -Z 20 -M 1e9 -k 21 $i --gzip -o $i.trimnoV.fq.gz;
done
