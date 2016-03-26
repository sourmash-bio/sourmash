../nullgraph/make-random-genome.py --name a -s 1 -l 1000 > genome-a.fa
../nullgraph/make-random-genome.py --name b -s 2 -l 1000 > genome-b.fa
cat genome-a.fa genome-b.fa > genomes-both.fa
../nullgraph/make-reads.py -C 50 genome-a.fa > reads-a.fa
../nullgraph/make-reads.py -C 50 genome-b.fa > reads-b.fa
echo '>both' > genomes-both.fa
cat genome-a.fa genome-b.fa | grep -v ^'>' >> genomes-both.fa
../nullgraph/make-reads.py -C 50 genomes-both.fa > reads-both.fa
