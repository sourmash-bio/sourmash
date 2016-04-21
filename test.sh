#! /bin/bash -exu
./sourmash compute test.fq -f
./sourmash compare test.fq.sig test.fq.sig
./sourmash clean -n -s test.fq.sig

