#! /usr/bin/env python
import screed
import sys
import string

for record in screed.open(sys.argv[1]):
    s = record.sequence.upper()
    s = s.replace('R', 'N')
    s = s.replace('W', 'N')
    s = s.replace('Y', 'N')
    if s.count('A') + s.count('G') + s.count('T') + s.count('C') + s.count('N') != len(s):
        s = s.replace('A', '')
        s = s.replace('C', '')
        s = s.replace('G', '')
        s = s.replace('T', '')
        s = s.replace('N', '')
        print s
 #    print '>%s\n%s' % (record.name, record.sequence)
