#! /usr/bin/env python
from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
import sourmash_fig
import argparse

import numpy
import scipy
import pylab
import scipy.cluster.hierarchy as sch
import os.path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('distances', help="output from 'sourmash compare'")
    parser.add_argument('--pdf', action='store_true')
    args = parser.parse_args()

    D_filename = args.distances
    labelfilename = D_filename + '.labels.txt'

    dendrogram_out = os.path.basename(D_filename) + '.dendro'
    if args.pdf:
        dendrogram_out += '.pdf'
    else:
        dendrogram_out += '.png'

    matrix_out = os.path.basename(D_filename) + '.matrix'
    if args.pdf:
        matrix_out += '.pdf'
    else:
        matrix_out += '.png'

    D = numpy.load(open(D_filename, 'rb'))
    labeltext = [ x.strip() for x in open(labelfilename) ]
    samples = [ (i, x) for i, x in enumerate(labeltext) ]

    fig = pylab.figure(figsize=(8,8))
    #ax1 = fig.add_axes([0.09,0.1,0.2,0.6])                                     
    Y = sch.linkage(D, method='single') # centroid                              
    Z1 = sch.dendrogram(Y, orientation='right', labels=labeltext)
    fig.show()
    fig.savefig(dendrogram_out)
    print('wrote', dendrogram_out)

    fig = sourmash_fig.plot_composite_matrix(D, labeltext, show_labels=False,
                                             show_indices=True)
    fig.savefig(matrix_out)
    print('wrote', matrix_out)

    for i, name in samples:
        print(i, '\t', name)

if __name__ == '__main__':
    main()
