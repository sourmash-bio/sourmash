#! /usr/bin/env python
from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
from sourmash_lib import fig as sourmash_fig
import argparse

import numpy
import scipy
import pylab
import scipy.cluster.hierarchy as sch
import os.path

def main():
    # set up cmd line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('distances', help="output from 'sourmash compare'")
    parser.add_argument('--pdf', action='store_true')
    parser.add_argument('--labels', action='store_true')
    parser.add_argument('--indices', action='store_false')
    parser.add_argument('--vmax', default=1.0, type=float)
    parser.add_argument('--vmin', default=0.0, type=float)
    args = parser.parse_args()

    # load files
    D_filename = args.distances
    labelfilename = D_filename + '.labels.txt'

    D = numpy.load(open(D_filename, 'rb'))
    labeltext = [ x.strip() for x in open(labelfilename) ]

    # build filenames, decide on PDF/PNG output
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

    ### make the dendrogram:
    fig = pylab.figure(figsize=(8,5))
    ax1 = fig.add_axes([0.1, 0.1, 0.7, 0.8])
    ax1.set_xticks([])
    ax1.set_yticks([])
    
    Y = sch.linkage(D, method='single') # cluster!
    Z1 = sch.dendrogram(Y, orientation='right', labels=labeltext)
    fig.savefig(dendrogram_out)
    print('wrote', dendrogram_out)

    ### make the dendrogram+matrix:
    fig = sourmash_fig.plot_composite_matrix(D, labeltext,
                                             show_labels=args.labels,
                                             show_indices=args.indices,
                                             vmin=args.vmin,
                                             vmax=args.vmax)
    fig.savefig(matrix_out)
    print('wrote', matrix_out)

    # print out sample numbering for FYI.
    for i, name in enumerate(labeltext):
        print(i, '\t', name)

if __name__ == '__main__':
    main()
