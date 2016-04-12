#! /usr/bin/env python
from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
import argparse
import numpy
import scipy
import pylab
import scipy.cluster.hierarchy as sch
from sklearn import metrics
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

    fig = pylab.figure(figsize=(8,8))
    ax1 = fig.add_axes([0.09,0.1,0.2,0.6])
    Y = sch.linkage(D, method='single') # centroid
    Z1 = sch.dendrogram(Y, orientation='right', labels=labeltext)
    ax1.set_xticks([])
    ax1.set_yticks([])

    # Compute and plot second dendrogram.
    ax2 = fig.add_axes([0.3,0.71,0.6,0.2])
    Y = sch.linkage(D, method='single')
    Z2 = sch.dendrogram(Y, labels=labeltext)
    ax2.set_xticks([])
    ax2.set_yticks([])

    # Plot distance matrix.
    axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])
    idx1 = Z1['leaves']
    idx2 = Z2['leaves']
    D = D[idx1,:]
    D = D[:,idx2]
    im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=pylab.cm.YlGnBu, vmin=0)
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

    # Plot colorbar.
    axcolor = fig.add_axes([0.91,0.1,0.02,0.6])
    pylab.colorbar(im, cax=axcolor)
    fig.show()
    fig.savefig(matrix_out)
    print('wrote', matrix_out)

    for i, name in samples:
        print(i, '\t', name)

if __name__ == '__main__':
    main()
