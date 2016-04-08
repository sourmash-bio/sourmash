#! /usr/bin/env python
import sourmash
import argparse
import screed
import sig
import scipy
import pylab
import scipy.cluster.hierarchy as sch
from sklearn import metrics

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('signatures', nargs='+')
    args = parser.parse_args()

    siglist = []
    for filename in args.signatures:
        data = open(filename).read()
        print('loading', filename)
        signature = sig.load_signature(data)
        siglist.append(signature)

    D = scipy.zeros([len(siglist), len(siglist)])
    
    i = 0
    labeltext = []
    samples = []
    for E in siglist:
        j = 0
        for E2 in siglist:
            D[i][j] = E.jaccard(E2)
            j += 1
            
        print('%d-%20s\t%s' % (i, E.name(), D[i , :,],))
        labeltext.append('%d-%s' % (i,E.name()))
        samples.append((i, E.name(), E.d['filename']))
        i += 1

    fig = pylab.figure(figsize=(5,8))
    #ax1 = fig.add_axes([0.09,0.1,0.2,0.6])
    Y = sch.linkage(D, method='single') # centroid
    Z1 = sch.dendrogram(Y, labels=labeltext)
    #ax1.set_xticks([])
    #ax1.set_yticks([])

    fig.show()
    fig.savefig('xxx.png')

    # XXX

    fig = pylab.figure(figsize=(8,8))
    #ax1 = fig.add_axes([0.09,0.1,0.2,0.6])                                     
    Y = sch.linkage(D, method='single') # centroid                              
    Z1 = sch.dendrogram(Y, orientation='right', labels=labeltext)
    fig.show()
    fig.savefig('xxx.png')

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
    im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=pylab.cm.YlGnBu)
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

    # Plot colorbar.
    axcolor = fig.add_axes([0.91,0.1,0.02,0.6])
    pylab.colorbar(im, cax=axcolor)
    fig.show()
    fig.savefig('dendrogram.png')

    for i, name, filename in samples:
        print(i, '\t', name)

if __name__ == '__main__':
    main()
    
