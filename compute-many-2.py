#! /usr/bin/env python
import sourmash
import argparse
import screed
from cPickle import load
import scipy
import pylab
import scipy.cluster.hierarchy as sch
from sklearn import metrics

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('dumpfile')
    args = parser.parse_args()

    labels = dict([ i.strip().split(' ', 2) for i in open('labels.txt') ])

    emins = load(open(args.dumpfile))
    estimators = []
    for (filename, mins) in emins:
        E = sourmash.Estimators()
        E._mins = mins
        estimators.append((filename, E))

    D = scipy.zeros([len(estimators), len(estimators)])
    
    i = 0
    for f, E in estimators:
        j = 0
        for f2, E2 in estimators:
            #print labels[f], labels[f2], E.jaccard(E2)
            D[i][j] = E.jaccard(E2)
            j += 1
            
        print '%20s\t%s' % (labels[f], D[i , :,],)
        i += 1


    fig = pylab.figure(figsize=(8,8))
    ax1 = fig.add_axes([0.09,0.1,0.2,0.6])
    Y = sch.linkage(D, method='single') # centroid
    Z1 = sch.dendrogram(Y, orientation='right')
    ax1.set_xticks([])
    ax1.set_yticks([])

    # Compute and plot second dendrogram.
    ax2 = fig.add_axes([0.3,0.71,0.6,0.2])
    Y = sch.linkage(D, method='single')
    Z2 = sch.dendrogram(Y)
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

    

if __name__ == '__main__':
    main()
    
