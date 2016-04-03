#! /usr/bin/env python2
import argparse
import sys
import sourmash
import screed
import scipy
import pylab
import scipy.cluster.hierarchy as sch
from sklearn import metrics

LEN_SCALING=10

# asymmetric mh based on size?
# correlate w/trinity family info


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('filename')
    args = parser.parse_args()

    x = []
    for n, record in enumerate(screed.open(args.filename)):
        if n % 100 == 0:
            print >>sys.stderr, '...', n
            if n >= 5000:
                break
        num_sigs = int(len(record.sequence) / LEN_SCALING)
        E = sourmash.Estimators(n=num_sigs, ksize=12)
        E.add_sequence(record.sequence)
        x.append(E)

    dists = []
    

    D = scipy.zeros([len(x), len(x)])
    for i in range(len(x)):
        for j in range(len(x)):
            v = x[i].jaccard(x[j])
            if v >= 0.1:
                dists.append(v)
            D[i, j] = v
            #print i, j, v

    pylab.clf()
    pylab.hist(dists, 100)
    pylab.savefig('hist.png')

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

    # break on 5% or greater similarity
    clusters = sch.fcluster(Y, 0.95, criterion='distance')
    print len(set(clusters))

    d = {}
    for k in clusters:
        d[k] = d.get(k, 0) + 1
    for k, v in sorted(d.items(), key=lambda x: -x[1]):
        if v == 1:
            break
        print k, v

    names_d = {}
    names_clusterid = {}
    cluster_n = 1
    names_cluster_assign = []

    fp = open('xxx.fa', 'w')
    for n, record in enumerate(screed.open(args.filename)):
        if n >= len(clusters):
            break
        print >>fp, '>%s %d\n%s' % (record.name.split()[0], clusters[n], record.sequence)

        k = record.name.split()[0].split('_')[1]
        c = clusters[n]

        x = names_d.get(k, [])
        if c not in x:
            x.append(c)
        names_d[k] = x

        if k in names_clusterid:
            names_cluster_assign.append(names_clusterid[k])
        else:
            names_clusterid[k] = cluster_n
            names_cluster_assign.append(cluster_n)
            cluster_n += 1

    for k in names_d:
        if len(names_d[k]) > 1:
            print k, names_d[k]

    ###

    print len(set(clusters))
    print len(set(names_cluster_assign))
    
    print metrics.adjusted_mutual_info_score(clusters, clusters)
    print metrics.adjusted_mutual_info_score(clusters, names_cluster_assign)


if __name__ == '__main__':
    main()
