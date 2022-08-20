#! /usr/bin/env python
"""
Make plots using the distance matrix+labels output by ``sourmash compare``.
"""
from .logging import error, notify
try:
    import numpy
    import pylab
    import scipy.cluster.hierarchy as sch
except (RuntimeError, ImportError):
    pass

def load_matrix_and_labels(basefile):
    """Load the comparison matrix and associated labels.

    Returns a square numpy matrix & list of labels.
    """
    D = numpy.load(open(basefile, 'rb'))
    labeltext = [x.strip() for x in open(basefile + '.labels.txt')]
    return (D, labeltext)


def plot_composite_matrix(D, labeltext, show_labels=True, show_indices=True,
                          vmax=1.0, vmin=0.0, force=False):
    """Build a composite plot showing dendrogram + distance matrix/heatmap.

    Returns a matplotlib figure."""
    if D.max() > 1.0 or D.min() < 0.0:
        error('This matrix doesn\'t look like a distance matrix - min value {}, max value {}', D.min(), D.max())
        if not force:
            raise ValueError("not a distance matrix")
        else:
            notify('force is set; scaling to [0, 1]')
            D -= D.min()
            D /= D.max()

    if show_labels:
        show_indices = True

    fig = pylab.figure(figsize=(11, 8))
    ax1 = fig.add_axes([0.09, 0.1, 0.2, 0.6])

    # plot dendrogram
    Y = sch.linkage(D, method='single')  # centroid

    dendrolabels = labeltext
    if not show_labels:
        dendrolabels = [str(i) for i in range(len(labeltext))]

    Z1 = sch.dendrogram(Y, orientation='left', labels=dendrolabels,
                        no_labels=not show_indices, get_leaves=True)
    ax1.set_xticks([])

    xstart = 0.45
    width = 0.45
    if not show_labels:
        xstart = 0.315
    scale_xstart = xstart + width + 0.01

    # re-order labels along rows, top to bottom
    idx1 = Z1['leaves']
    reordered_labels = [ labeltext[i] for i in idx1 ]

    # reorder D by the clustering in the dendrogram
    D = D[idx1, :]
    D = D[:, idx1]

    # show matrix
    axmatrix = fig.add_axes([xstart, 0.1, width, 0.6])

    im = axmatrix.matshow(D, aspect='auto', origin='lower',
                          cmap=pylab.cm.YlGnBu, vmin=vmin, vmax=vmax)
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

    # Plot colorbar.
    axcolor = fig.add_axes([scale_xstart, 0.1, 0.02, 0.6])
    pylab.colorbar(im, cax=axcolor)

    return fig, reordered_labels, D
