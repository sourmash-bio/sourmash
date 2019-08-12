"""
Compute nearest neighbor graphs on SBTs for clustering and 2d embedding
"""

from itertools import groupby
import math

import numpy as np
import scipy.sparse
from umap.umap_ import smooth_knn_dist, compute_membership_strengths, spectral_layout

from .sbt import Leaf
from .sbtmh import SigLeaf



### --- Nearest Neighbor Functionality ---

def get_leaves(tree):
    for i, node in tree.nodes.items():
        if isinstance(node, SigLeaf) or isinstance(node, Leaf):
            yield i, node


def nearest_neighbors(tree, n_neighbors, ignore_abundance, downsample,
                      min_similarity=0.0):
    """Yield leaf1, leaf2, 1-similarity"""
    n_parent_levels = math.log2(n_neighbors) + 1

    for position1, leaf1 in get_leaves(tree):
        n = 1
        upper_internal_node = tree.parent(position1)
        while n < n_parent_levels:
            upper_internal_node = tree.parent(upper_internal_node.pos)
            n += 1
        #         print("upper_internal_node:", upper_internal_node)
        leaves = tree.get_leaves_under(upper_internal_node.pos)

        similarities = []
        for leaf2 in leaves:
            if leaf2 == leaf1:
                continue
            similarity = leaf1.data.similarity(leaf2.data,
                                              ignore_abundance=ignore_abundance,
                                              downsample=downsample)
            if similarity > min_similarity:
                similarities.append(
                    [leaf1.data.name(), leaf2.data.name(), similarity])
        adjacent = sorted(similarities, key=lambda x: x[1])[-n_neighbors:]
        for adjacency in adjacent:
            yield adjacency


def similarity_adjacency_to_knn(adjacencies, tree):
    """Convert adjacency list to k-nearest neighbor indices and distances

    "distance" = dissimilarity = 1-similarity

    :param adjacencies:
    :param tree:
    :return:
    """
    # Make arbitrary index integers for each leaf
    leaf_to_index = dict(
        (node.data.name(), i) for i, node in enumerate(get_leaves(tree)))
    knn_indices = []
    knn_dists = []

    for u, items in groupby(adjacencies, key=lambda x: x[0]):
        knn_indices_line = []
        knn_dists_line = []
        for u, v, similarity in items:
            knn_indices_line.append(leaf_to_index[v])

            # Dissimilarity = 1-similarity
            knn_dists_line.append(1 - similarity)
        knn_indices.append(knn_indices_line)
        knn_dists.append(knn_dists_line)
    knn_indices = np.array(knn_indices)
    knn_dists = np.array(knn_dists)

    return knn_indices, knn_dists, leaf_to_index


def fuzzy_simplicial_set(knn_indices, knn_dists, n_neighbors,
                         local_connectivity=1,
                         set_op_mix_ratio=1):
    """Compute probabalistic existence of edges based on distances

    Adapted from `fuzzy_simplicial_set` from umap-learn:
    https://github.com/lmcinnes/umap/blob/master//umap/umap_.py#L474

    :param knn_indices:
    :param knn_dists:
    :param n_neighbors:
    :param local_connectivity:
    :param set_op_mix_ratio:
    :return:
    """
    sigmas, rhos = smooth_knn_dist(
        knn_dists, n_neighbors, local_connectivity=local_connectivity
    )

    rows, cols, vals = compute_membership_strengths(
        knn_indices, knn_dists, sigmas, rhos
    )

    result = scipy.sparse.coo_matrix(
        (vals, (rows, cols)),
        shape=(knn_dists.shape[0], knn_dists.shape[0])
    )
    result.eliminate_zeros()

    transpose = result.transpose()

    prod_matrix = result.multiply(transpose)

    result = (
            set_op_mix_ratio * (result + transpose - prod_matrix)
            + (1.0 - set_op_mix_ratio) * prod_matrix
    )

    result.eliminate_zeros()
    return result


def make_search_graph(knn_indices, knn_dists):
    _search_graph = scipy.sparse.lil_matrix(
        (knn_indices.shape[0], knn_indices.shape[0]), dtype=np.int8
    )
    _search_graph.rows = knn_indices
    _search_graph.data = (knn_dists != 0).astype(np.int8)
    _search_graph = _search_graph.maximum(
        _search_graph.transpose()
    ).tocsr()
    return _search_graph



def simplicial_set_embedding(sbt, graph, n_components, random_state,
                            ignore_abundance, downsample,
                             init='spectral', n_epochs=0, ):
    graph = graph.tocoo()
    graph.sum_duplicates()
    n_vertices = graph.shape[1]

    if n_epochs <= 0:
        # For smaller datasets we can use more epochs
        if graph.shape[0] <= 10000:
            n_epochs = 500
        else:
            n_epochs = 200

    graph.data[graph.data < (graph.data.max() / float(n_epochs))] = 0.0
    graph.eliminate_zeros()

    if isinstance(init, str) and init == "random":
        embedding = random_state.uniform(
            low=-10.0, high=10.0, size=(graph.shape[0], n_components)
        ).astype(np.float32)
    elif isinstance(init, str) and init == "spectral":
        # We add a little noise to avoid local minima for optimization to come
        initialisation = spectral_layout(
            sbt,
            graph,
            n_components,
            random_state,
            ignore_abundance=ignore_abundance,
            downsample=downsample,
        )
        expansion = 10.0 / np.abs(initialisation).max()
        embedding = (initialisation * expansion).astype(
            np.float32
        ) + random_state.normal(
            scale=0.0001, size=[graph.shape[0], n_components]
        ).astype(
            np.float32
        )


def spectral_layout(sbt, graph, dim, random_state, ):
    n_samples = graph.shape[0]
    n_connected_components, labels = scipy.sparse.csgraph.connected_components(graph)
    n_connected_components



def multi_connected_component_layout(sbt, graph, n_components, labels, dim,
                                     random_state):
    if n_components > 2 * dim:
        meta_embedding = component_layout(
            data,
            n_components,
            component_labels,
            dim,
            metric=metric,
            metric_kwds=metric_kwds,
        )
    else:
        k = int(np.ceil(n_components / 2.0))
        base = np.hstack([np.eye(k), np.zeros((k, dim - k))])
        meta_embedding = np.vstack([base, -base])[:n_components]



def component_layout(sbt, n_connected_components, component_labels, dim,
                     ignore_abundance=False, downsample=False):
    """Provide a layout relating the separate connected components. This is done
    by taking the centroid of each component and then performing a spectral embedding
    of the centroids.
    Parameters
    ----------
    data: array of shape (n_samples, n_features)
        The source data -- required so we can generate centroids for each
        connected component of the graph.
    n_components: int
        The number of distinct components to be layed out.
    component_labels: array of shape (n_samples)
        For each vertex in the graph the label of the component to
        which the vertex belongs.
    dim: int
        The chosen embedding dimension.
    metric: string or callable (optional, default 'euclidean')
        The metric used to measure distances among the source data points.
    metric_kwds: dict (optional, default {})
        Keyword arguments to be passed to the metric function.
    Returns
    -------
    component_embedding: array of shape (n_components, dim)
        The ``dim``-dimensional embedding of the ``n_components``-many
        connected components.
    """

    component_centroids = np.empty((n_connected_components, data.shape[1]), dtype=np.float64)

    for label in range(n_components):
        component_centroids[label] = data[component_labels == label].mean(axis=0)

    distance_matrix = pairwise_distances(
        component_centroids, metric=metric, **metric_kwds
    )
    affinity_matrix = np.exp(-distance_matrix ** 2)

    component_embedding = SpectralEmbedding(
        n_components=dim, affinity="precomputed"
    ).fit_transform(affinity_matrix)
    component_embedding /= component_embedding.max()

    return component_embedding
