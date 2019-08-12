"""
Compute nearest neighbor graphs on SBTs for clustering and 2d embedding
"""
import copy
from itertools import groupby
import math
from warnings import warn

import numpy as np
import scipy.sparse
from sklearn.manifold import SpectralEmbedding
from umap.umap_ import smooth_knn_dist, compute_membership_strengths, spectral_layout

from .compare import compare_all_pairs
from .sbt import Leaf
from .sbtmh import SigLeaf



### --- Nearest Neighbor Functionality ---

def get_leaves(tree):
    for i, node in tree.nodes.items():
        if isinstance(node, (Leaf, SigLeaf)):
            yield i, node


def nearest_neighbors(tree, n_neighbors, ignore_abundance, downsample,
                      min_similarity=0.0):
    """Yield leaf1, leaf2, 1-similarity"""
    n_parent_levels = math.log2(n_neighbors) + 1

    adjacencies = []

    for position1, leaf1 in get_leaves(tree):
        n = 1
        upper_internal_node = tree.parent(position1)
        while n < n_parent_levels:
            upper_internal_node = tree.parent(upper_internal_node.pos)
            n += 1
        #         print("upper_internal_node:", upper_internal_node)
        leaves = tree.leaves_under(upper_internal_node.pos)

        similarities = []
        for leaf2 in leaves:
            if leaf2 == leaf1:
                continue
            similarity = leaf1.data.similarity(leaf2.data,
                                              ignore_abundance=ignore_abundance,
                                              downsample=downsample)
            # if similarity > min_similarity:
            similarities.append(
                [leaf1.data.name(), leaf2.data.name(), similarity])
        adjacent = sorted(similarities, key=lambda x: x[1])[-n_neighbors:]
        adjacencies.extend(adjacent)
    return adjacencies


def similarity_adjacency_to_knn(adjacencies, tree):
    """Convert adjacency list to k-nearest neighbor indices and distances

    "distance" = dissimilarity = 1-similarity

    :param adjacencies:
    :param tree:
    :return:
    """
    # Make arbitrary index integers for each leaf
    leaf_to_position = dict(
        (node.data.name(), position) for position, node in tree.nodes.items()
        if isinstance(node, Leaf))
    leaf_to_index = dict((name, i) for i, name in enumerate(leaf_to_position.keys()))
    if len(leaf_to_index) > 0:
        raise ValueError("leaf_to_index dictionary not properly constructed! "
                         "Length is zero :(")
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
    """Create graph for searching (??)

    From UMAP
    https://github.com/lmcinnes/umap/blob/834184f9c0455f26db13ab148c0abd2d3767d968//umap/umap_.py#L1451

    :param knn_indices:
    :param knn_dists:
    :return:
    """
    _search_graph = scipy.sparse.lil_matrix(
        (knn_indices.shape[0], knn_indices.shape[0]), dtype=np.int8
    )
    _search_graph.rows = knn_indices
    _search_graph.data = (knn_dists != 0).astype(np.int8)
    _search_graph = _search_graph.maximum(
        _search_graph.transpose()
    ).tocsr()
    return _search_graph


def _clean_graph(graph, n_epochs):
    graph = graph.tocoo()
    graph.sum_duplicates()

    if n_epochs <= 0:
        # For smaller datasets we can use more epochs
        if graph.shape[0] <= 10000:
            n_epochs = 500
        else:
            n_epochs = 200

    graph.data[graph.data < (graph.data.max() / float(n_epochs))] = 0.0
    graph.eliminate_zeros()
    return graph, n_epochs

#
# def simplicial_set_embedding(sbt, graph, n_components, random_state,
#                             ignore_abundance, downsample,
#                              init='spectral', n_epochs=0, ):
#     graph, n_epochs = _clean_graph(graph, n_epochs)
#
#     if isinstance(init, str) and init == "random":
#         embedding = random_state.uniform(
#             low=-10.0, high=10.0, size=(graph.shape[0], n_components)
#         ).astype(np.float32)
#     elif isinstance(init, str) and init == "spectral":
#         # We add a little noise to avoid local minima for optimization to come
#         initialisation = spectral_layout(
#             sbt,
#             graph,
#             n_components,
#             random_state,
#             ignore_abundance=ignore_abundance,
#             downsample=downsample,
#         )
#         expansion = 10.0 / np.abs(initialisation).max()
#         embedding = (initialisation * expansion).astype(
#             np.float32
#         ) + random_state.normal(
#             scale=0.0001, size=[graph.shape[0], n_components]
#         ).astype(
#             np.float32
#         )
#
#
# def spectral_layout(sbt, graph, dim, random_state, ):
#     n_samples = graph.shape[0]
#     n_connected_components, labels = scipy.sparse.csgraph.connected_components(graph)
#     n_connected_components
#
#
#
#
#
def multi_connected_component_layout(tree, graph, n_connected_components,
                                     component_labels, dim, random_state,
                                     ignore_abundance, downsample):
    result = np.empty((graph.shape[0], dim), dtype=np.float32)
    if n_connected_components > 2 * dim:
        meta_embedding = component_layout(
            tree,
            n_connected_components,
            component_labels,
            dim,
            ignore_abundance,
            downsample,
        )
    else:
        k = int(np.ceil(n_connected_components / 2.0))
        base = np.hstack([np.eye(k), np.zeros((k, dim - k))])
        meta_embedding = np.vstack([base, -base])[:n_connected_components]

    for label in range(n_connected_components):
        component_graph = graph.tocsr()[component_labels == label, :].tocsc()
        component_graph = component_graph[:, component_labels == label].tocoo()

        distances = pairwise_distances([meta_embedding[label]], meta_embedding)
        data_range = distances[distances > 0.0].min() / 2.0

        if component_graph.shape[0] < 2 * dim:
            result[component_labels == label] = (
                    random_state.uniform(
                        low=-data_range,
                        high=data_range,
                        size=(component_graph.shape[0], dim),
                    )
                    + meta_embedding[label]
            )
            continue

        diag_data = np.asarray(component_graph.sum(axis=0))
        # standard Laplacian
        # D = scipy.sparse.spdiags(diag_data, 0, graph.shape[0], graph.shape[0])
        # L = D - graph
        # Normalized Laplacian
        I = scipy.sparse.identity(component_graph.shape[0], dtype=np.float64)
        D = scipy.sparse.spdiags(
            1.0 / np.sqrt(diag_data),
            0,
            component_graph.shape[0],
            component_graph.shape[0],
        )
        L = I - D * component_graph * D

        k = dim + 1
        num_lanczos_vectors = max(2 * k + 1,
                                  int(np.sqrt(component_graph.shape[0])))
        try:
            eigenvalues, eigenvectors = scipy.sparse.linalg.eigsh(
                L,
                k,
                which="SM",
                ncv=num_lanczos_vectors,
                tol=1e-4,
                v0=np.ones(L.shape[0]),
                maxiter=graph.shape[0] * 5,
            )
            order = np.argsort(eigenvalues)[1:k]
            component_embedding = eigenvectors[:, order]
            expansion = data_range / np.max(np.abs(component_embedding))
            component_embedding *= expansion
            result[component_labels == label] = (
                    component_embedding + meta_embedding[label]
            )
        except scipy.sparse.linalg.ArpackError:
            warn(
                "WARNING: spectral initialisation failed! The eigenvector solver\n"
                "failed. This is likely due to too small an eigengap. Consider\n"
                "adding some noise or jitter to your data.\n\n"
                "Falling back to random initialisation!"
            )
            result[component_labels == label] = (
                    random_state.uniform(
                        low=-data_range,
                        high=data_range,
                        size=(component_graph.shape[0], dim),
                    )
                    + meta_embedding[label]
            )

    return result

def component_layout(tree, n_connected_components, component_labels, dim,
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

    centroid_signatures = []

    leaf_to_index = dict((node.data.name(), i) for i, (position, node) in
                         enumerate(get_leaves(tree)))
    leaf_to_position = dict(
        (node.data.name(), position) for position, node in get_leaves(tree))
    index_to_leaf = dict(zip(leaf_to_index.values(), leaf_to_index.keys()))

    for label in range(n_connected_components):
        graph_indexes = np.where(component_labels == label)[0]
        signatures = [tree.nodes.get(leaf_to_position[index_to_leaf[i]]).data
                      for i in graph_indexes]

        first_sig = signatures[0]
        merged_minhashes = copy.copy(first_sig.minhash)
        for sig in signatures[1:]:
            merged_minhashes = merged_minhashes.merge(sig.minhash)
        centroid_signatures.append(merged_minhashes)

    affinity_matrix = compare_all_pairs(centroid_signatures,
                                        ignore_abundance=ignore_abundance,
                                        downsample=downsample)

    component_embedding = SpectralEmbedding(
        n_components=dim, affinity="precomputed"
    ).fit_transform(affinity_matrix)
    component_embedding /= component_embedding.max()

    return component_embedding
