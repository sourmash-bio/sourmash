from __future__ import print_function
from itertools import groupby
from warnings import warn

from sklearn.base import BaseEstimator
from sklearn.utils import check_random_state, check_array
from sklearn.manifold import SpectralEmbedding
from sklearn.metrics import pairwise_distances
from sklearn.neighbors import KDTree
from sklearn.preprocessing import normalize

try:
    import joblib
except ImportError:
    # sklearn.externals.joblib is deprecated in 0.21, will be removed in 0.23
    from sklearn.externals import joblib

import numpy as np
import scipy.sparse
import scipy.sparse.csgraph

import umap.distances as dist


from umap.utils import deheap_sort, submatrix, ts
from umap.nndescent import (
    make_initialisations,
    make_initialized_nnd_search,
    initialise_search,
)
from umap.umap_ import (
    smooth_knn_dist,
    compute_membership_strengths,
    make_epochs_per_sample,
    optimize_layout,
    find_ab_params,
    general_simplicial_set_intersection,
    reset_local_connectivity,
    init_transform,
    categorical_simplicial_set_intersection,
    fuzzy_simplicial_set
)

from .compare import compare_all_pairs
from .logging import notify
from .signature import SourmashSignature
from .sbt import Leaf

import locale

locale.setlocale(locale.LC_NUMERIC, "C")

INT32_MIN = np.iinfo(np.int32).min + 1
INT32_MAX = np.iinfo(np.int32).max - 1

SMOOTH_K_TOLERANCE = 1e-5
MIN_K_DIST_SCALE = 1e-3
NPY_INFINITY = np.inf



def fuzzy_simplicial_set_from_knn(knn_indices, knn_dists, n_neighbors,
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


def simplicial_set_embedding(
    tree,
    graph,
    n_components,
    initial_alpha,
    a,
    b,
    gamma,
    negative_sample_rate,
    n_epochs,
    init,
    random_state,
    ignore_abundance,
    downsample,
    verbose,
):
    """Perform a fuzzy simplicial set embedding, using a specified
    initialisation method and then minimizing the fuzzy set cross entropy
    between the 1-skeletons of the high and low dimensional fuzzy simplicial
    sets.
    Parameters
    ----------
    data: array of shape (n_samples, n_features)
        The source data to be embedded by UMAP.
    graph: sparse matrix
        The 1-skeleton of the high dimensional fuzzy simplicial set as
        represented by a graph for which we require a sparse matrix for the
        (weighted) adjacency matrix.
    n_components: int
        The dimensionality of the euclidean space into which to embed the data.
    initial_alpha: float
        Initial learning rate for the SGD.
    a: float
        Parameter of differentiable approximation of right adjoint functor
    b: float
        Parameter of differentiable approximation of right adjoint functor
    gamma: float
        Weight to apply to negative samples.
    negative_sample_rate: int (optional, default 5)
        The number of negative samples to select per positive sample
        in the optimization process. Increasing this value will result
        in greater repulsive force being applied, greater optimization
        cost, but slightly more accuracy.
    n_epochs: int (optional, default 0)
        The number of training epochs to be used in optimizing the
        low dimensional embedding. Larger values result in more accurate
        embeddings. If 0 is specified a value will be selected based on
        the size of the input dataset (200 for large datasets, 500 for small).
    init: string
        How to initialize the low dimensional embedding. Options are:
            * 'spectral': use a spectral embedding of the fuzzy 1-skeleton
            * 'random': assign initial embedding positions at random.
            * A numpy array of initial embedding positions.
    random_state: numpy RandomState or equivalent
        A state capable being used as a numpy random state.
    metric: string
        The metric used to measure distance in high dimensional space; used if
        multiple connected components need to be layed out.
    metric_kwds: dict
        Key word arguments to be passed to the metric function; used if
        multiple connected components need to be layed out.
    verbose: bool (optional, default False)
        Whether to report information on the current progress of the algorithm.
    Returns
    -------
    embedding: array of shape (n_samples, n_components)
        The optimized of ``graph`` into an ``n_components`` dimensional
        euclidean space.
    """
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
            tree,
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
    else:
        init_data = np.array(init)
        if len(init_data.shape) == 2:
            if np.unique(init_data, axis=0).shape[0] < init_data.shape[0]:
                tree = KDTree(init_data)
                dist, ind = tree.query(init_data, k=2)
                nndist = np.mean(dist[:, 1])
                embedding = init_data + random_state.normal(
                    scale=0.001 * nndist, size=init_data.shape
                ).astype(np.float32)
            else:
                embedding = init_data

    epochs_per_sample = make_epochs_per_sample(graph.data, n_epochs)

    head = graph.row
    tail = graph.col

    rng_state = random_state.randint(INT32_MIN, INT32_MAX, 3).astype(np.int64)
    embedding = optimize_layout(
        embedding,
        embedding,
        head,
        tail,
        n_epochs,
        n_vertices,
        epochs_per_sample,
        a,
        b,
        rng_state,
        gamma,
        initial_alpha,
        negative_sample_rate,
        verbose=verbose,
    )

    return embedding



def spectral_layout(tree, graph, dim, random_state, ignore_abundance,
                    downsample):
    """Given a graph compute the spectral embedding of the graph. This is
    simply the eigenvectors of the laplacian of the graph. Here we use the
    normalized laplacian.

    Adapted from UMAP
    https://github.com/lmcinnes/umap/blob/master/umap/spectral.py#L199

    Parameters
    ----------
    data: array of shape (n_samples, n_features)
        The source data
    graph: sparse matrix
        The (weighted) adjacency matrix of the graph as a sparse matrix.
    dim: int
        The dimension of the space into which to embed.
    random_state: numpy RandomState or equivalent
        A state capable being used as a numpy random state.
    Returns
    -------
    embedding: array of shape (n_vertices, dim)
        The spectral embedding of the graph.
    """
    n_samples = graph.shape[0]
    n_connected_components, labels = scipy.sparse.csgraph.connected_components(graph)

    if n_connected_components > 1:
        warn(
            "Embedding a total of {} separate connected components using meta-embedding (experimental)".format(
                n_connected_components
            )
        )
        return multi_connected_component_layout(
            tree,
            graph,
            n_connected_components,
            labels,
            dim,
            random_state,
            ignore_abundance=ignore_abundance,
            downsample=downsample
        )

    diag_data = np.asarray(graph.sum(axis=0))
    # standard Laplacian
    # D = scipy.sparse.spdiags(diag_data, 0, graph.shape[0], graph.shape[0])
    # L = D - graph
    # Normalized Laplacian
    I = scipy.sparse.identity(graph.shape[0], dtype=np.float64)
    D = scipy.sparse.spdiags(
        1.0 / np.sqrt(diag_data), 0, graph.shape[0], graph.shape[0]
    )
    L = I - D * graph * D

    k = dim + 1
    num_lanczos_vectors = max(2 * k + 1, int(np.sqrt(graph.shape[0])))
    try:
        if L.shape[0] < 2000000:
            eigenvalues, eigenvectors = scipy.sparse.linalg.eigsh(
                L,
                k,
                which="SM",
                ncv=num_lanczos_vectors,
                tol=1e-4,
                v0=np.ones(L.shape[0]),
                maxiter=graph.shape[0] * 5,
            )
        else:
            eigenvalues, eigenvectors = scipy.sparse.linalg.lobpcg(
                L, random_state.normal(size=(L.shape[0], k)), largest=False, tol=1e-8
            )
        order = np.argsort(eigenvalues)[1:k]
        return eigenvectors[:, order]
    except scipy.sparse.linalg.ArpackError:
        warn(
            "WARNING: spectral initialisation failed! The eigenvector solver\n"
            "failed. This is likely due to too small an eigengap. Consider\n"
            "adding some noise or jitter to your data.\n\n"
            "Falling back to random initialisation!"
        )
        return random_state.uniform(low=-10.0, high=10.0, size=(graph.shape[0], dim))


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
                "WARNING: spectral initialisation failed! The eigenvector "
                "solver\nfailed. This is likely due to too small an eigengap."
                " Consider\nadding some noise or jitter to your data.\n\n"
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
    notify("Merging signatures in connected components", end='\r')
    centroid_signatures = merge_connected_component_signatures(
        tree, n_connected_components, component_labels)

    notify("Computing similarity matrix", end='\r')
    affinity_matrix = compare_all_pairs(centroid_signatures,
                                        ignore_abundance=ignore_abundance,
                                        downsample=downsample)

    component_embedding = SpectralEmbedding(
        n_components=dim, affinity="precomputed"
    ).fit_transform(affinity_matrix)
    component_embedding /= component_embedding.max()

    return component_embedding


def merge_connected_component_signatures(tree, n_connected_components,
                                         component_labels):
    centroid_signatures = [None for _ in range(n_connected_components)]
    leaf_to_index = dict((node.data.name(), i) for i, (position, node) in
                         enumerate(tree._leaves()))
    leaf_to_position = dict(
        (node.data.name(), position) for position, node in tree._leaves())
    index_to_leaf = dict(zip(leaf_to_index.values(), leaf_to_index.keys()))

    for label in range(n_connected_components):
        graph_indexes = np.where(component_labels == label)[0]
        signatures = [tree.nodes.get(leaf_to_position[index_to_leaf[i]]).data
                      for i in graph_indexes]

        first_sig = signatures[0]
        merged_signatures = SourmashSignature(first_sig.minhash, name=label)
        for s in signatures[1:]:
            merged_signatures.minhash = merged_signatures.minhash.merge(
                s.minhash)
        centroid_signatures[label] = merged_signatures
    return centroid_signatures

class UMAP(BaseEstimator):
    """Uniform Manifold Approximation and Projection
    Finds a low dimensional embedding of the data that approximates
    an underlying manifold.
    Parameters
    ----------
    n_neighbors: float (optional, default 15)
        The size of local neighborhood (in terms of number of neighboring
        sample points) used for manifold approximation. Larger values
        result in more global views of the manifold, while smaller
        values result in more local data being preserved. In general
        values should be in the range 2 to 100.
    n_components: int (optional, default 2)
        The dimension of the space to embed into. This defaults to 2 to
        provide easy visualization, but can reasonably be set to any
        integer value in the range 2 to 100.
    n_epochs: int (optional, default None)
        The number of training epochs to be used in optimizing the
        low dimensional embedding. Larger values result in more accurate
        embeddings. If None is specified a value will be selected based on
        the size of the input dataset (200 for large datasets, 500 for small).
    learning_rate: float (optional, default 1.0)
        The initial learning rate for the embedding optimization.
    init: string (optional, default 'spectral')
        How to initialize the low dimensional embedding. Options are:
            * 'spectral': use a spectral embedding of the fuzzy 1-skeleton
            * 'random': assign initial embedding positions at random.
            * A numpy array of initial embedding positions.
    min_dist: float (optional, default 0.1)
        The effective minimum distance between embedded points. Smaller values
        will result in a more clustered/clumped embedding where nearby points
        on the manifold are drawn closer together, while larger values will
        result on a more even dispersal of points. The value should be set
        relative to the ``spread`` value, which determines the scale at which
        embedded points will be spread out.
    spread: float (optional, default 1.0)
        The effective scale of embedded points. In combination with ``min_dist``
        this determines how clustered/clumped the embedded points are.
    set_op_mix_ratio: float (optional, default 1.0)
        Interpolate between (fuzzy) union and intersection as the set operation
        used to combine local fuzzy simplicial sets to obtain a global fuzzy
        simplicial sets. Both fuzzy set operations use the product t-norm.
        The value of this parameter should be between 0.0 and 1.0; a value of
        1.0 will use a pure fuzzy union, while 0.0 will use a pure fuzzy
        intersection.
    local_connectivity: int (optional, default 1)
        The local connectivity required -- i.e. the number of nearest
        neighbors that should be assumed to be connected at a local level.
        The higher this value the more connected the manifold becomes
        locally. In practice this should be not more than the local intrinsic
        dimension of the manifold.
    repulsion_strength: float (optional, default 1.0)
        Weighting applied to negative samples in low dimensional embedding
        optimization. Values higher than one will result in greater weight
        being given to negative samples.
    negative_sample_rate: int (optional, default 5)
        The number of negative samples to select per positive sample
        in the optimization process. Increasing this value will result
        in greater repulsive force being applied, greater optimization
        cost, but slightly more accuracy.
    transform_queue_size: float (optional, default 4.0)
        For transform operations (embedding new points using a trained model_
        this will control how aggressively to search for nearest neighbors.
        Larger values will result in slower performance but more accurate
        nearest neighbor evaluation.
    a: float (optional, default None)
        More specific parameters controlling the embedding. If None these
        values are set automatically as determined by ``min_dist`` and
        ``spread``.
    b: float (optional, default None)
        More specific parameters controlling the embedding. If None these
        values are set automatically as determined by ``min_dist`` and
        ``spread``.
    random_state: int, RandomState instance or None, optional (default: None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.
    metric_kwds: dict (optional, default None)
        Arguments to pass on to the metric, such as the ``p`` value for
        Minkowski distance. If None then no arguments are passed on.
    angular_rp_forest: bool (optional, default False)
        Whether to use an angular random projection forest to initialise
        the approximate nearest neighbor search. This can be faster, but is
        mostly on useful for metric that use an angular style distance such
        as cosine, correlation etc. In the case of those metrics angular forests
        will be chosen automatically.
    target_n_neighbors: int (optional, default -1)
        The number of nearest neighbors to use to construct the target simplcial
        set. If set to -1 use the ``n_neighbors`` value.
    target_metric: string or callable (optional, default 'categorical')
        The metric used to measure distance for a target array is using supervised
        dimension reduction. By default this is 'categorical' which will measure
        distance in terms of whether categories match or are different. Furthermore,
        if semi-supervised is required target values of -1 will be trated as
        unlabelled under the 'categorical' metric. If the target array takes
        continuous values (e.g. for a regression problem) then metric of 'l1'
        or 'l2' is probably more appropriate.
    target_weight: float (optional, default 0.5)
        weighting factor between data topology and target topology. A value of
        0.0 weights entirely on data, a value of 1.0 weights entirely on target.
        The default of 0.5 balances the weighting equally between data and target.
    transform_seed: int (optional, default 42)
        Random seed used for the stochastic aspects of the transform operation.
        This ensures consistency in transform operations.
    verbose: bool (optional, default False)
        Controls verbosity of logging.
    """

    def __init__(
        self,
        n_neighbors=15,
        n_components=2,
        n_epochs=None,
        learning_rate=1.0,
        init="spectral",
        min_dist=0.1,
        spread=1.0,
            metric='cosine',
            metric_kwds=None,
        set_op_mix_ratio=1.0,
        local_connectivity=1.0,
        repulsion_strength=1.0,
        negative_sample_rate=5,
        transform_queue_size=4.0,
        a=None,
        b=None,
        random_state=None,
        angular_rp_forest=False,
        target_n_neighbors=-1,
        target_metric="categorical",
        target_metric_kwds=None,
        target_weight=0.5,
        transform_seed=42,
        verbose=False,
        ignore_abundance=False,
        downsample=False,
    ):

        self.n_neighbors = n_neighbors
        self.n_epochs = n_epochs
        self.init = init
        self.n_components = n_components
        self.repulsion_strength = repulsion_strength
        self.learning_rate = learning_rate

        self.metric = metric
        self.metric_kwds = metric_kwds

        self.spread = spread
        self.min_dist = min_dist
        self.set_op_mix_ratio = set_op_mix_ratio
        self.local_connectivity = local_connectivity
        self.negative_sample_rate = negative_sample_rate
        self.random_state = random_state
        self.angular_rp_forest = angular_rp_forest
        self.transform_queue_size = transform_queue_size
        self.target_n_neighbors = target_n_neighbors
        self.target_metric = target_metric
        self.target_metric_kwds = target_metric_kwds
        self.target_weight = target_weight
        self.transform_seed = transform_seed
        self.verbose = verbose

        self.a = a
        self.b = b

        # Signature comparison parameters
        self.ignore_abundance = ignore_abundance
        self.downsample = downsample

    def _validate_parameters(self):
        if self.set_op_mix_ratio < 0.0 or self.set_op_mix_ratio > 1.0:
            raise ValueError("set_op_mix_ratio must be between 0.0 and 1.0")
        if self.repulsion_strength < 0.0:
            raise ValueError("repulsion_strength cannot be negative")
        if self.min_dist > self.spread:
            raise ValueError("min_dist must be less than or equal to spread")
        if self.min_dist < 0.0:
            raise ValueError("min_dist must be greater than 0.0")
        if not isinstance(self.init, str) and not isinstance(self.init, np.ndarray):
            raise ValueError("init must be a string or ndarray")
        if isinstance(self.init, str) and self.init not in ("spectral", "random"):
            raise ValueError('string init values must be "spectral" or "random"')
        if (
            isinstance(self.init, np.ndarray)
            and self.init.shape[1] != self.n_components
        ):
            raise ValueError("init ndarray must match n_components value")
        if not isinstance(self.metric, str) and not callable(self.metric):
            raise ValueError("metric must be string or callable")
        if self.negative_sample_rate < 0:
            raise ValueError("negative sample rate must be positive")
        if self._initial_alpha < 0.0:
            raise ValueError("learning_rate must be positive")
        if self.n_neighbors < 2:
            raise ValueError("n_neighbors must be greater than 2")
        if self.target_n_neighbors < 2 and self.target_n_neighbors != -1:
            raise ValueError("target_n_neighbors must be greater than 2")
        if not isinstance(self.n_components, int):
            raise ValueError("n_components must be an int")
        if self.n_components < 1:
            raise ValueError("n_components must be greater than 0")
        if self.n_epochs is not None and (
            self.n_epochs <= 10 or not isinstance(self.n_epochs, int)
        ):
            raise ValueError("n_epochs must be a positive integer " "larger than 10")

    def fit(self, tree, y=None):
        """Fit X into an embedded space.
        Optionally use y for supervised dimension reduction.
        Parameters
        ----------
        X : array, shape (n_samples, n_features) or (n_samples, n_samples)
            If the metric is 'precomputed' X must be a square distance
            matrix. Otherwise it contains a sample per row. If the method
            is 'exact', X may be a sparse matrix of type 'csr', 'csc'
            or 'coo'.
        y : array, shape (n_samples)
            A target array for supervised dimension reduction. How this is
            handled is determined by parameters UMAP was instantiated with.
            The relevant attributes are ``target_metric`` and
            ``target_metric_kwds``.
        """

        self._raw_tree = tree

        # Handle all the optional arguments, setting default
        if self.a is None or self.b is None:
            self._a, self._b = find_ab_params(self.spread, self.min_dist)
        else:
            self._a = self.a
            self._b = self.b

        if self.metric_kwds is not None:
            self._metric_kwds = self.metric_kwds
        else:
            self._metric_kwds = {}

        if self.target_metric_kwds is not None:
            self._target_metric_kwds = self.target_metric_kwds
        else:
            self._target_metric_kwds = {}

        if isinstance(self.init, np.ndarray):
            init = check_array(self.init, dtype=np.float32, accept_sparse=False)
        else:
            init = self.init

        self._initial_alpha = self.learning_rate

        self._validate_parameters()

        if self.verbose:
            print(str(self))

        # Error check n_neighbors based on data size
        self._n_leaves = sum(1 for _ in tree.leaves())
        if self._n_leaves == 0:
            raise ValueError("Tree has zero leaves")

        if self._n_leaves <= self.n_neighbors:
            if self._n_leaves == 1:
                self.embedding_ = np.zeros(
                    (1, self.n_components)
                )  # needed to sklearn comparability
                return self

            warn(
                "n_neighbors is larger than the dataset size; truncating to "
                "n_leaves - 1"
            )
            self._n_neighbors = self._n_leaves - 1
        else:
            self._n_neighbors = self.n_neighbors

        # if scipy.sparse.isspmatrix_csr(X):
        #     if not X.has_sorted_indices:
        #         X.sort_indices()
        #     self._sparse_data = True
        # else:
        #     self._sparse_data = False

        random_state = check_random_state(self.random_state)

        if self.verbose:
            print("Construct fuzzy simplicial set")

        # # Handle small cases efficiently by computing all distances
        # if self._n_leaves < 500:
        #     self._small_data = True
        #     signatures = [node.data for node in tree.nodes.values() if isinstance(node, Leaf)]
        #
        #     dmat = compare_all_pairs(signatures,
        #                              ignore_abundance=self.ignore_abundance,
        #                              downsample=self.downsample)
        #     self.graph_ = fuzzy_simplicial_set(
        #         dmat,
        #         self._n_neighbors,
        #         random_state,
        #         "precomputed",
        #         self._metric_kwds,
        #         None,
        #         None,
        #         self.angular_rp_forest,
        #         self.set_op_mix_ratio,
        #         self.local_connectivity,
        #         self.verbose,
        #     )
        # else:
        self._small_data = False
        # Standard case
        adjacencies = tree.nearest_neighbors(self._n_neighbors,
                                             self.ignore_abundance,
                                             self.downsample)
        (self._knn_indices, self._knn_dists, self._rp_forest) = \
            tree.similarity_adjacency_to_knn(adjacencies)

        self.graph_ = fuzzy_simplicial_set_from_knn(
            self._knn_indices,
            self._knn_dists,
            self.n_neighbors,
            self.set_op_mix_ratio,
            self.local_connectivity,
        )
        self._search_graph = make_search_graph(self._knn_indices,
                                               self._knn_dists)

        if y is not None:
            if self._n_leaves != len(y):
                raise ValueError(
                    "Number of leaves = {n_leaves}, length of y = {len_y}, while it must be equal.".format(
                        n_leaves=self._n_leaves, len_y=len(y)
                    )
                )
            y_ = check_array(y, ensure_2d=False)
            if self.target_metric == "categorical":
                if self.target_weight < 1.0:
                    far_dist = 2.5 * (1.0 / (1.0 - self.target_weight))
                else:
                    far_dist = 1.0e12
                self.graph_ = categorical_simplicial_set_intersection(
                    self.graph_, y_, far_dist=far_dist
                )
            else:
                if self.target_n_neighbors == -1:
                    target_n_neighbors = self._n_neighbors
                else:
                    target_n_neighbors = self.target_n_neighbors

                # Handle the small case as precomputed as before
                if y.shape[0] < 4096:
                    ydmat = pairwise_distances(
                        y_[np.newaxis, :].T,
                        metric=self.target_metric,
                        **self._target_metric_kwds
                    )
                    target_graph = fuzzy_simplicial_set(
                        ydmat,
                        target_n_neighbors,
                        random_state,
                        "precomputed",
                        None,
                        None,
                        None,
                        False,
                        1.0,
                        1.0,
                        False,
                    )
                else:
                    # Standard case
                    target_graph = fuzzy_simplicial_set(
                        y_[np.newaxis, :].T,
                        target_n_neighbors,
                        random_state,
                        self.target_metric,
                        self._target_metric_kwds,
                        None,
                        None,
                        False,
                        1.0,
                        1.0,
                        False,
                    )
                # product = self.graph_.multiply(target_graph)
                # # self.graph_ = 0.99 * product + 0.01 * (self.graph_ +
                # #                                        target_graph -
                # #                                        product)
                # self.graph_ = product
                self.graph_ = general_simplicial_set_intersection(
                    self.graph_, target_graph, self.target_weight
                )
                self.graph_ = reset_local_connectivity(self.graph_)

        if self.n_epochs is None:
            n_epochs = 0
        else:
            n_epochs = self.n_epochs

        if self.verbose:
            print(ts(), "Construct embedding")

        self.embedding_ = simplicial_set_embedding(
            self._raw_tree,
            self.graph_,
            self.n_components,
            self._initial_alpha,
            self._a,
            self._b,
            self.repulsion_strength,
            self.negative_sample_rate,
            n_epochs,
            init,
            random_state,
            self.metric,
            self._metric_kwds,
            self.verbose,
        )

        if self.verbose:
            print(ts() + " Finished embedding")

        # self._input_hash = joblib.hash(self._raw_tree)

        return self

    def fit_transform(self, tree, y=None):
        """Fit X into an embedded space and return that transformed
        output.
        Parameters
        ----------
        X : array, shape (n_samples, n_features) or (n_samples, n_samples)
            If the metric is 'precomputed' X must be a square distance
            matrix. Otherwise it contains a sample per row.
        y : array, shape (n_samples)
            A target array for supervised dimension reduction. How this is
            handled is determined by parameters UMAP was instantiated with.
            The relevant attributes are ``target_metric`` and
            ``target_metric_kwds``.
        Returns
        -------
        X_new : array, shape (n_samples, n_components)
            Embedding of the training data in low-dimensional space.
        """
        self.fit(tree, y)
        return self.embedding_

    def transform(self, tree):
        """Transform X into the existing embedded space and return that
        transformed output.
        Parameters
        ----------
        X : array, shape (n_samples, n_features)
            New data to be transformed.
        Returns
        -------
        X_new : array, shape (n_samples, n_components)
            Embedding of the new data in low-dimensional space.
        """
        # If we fit just a single instance then error
        if self.embedding_.shape[0] == 1:
            raise ValueError(
                "Transform unavailable when model was fit with"
                "only a single data sample."
            )
        # If we just have the original input then short circuit things
        X = check_array(X, dtype=np.float32, accept_sparse="csr")
        x_hash = joblib.hash(X)
        if x_hash == self._input_hash:
            return self.embedding_

        if self._sparse_data:
            raise ValueError("Transform not available for sparse input.")
        elif self.metric == "precomputed":
            raise ValueError(
                "Transform  of new data not available for " "precomputed metric."
            )

        X = check_array(X, dtype=np.float32, order="C")
        random_state = check_random_state(self.transform_seed)
        rng_state = random_state.randint(INT32_MIN, INT32_MAX, 3).astype(np.int64)

        if self._small_data:
            dmat = pairwise_distances(
                X, self._raw_data, metric=self.metric, **self._metric_kwds
            )
            indices = np.argpartition(dmat, self._n_neighbors)[:, : self._n_neighbors]
            dmat_shortened = submatrix(dmat, indices, self._n_neighbors)
            indices_sorted = np.argsort(dmat_shortened)
            indices = submatrix(indices, indices_sorted, self._n_neighbors)
            dists = submatrix(dmat_shortened, indices_sorted, self._n_neighbors)
        else:
            init = initialise_search(
                self._rp_forest,
                self._raw_data,
                X,
                int(self._n_neighbors * self.transform_queue_size),
                self._random_init,
                self._tree_init,
                rng_state,
            )
            result = self._search(
                self._raw_data,
                self._search_graph.indptr,
                self._search_graph.indices,
                init,
                X,
            )

            indices, dists = deheap_sort(result)
            indices = indices[:, : self._n_neighbors]
            dists = dists[:, : self._n_neighbors]

        adjusted_local_connectivity = max(0, self.local_connectivity - 1.0)
        sigmas, rhos = smooth_knn_dist(
            dists, self._n_neighbors, local_connectivity=adjusted_local_connectivity
        )

        rows, cols, vals = compute_membership_strengths(indices, dists, sigmas, rhos)

        graph = scipy.sparse.coo_matrix(
            (vals, (rows, cols)), shape=(X.shape[0], self._raw_data.shape[0])
        )

        # This was a very specially constructed graph with constant degree.
        # That lets us do fancy unpacking by reshaping the csr matrix indices
        # and data. Doing so relies on the constant degree assumption!
        csr_graph = normalize(graph.tocsr(), norm="l1")
        inds = csr_graph.indices.reshape(X.shape[0], self._n_neighbors)
        weights = csr_graph.data.reshape(X.shape[0], self._n_neighbors)
        embedding = init_transform(inds, weights, self.embedding_)

        if self.n_epochs is None:
            # For smaller datasets we can use more epochs
            if graph.shape[0] <= 10000:
                n_epochs = 100
            else:
                n_epochs = 30
        else:
            n_epochs = self.n_epochs // 3.0

        graph.data[graph.data < (graph.data.max() / float(n_epochs))] = 0.0
        graph.eliminate_zeros()

        epochs_per_sample = make_epochs_per_sample(graph.data, n_epochs)

        head = graph.row
        tail = graph.col

        embedding = optimize_layout(
            embedding,
            self.embedding_.astype(np.float32, copy=False),  # Fix #179
            head,
            tail,
            n_epochs,
            graph.shape[1],
            epochs_per_sample,
            self._a,
            self._b,
            rng_state,
            self.repulsion_strength,
            self._initial_alpha,
            self.negative_sample_rate,
            verbose=self.verbose,
        )

        return embedding