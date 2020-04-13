"""
Compute nearest neighbor graphs on SBTs for clustering and 2d embedding
"""
import math

import numpy as np

from .sbt import Leaf
from .sbtmh import SigLeaf

INT32_MIN = np.iinfo(np.int32).min + 1
INT32_MAX = np.iinfo(np.int32).max - 1


# def get_leaves(tree):
#     for i, node in tree.nodes.items():
#         if isinstance(node, (Leaf, SigLeaf)):
#             yield i, node
#
#
# def nearest_neighbors(tree, n_neighbors, ignore_abundance, downsample):
#     """Yield leaf1, leaf2, similarity"""
#     n_parent_levels = math.log2(n_neighbors) + 1
#
#     adjacencies = []
#
#     for position1, leaf1 in get_leaves(tree):
#         n = 1
#         upper_internal_node = tree.parent(position1)
#         while n < n_parent_levels:
#             upper_internal_node = tree.parent(upper_internal_node.pos)
#             n += 1
#         #         print("upper_internal_node:", upper_internal_node)
#         leaves = tree.leaves_under(upper_internal_node.pos)
#
#         similarities = []
#         for leaf2 in leaves:
#             if leaf2 == leaf1:
#                 continue
#             similarity = leaf1.data.similarity(
#                 leaf2.data, ignore_abundance=ignore_abundance,
#                 downsample=downsample)
#             # Don't filter for minimum similarity as some samples are bad
#             # and don't have many neighbors, and arrays containing lists
#             # of multiple sizes won't cast to int/float in numpy, so
#             # similarity_adjacency_to_knn won't work on the adjacencies
#             similarities.append(
#                 [leaf1.data.name(), leaf2.data.name(), similarity])
#         adjacent = sorted(similarities, key=lambda x: x[1])[-n_neighbors:]
#         adjacencies.extend(adjacent)
#     return adjacencies


