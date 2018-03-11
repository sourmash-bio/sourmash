=======================
``sourmash`` Python API
=======================

The primary programmatic way of interacting with ``sourmash`` is via
its Python API.  (The core MinHash Python API closely mirrors the
underlying C++ code, but for now this is undocumented.)

.. contents::
   :depth: 2

``MinHash``: basic MinHash sketch functionality
===============================================

.. automodule:: sourmash
   :members:

``SourmashSignature``: save and load MinHash sketches in JSON
=============================================================

.. automodule:: sourmash.signature
   :members:

``SBT``: save and load Sequence Bloom Trees in JSON
=============================================================

.. automodule:: sourmash.sbt
   :members: GraphFactory, Node, NodePos, SBT, Leaf
   :undoc-members:

``sourmash.fig``: make plots and figures
============================================

.. automodule:: sourmash.fig
   :members:
