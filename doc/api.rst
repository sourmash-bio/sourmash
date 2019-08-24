=======================
``sourmash`` Python API
=======================

The primary programmatic way of interacting with ``sourmash`` is via
its Python API.  Please also see `examples of using the API <api-example.html>`_.

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
