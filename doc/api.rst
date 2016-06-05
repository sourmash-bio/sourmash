=======================
``sourmash`` Python API
=======================

The primary programmatic way of interacting with ``sourmash`` is via
its Python API.  (The core MinHash Python API closely mirrors the
underlying C++ code, but for now this is undocumented.)

.. contents::
   :depth: 2

``Estimators``: basic MinHash sketch functionality
==================================================

.. automodule:: sourmash_lib
   :members:

``SourmashSignature``: save and load MinHash sketches in YAML
=============================================================

.. automodule:: sourmash_lib.signature
   :members:

``sourmash_lib.fig``: make plots and figures
============================================

.. automodule:: sourmash_lib.fig
   :members:

