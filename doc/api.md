# `sourmash` Python API

The primary programmatic way of interacting with `sourmash` is via
its Python API.

**Please also see [examples of using the API](api-example.md).**

```{contents}
   :depth: 2
```

## `MinHash`: basic MinHash sketch functionality

```{eval-rst}
.. autoclass:: sourmash.MinHash
   :members:

   .. automethod:: __init__
```

## `SourmashSignature`: save and load MinHash sketches in JSON

```{eval-rst}
.. automodule:: sourmash.signature
   :members:
```

## `SBT`: save and load Sequence Bloom Trees in JSON

```{eval-rst}
.. automodule:: sourmash.sbt
   :members: GraphFactory, Node, NodePos, SBT, Leaf
   :undoc-members:
```

## `sourmash.fig`: make plots and figures

```{eval-rst}
.. automodule:: sourmash.fig
   :members:
```
