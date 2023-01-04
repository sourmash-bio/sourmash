## sourmash plugins via Python entry points

As of version 4.7.0, sourmash has experimental support for Python
plugins to load and save signatures in different ways (e.g. file
formats, RPC servers, databases, etc.).  This support is provided via
the "entry points" mechanism supplied by
[`importlib.metadata`](https://docs.python.org/3/library/importlib.metadata.html)
and documented
[here](https://setuptools.pypa.io/en/latest/userguide/entry_point.html).

```{note}
Note: The plugin API is _not_ finalized or subject to semantic
versioning just yet!  Please subscribe to
[sourmash#1353](https://github.com/sourmash-bio/sourmash/issues/1353)
if you want to keep up to date on plugin support.
```

You can define entry points in the `pyproject.toml` file
like so:

```
[project.entry-points."sourmash.load_from"]
a_reader = "module_name:load_sketches"

[project.entry-points."sourmash.save_to"]
a_writer = "module_name:SaveSignatures_WriteFile"
```

Here, `module_name` should be the name of the module to import.
`load_sketches` should be a function that takes a location along with
arbitrary keyword arguments. `SaveSignatures_WriteFile` should be a
class that subclasses `BaseSave_SignaturesToLocation`. See the
`sourmash.save_load` module for the saving and loading code already
included in sourmash.

Note that if the function or class has a `priority` attribute, this will
be used to determine the order in which the plugins are called.

## Examples

Some beta plugins are available as examples:

* [sourmash-bio/sourmash_plugin_load_urls](https://github.com/sourmash-bio/sourmash_plugin_load_urls) - load signatures and CSV manifests via [fsspec](https://filesystem-spec.readthedocs.io/).
* [sourmash-bio/sourmash_plugin_avro](https://github.com/sourmash-bio/sourmash_plugin_avro) - use [Apache Avro](https://avro.apache.org/) as a serialization format.

