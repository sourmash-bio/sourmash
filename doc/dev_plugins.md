# sourmash plugins via Python entry points

As of version 4.7.0, sourmash has experimental support for Python
plugins to load and save signatures in different ways (e.g. file
formats, RPC servers, databases, etc.) and to run additional commands
via the command-line.  This support is provided via
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

[project.entry-points."sourmash.cli_script"]
new_cli = "module_name:Command_NewCommand"
```

Here, `module_name` should be the name of the module to import.

* `load_sketches` should be a function that takes a location along with
arbitrary keyword arguments and returns an `Index` object
(e.g. `LinearIndex` for a collection of in-memory
signatures).
* `SaveSignatures_WriteFile` should be a class that
subclasses `BaseSave_SignaturesToLocation` and implements its own
mechanisms of saving signatures. See the `sourmash.save_load` module
for saving and loading code already used in sourmash.
* `Command_NewCommand` should be a class that subclasses
  `plugins.CommandLinePlugin` and provides an `__init__` and
  `main` method.

Note that if the reader function or writer class has a `priority`
attribute, this will be used to determine the order in which the
plugins are called.

The `name` attribute of the plugin (`a_reader`, `a_writer`, and `new_cli` in
`pyproject.toml`, above) is only used in debugging.

You can provide zero or more plugins, and you can define just a reader, or
just a writer, or just a CLI plugin.

## Templates and examples

If you want to create your own plug-in, you can start with the
[sourmash_plugin_template](https://github.com/sourmash-bio/sourmash_plugin_template) repo.

Some (early stage) plugins are also available as examples:

* [sourmash-bio/sourmash_plugin_load_urls](https://github.com/sourmash-bio/sourmash_plugin_load_urls) - load signatures and CSV manifests via [fsspec](https://filesystem-spec.readthedocs.io/).
* [sourmash-bio/sourmash_plugin_avro](https://github.com/sourmash-bio/sourmash_plugin_avro) - use [Apache Avro](https://avro.apache.org/) as a serialization format.

## Debugging plugins

`sourmash info -v` will list all detected plugins of any known kind.

`sourmash sig cat <input sig> -o <output sig>` is a simple way to
invoke a `save_to` plugin. Use `-d` to turn on debugging output.

`sourmash sig describe <input location>` is a simple way to invoke
a `load_from` plugin. Use `-d` to turn on debugging output.

`sourmash scripts` will list available command-line plugins.

## Semantic versioning and listing sourmash as a dependency

Plugins should generally list sourmash as a dependency for installation.

Once plugins are officially supported by sourmash, the plugin API will
be under [semantic versioning constraints](https://semver.org/). That
means that you should constrain plugins to depend on sourmash only up
to the next major version, e.g. sourmash v5.

Specifically, we suggest placing something like:
```
dependencies = ['sourmash>=4.8.0,<5']
```
in your `pyproject.toml` file.
