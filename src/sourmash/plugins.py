"""
Support for plugins to sourmash via importlib.metadata entrypoints.

Plugin entry point names:
* 'sourmash.load_from' - Index class loading.
* 'sourmash.save_to' - Signature saving.
* 'sourmash.picklist_filters' - extended Picklist functionality.

CTB TODO:

* is there a way to provide attributes like 'priority' in the pyproject.toml
  of the plugin?
* consider using something other than 'name' for loader fn name. Maybe __doc__?
"""

DEFAULT_LOAD_FROM_PRIORITY = 99

from .logging import error, debug_literal

# cover for older versions of Python that don't support selection on load
# (the 'group=' below).
from importlib.metadata import entry_points

# load 'load_from' entry points.
try:
    _plugin_load_from = entry_points(group='sourmash.load_from')
except TypeError:
    from importlib_metadata import entry_points
    _plugin_load_from = entry_points(group='sourmash.load_from')


def get_load_from_functions():
    "Load the 'load_from' plugins and yield tuples (priority, name, fn)."
    debug_literal(f"load_from plugins: {_plugin_load_from}")

    # Load each plugin,
    for plugin in _plugin_load_from:
        loader_fn = plugin.load()

        # get 'priority' if it is available
        priority = getattr(loader_fn, 'priority', DEFAULT_LOAD_FROM_PRIORITY)

        # retrieve name (which is specified by plugin?)
        name = plugin.name
        debug_literal(f"plugins.load_from_functions: got '{name}', priority={priority}")
        yield priority, name, loader_fn