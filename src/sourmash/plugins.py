"""
Support for plugins to sourmash via importlib.metadata entrypoints.

Plugin entry point names:
* 'sourmash.load_from' - Index class loading.
* 'sourmash.save_to' - Signature saving.
* 'sourmash.picklist_filters' - extended Picklist functionality.

CTB TODO:

* consider using something other than 'name' for loader fn name. Maybe __doc__?
* try implement picklist plugin?
"""

DEFAULT_LOAD_FROM_PRIORITY = 99
DEFAULT_SAVE_TO_PRIORITY = 99

from .logging import debug_literal

# cover for older versions of Python that don't support selection on load
# (the 'group=' below).
from importlib.metadata import entry_points

# load 'load_from' entry points. NOTE: this executes on import of this module.
try:
    _plugin_load_from = entry_points(group='sourmash.load_from')
except TypeError:
    from importlib_metadata import entry_points
    _plugin_load_from = entry_points(group='sourmash.load_from')

# load 'save_to' entry points as well.
_plugin_save_to = entry_points(group='sourmash.save_to')


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


def get_save_to_functions():
    "Load the 'save_to' plugins and yield tuples (priority, fn)."
    debug_literal(f"save_to plugins: {_plugin_save_to}")

    # Load each plugin,
    for plugin in _plugin_save_to:
        save_cls = plugin.load()

        # get 'priority' if it is available
        priority = getattr(save_cls, 'priority', DEFAULT_SAVE_TO_PRIORITY)

        # retrieve name (which is specified by plugin?)
        name = plugin.name
        debug_literal(f"plugins.save_to_functions: got '{name}', priority={priority}")
        yield priority, save_cls


def get_cli_scripts_descriptions():
    "Get the descriptions for command-line plugins."
    plugin_list = entry_points(group='sourmash.cli_script')
    for plugin in plugin_list:
        name = plugin.name
        script_cls = plugin.load()
        yield script_cls.helpstring


def add_cli_scripts(parser):
    "Configure parsing for command-line plugins."
    plugin_list = entry_points(group='sourmash.cli_script')
    for plugin in plugin_list:
        name = plugin.name
        script_cls = plugin.load()
        subparser = parser.add_parser(script_cls.command)
        debug_literal(f"cls_script plugin '{name}' adding command '{script_cls.command}'")
        obj = script_cls(subparser)

    return 0
