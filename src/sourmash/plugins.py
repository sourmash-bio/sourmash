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

import itertools

from .logging import debug_literal, error, notify

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

# aaaaand CLI entry points:
_plugin_cli_cache = None

###

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


def get_cli_script_plugins():
    global _plugin_cli_cache

    if _plugin_cli_cache is None:
        _plugin_cli_cache = []
        for plugin in entry_points(group='sourmash.cli_script'):
            name = plugin.name
            mod = plugin.module
            script_cls = plugin.load()

            command = getattr(script_cls, 'command', None)
            if command is None:
                error(f"ERROR: no command provided by cli_script plugin '{name}' from {mod}; skipping")

            _plugin_cli_cache.append(plugin)

    return _plugin_cli_cache


def get_cli_scripts_descriptions():
    "Build the descriptions for command-line plugins."
    for plugin in get_cli_script_plugins():
        name = plugin.name
        script_cls = plugin.load()
        command = getattr(script_cls, 'command')
        description = getattr(script_cls, 'description',
                              f"(no description provided by plugin '{name}')")
        yield f"sourmash scripts {command:16s} - {description}"


def add_cli_scripts(parser):
    "Configure parsing for command-line plugins."
    d = {}

    # @CTB: factor out common code
    for plugin in get_cli_script_plugins():
        name = plugin.name
        script_cls = plugin.load()

        subparser = parser.add_parser(script_cls.command)
        debug_literal(f"cls_script plugin '{name}' adding command '{script_cls.command}'")
        obj = script_cls(subparser)
        d[script_cls.command] = obj

    return d


def list_all_plugins():
    plugins = itertools.chain(_plugin_load_from,
                              _plugin_save_to,
                              _plugin_cli_cache)
    plugins = list(plugins)

    if not plugins:
        notify("\n(no plugins detected)\n")

    notify("")
    notify("the following plugins are installed:")
    notify("")
    notify(f"{'plugin type':<20s} {'from python module':<30s} {'v':<5s} {'entry point name':<20s}")
    notify(f"{'-'*20} {'-'*30} {'-'*5} {'-'*20}")

    for plugin in plugins:
        name = plugin.name
        mod = plugin.module
        version = plugin.dist.version
        group = plugin.group

        notify(f"{group:<20s} {mod:<30s} {version:<5s} {name:<20s}")
