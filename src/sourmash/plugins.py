"""
Support for plugins to sourmash via importlib.metadata entrypoints.

Plugin entry point names:
* 'sourmash.load_from' - Index class loading.
* 'sourmash.save_to' - Signature saving.
* 'sourmash.cli_script' - command-line extension.

CTB TODO:

* consider using something other than 'name' for loader fn name. Maybe __doc__?
* try implement picklist plugin?
"""

DEFAULT_LOAD_FROM_PRIORITY = 99
DEFAULT_SAVE_TO_PRIORITY = 99

import itertools

from .logging import (debug_literal, error, notify, set_quiet)

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
_plugin_cli = entry_points(group='sourmash.cli_script')
_plugin_cli_once = False

###

def get_load_from_functions():
    "Load the 'load_from' plugins and yield tuples (priority, name, fn)."
    debug_literal(f"load_from plugins: {_plugin_load_from}")

    # Load each plugin,
    for plugin in _plugin_load_from:
        try:
            loader_fn = plugin.load()
        except (ModuleNotFoundError, AttributeError) as e:
            debug_literal(f"plugins.load_from_functions: got error loading {plugin.name}: {str(e)}")
            continue

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
        try:
            save_cls = plugin.load()
        except (ModuleNotFoundError, AttributeError) as e:
            debug_literal(f"plugins.load_from_functions: got error loading {plugin.name}: {str(e)}")
            continue

        # get 'priority' if it is available
        priority = getattr(save_cls, 'priority', DEFAULT_SAVE_TO_PRIORITY)

        # retrieve name (which is specified by plugin?)
        name = plugin.name
        debug_literal(f"plugins.save_to_functions: got '{name}', priority={priority}")
        yield priority, save_cls


class CommandLinePlugin:
    """
    Provide some minimal common CLI functionality - -q and -d.

    Subclasses should call super().__init__(parser) and super().main(args).
    """
    command = None
    description = None

    def __init__(self, parser):
        parser.add_argument(
            '-q', '--quiet', action='store_true',
            help='suppress non-error output'
        )
        parser.add_argument(
            '-d', '--debug', action='store_true',
            help='provide debugging output'
        )

    def main(self, args):
        set_quiet(args.quiet, args.debug)


def get_cli_script_plugins():
    global _plugin_cli_once

    x = []
    for plugin in _plugin_cli:
        name = plugin.name
        mod = plugin.module
        try:
            script_cls = plugin.load()
        except (ModuleNotFoundError, AttributeError):
            if _plugin_cli_once is False:
                error(f"ERROR: cannot find or load module for cli_script plugin '{name}'")
            continue

        command = getattr(script_cls, 'command', None)
        if command is None:
            # print error message only once...
            if _plugin_cli_once is False:
                error(f"ERROR: no command provided by cli_script plugin '{name}' from {mod}; skipping")
        else:
            x.append(plugin)

    _plugin_cli_once = True
    return x


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
                              _plugin_cli)
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
