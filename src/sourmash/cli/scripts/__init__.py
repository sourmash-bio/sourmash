"""Provide a mechanism to add CLI plugins to sourmash.

See https://sourmash.readthedocs.io/en/latest/dev_plugins.html for docs,
src/sourmash/plugins.py for core sourmash implementation code, and
https://github.com/sourmash-bio/sourmash_plugin_template for a template repo
for making new plugins.
"""

# CTB TODO:
# * provide suggestions for documentation & metadata for authors:
# * provide guidance on how to test your CLI plugin at the CLI
#   (minimal testing regime: sourmash scripts, look for description etc.)

import argparse
import sourmash

# Here, we decorate this module with the various extension objects
# e.g. 'sourmash scripts foo' will look up attribute 'scripts.foo'
# and we will return the extension class object, which will then
# be run by sourmash.__main__. This dictionary is loaded below
# by sourmash.plugins.add_cli_scripts.
_extension_dict = {}

def __getattr__(name):
    if name in _extension_dict:
        return _extension_dict[name]
    raise AttributeError(name)

def subparser(subparsers):
    subparser = subparsers.add_parser('scripts',
                                      usage=argparse.SUPPRESS,
                                      formatter_class=argparse.RawDescriptionHelpFormatter,
                                      aliases=['ext'])

    # get individual help strings:
    descrs = list(sourmash.plugins.get_cli_scripts_descriptions())
    if descrs:
        description = "\n".join(descrs)
    else:
        description = "(No script plugins detected!)"

    s = subparser.add_subparsers(title="available plugin/extension commands",
                                 dest='subcmd',
                                 metavar='subcmd',
                                 help=argparse.SUPPRESS,
                                 description=description)

    _extension_dict.update(sourmash.plugins.add_cli_scripts(s))
