"""Provide a mechanism to add CLI plugins to sourmash.

See https://sourmash.readthedocs.io/en/latest/dev_plugins.html for docs,
src/sourmash/plugins.py for core sourmash implementation code, and
https://github.com/sourmash-bio/sourmash_plugin_template for a template repo
for making new plugins.
"""

# CTB TODO:
# * provide suggestions for URLs etc.
# * provide guidance on how to test your CLI plugin at the CLI
#   (sourmash scripts, look for description etc.)
# * is there any reason to provide a callback mechanism of any kind?
# * provide default command class that implements -q and -d

import argparse
import sourmash

# decorate this module with the various extension objects
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
                                      formatter_class=argparse.RawDescriptionHelpFormatter)

    # get individual help strings:
    descrs = list(sourmash.plugins.get_cli_scripts_descriptions())
    s = subparser.add_subparsers(title="available plugin/extension commands",
                                 dest='subcmd',
                                 metavar='subcmd',
                                 help=argparse.SUPPRESS,
                                 description="\n".join(descrs))

    _extension_dict.update(sourmash.plugins.add_cli_scripts(s))
