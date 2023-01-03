"""
Test the plugin framework in sourmash.plugins, which uses importlib.metadata
entrypoints.
"""

import pytest
import sourmash

import sourmash_tst_utils as utils
from sourmash import plugins
from sourmash.index import LinearIndex


class FakeEntryPoint:
    """
    A class that stores a name and an object to be returned on 'load()'.
    Mocks the EntryPoint class used by importlib.metadata.
    """
    def __init__(self, name, load_obj):
        self.name = name
        self.load_obj = load_obj

    def load(self):
        return self.load_obj

#
# Test basic features of the load_from plugin.
#

class Test_EntryPointBasics_LoadFrom:
    def get_some_sigs(self, location, *args, **kwargs):
        ss2 = utils.get_test_data('2.fa.sig')
        ss47 = utils.get_test_data('47.fa.sig')
        ss63 = utils.get_test_data('63.fa.sig')

        sig2 = sourmash.load_one_signature(ss2, ksize=31)
        sig47 = sourmash.load_one_signature(ss47, ksize=31)
        sig63 = sourmash.load_one_signature(ss63, ksize=31)

        lidx = LinearIndex([sig2, sig47, sig63], location)

        return lidx
    get_some_sigs.priority = 1
        
    def setup_method(self):
        self.saved_plugins = plugins._plugin_load_from
        plugins._plugin_load_from = [FakeEntryPoint('test_load', self.get_some_sigs)]

    def teardown_method(self):
        plugins._plugin_load_from = self.saved_plugins

    def test_load_1(self, runtmp):
        ps = list(plugins.get_load_from_functions())
        assert len(ps) == 1

        fake_location = runtmp.output('passed-through location')
        idx = sourmash.load_file_as_index(fake_location)
        print(idx, idx.location)

        assert len(idx) == 3
        assert idx.location == fake_location
