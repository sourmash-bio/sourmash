"""
Test the plugin framework in sourmash.plugins, which uses importlib.metadata
entrypoints.
"""

import pytest
import sourmash

import sourmash_tst_utils as utils
from sourmash import plugins


class FakeEntryPoint:
    def __init__(self, name, load_obj):
        self.name = name
        self.load_obj = load_obj

    def load(self):
        return self.load_obj

class Test_EntryPointBasics_LoadFrom:
    def method_was_called(self):
        self.call_flag = True
        
    def setup_method(self):
        self.saved_plugins = plugins._plugin_load_from
        plugins._plugin_load_from = [FakeEntryPoint('test_load', self.method_was_called)]

    def teardown_method(self):
        plugins._plugin_load_from = self.saved_plugins
    
    def test_load_1(runtmp):
        ps = list(plugins.get_load_from_functions())
        print()
        assert len(ps) == 1

    def test_load_2(runtmp):
        print('hello, world, 2')
