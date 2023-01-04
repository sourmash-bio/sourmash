"""
Test the plugin framework in sourmash.plugins, which uses importlib.metadata
entrypoints.

CTB TODO:
* check name?
"""

import pytest
import sourmash
from sourmash.logging import set_quiet

import sourmash_tst_utils as utils
from sourmash import plugins
from sourmash.index import LinearIndex
from sourmash.save_load import (Base_SaveSignaturesToLocation,
                                SaveSignaturesToLocation)


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
# Test basic features of the load_from plugin hook.
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

    def test_load_1(self):
        ps = list(plugins.get_load_from_functions())
        assert len(ps) == 1

    def test_load_2(self, runtmp):
        fake_location = runtmp.output('passed-through location')
        idx = sourmash.load_file_as_index(fake_location)
        print(idx, idx.location)

        assert len(idx) == 3
        assert idx.location == fake_location


class Test_EntryPoint_LoadFrom_Priority:
    def get_some_sigs(self, location, *args, **kwargs):
        ss2 = utils.get_test_data('2.fa.sig')
        ss47 = utils.get_test_data('47.fa.sig')
        ss63 = utils.get_test_data('63.fa.sig')

        sig2 = sourmash.load_one_signature(ss2, ksize=31)
        sig47 = sourmash.load_one_signature(ss47, ksize=31)
        sig63 = sourmash.load_one_signature(ss63, ksize=31)

        lidx = LinearIndex([sig2, sig47, sig63], location)

        return lidx
    get_some_sigs.priority = 5

    def set_called_flag_1(self, location, *args, **kwargs):
        # high priority 1, raise ValueError
        print('setting flag 1')
        self.was_called_flag_1 = True
        raise ValueError
    set_called_flag_1.priority = 1

    def set_called_flag_2(self, location, *args, **kwargs):
        # high priority 2, return None
        print('setting flag 2')
        self.was_called_flag_2 = True

        return None
    set_called_flag_2.priority = 2

    def set_called_flag_3(self, location, *args, **kwargs):
        # lower priority 10, should not be called
        print('setting flag 3')
        self.was_called_flag_3 = True

        return None
    set_called_flag_3.priority = 10

    def setup_method(self):
        self.saved_plugins = plugins._plugin_load_from
        plugins._plugin_load_from = [
            FakeEntryPoint('test_load', self.get_some_sigs),
            FakeEntryPoint('test_load_2', self.set_called_flag_1),
            FakeEntryPoint('test_load_3', self.set_called_flag_2),
            FakeEntryPoint('test_load_4', self.set_called_flag_3)
            ]
        self.was_called_flag_1 = False
        self.was_called_flag_2 = False
        self.was_called_flag_3 = False

    def teardown_method(self):
        plugins._plugin_load_from = self.saved_plugins

    def test_load_1(self):
        ps = list(plugins.get_load_from_functions())
        assert len(ps) == 4

        assert not self.was_called_flag_1
        assert not self.was_called_flag_2
        assert not self.was_called_flag_3

    def test_load_2(self, runtmp):
        fake_location = runtmp.output('passed-through location')
        idx = sourmash.load_file_as_index(fake_location)
        print(idx, idx.location)

        assert len(idx) == 3
        assert idx.location == fake_location

        assert self.was_called_flag_1
        assert self.was_called_flag_2
        assert not self.was_called_flag_3


#
# Test basic features of the save_to plugin hook.
#

class FakeSaveClass(Base_SaveSignaturesToLocation):
    """
    A fake save class that just records what was sent to it.
    """
    priority = 50

    def __init__(self, location):
        super().__init__(location)
        self.keep = []

    @classmethod
    def matches(cls, location):
        if location:
            return location.endswith('.this-is-a-test')

    def add(self, ss):
        super().add(ss)
        self.keep.append(ss)


class FakeSaveClass_HighPriority(FakeSaveClass):
    priority = 1


class Test_EntryPointBasics_SaveTo:
    # test the basics
    def setup_method(self):
        self.saved_plugins = plugins._plugin_save_to
        plugins._plugin_save_to = [FakeEntryPoint('test_save', FakeSaveClass)]

    def teardown_method(self):
        plugins._plugin_save_to = self.saved_plugins

    def test_save_1(self):
        ps = list(plugins.get_save_to_functions())
        print(ps)
        assert len(ps) == 1

    def test_save_2(self, runtmp):
        # load some signatures to save
        ss2 = utils.get_test_data('2.fa.sig')
        ss47 = utils.get_test_data('47.fa.sig')
        ss63 = utils.get_test_data('63.fa.sig')

        sig2 = sourmash.load_one_signature(ss2, ksize=31)
        sig47 = sourmash.load_one_signature(ss47, ksize=31)
        sig63 = sourmash.load_one_signature(ss63, ksize=31)

        # build a fake location that matches the FakeSaveClass
        # extension
        fake_location = runtmp.output('out.this-is-a-test')

        # this should use the plugin architecture to return an object
        # of type FakeSaveClass, with the three signatures in it.
        x = SaveSignaturesToLocation(fake_location)
        with x as save_sig:
            save_sig.add(sig2)
            save_sig.add(sig47)
            save_sig.add(sig63)

        print(len(x))
        print(x.keep)

        assert isinstance(x, FakeSaveClass)
        assert x.keep == [sig2, sig47, sig63]


class Test_EntryPointPriority_SaveTo:
    # test that priority is observed

    def setup_method(self):
        self.saved_plugins = plugins._plugin_save_to
        plugins._plugin_save_to = [
            FakeEntryPoint('test_save', FakeSaveClass),
            FakeEntryPoint('test_save2', FakeSaveClass_HighPriority),
        ]

    def teardown_method(self):
        plugins._plugin_save_to = self.saved_plugins

    def test_save_1(self):
        ps = list(plugins.get_save_to_functions())
        print(ps)
        assert len(ps) == 2

    def test_save_2(self, runtmp):
        # load some signatures to save
        ss2 = utils.get_test_data('2.fa.sig')
        ss47 = utils.get_test_data('47.fa.sig')
        ss63 = utils.get_test_data('63.fa.sig')

        sig2 = sourmash.load_one_signature(ss2, ksize=31)
        sig47 = sourmash.load_one_signature(ss47, ksize=31)
        sig63 = sourmash.load_one_signature(ss63, ksize=31)

        # build a fake location that matches the FakeSaveClass
        # extension
        fake_location = runtmp.output('out.this-is-a-test')

        # this should use the plugin architecture to return an object
        # of type FakeSaveClass, with the three signatures in it.
        x = SaveSignaturesToLocation(fake_location)
        with x as save_sig:
            save_sig.add(sig2)
            save_sig.add(sig47)
            save_sig.add(sig63)

        print(len(x))
        print(x.keep)

        assert isinstance(x, FakeSaveClass_HighPriority)
        assert x.keep == [sig2, sig47, sig63]
        assert x.priority == 1
