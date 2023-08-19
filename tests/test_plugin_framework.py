"""
Test the plugin framework in sourmash.plugins, which uses importlib.metadata
entrypoints.
"""

import sys
import pytest
import collections

import sourmash
from sourmash.logging import set_quiet

import sourmash_tst_utils as utils
from sourmash import plugins
from sourmash.index import LinearIndex
from sourmash.save_load import (Base_SaveSignaturesToLocation,
                                SaveSignaturesToLocation)


_Dist = collections.namedtuple('_Dist', ['version'])
class FakeEntryPoint:
    """
    A class that stores a name and an object to be returned on 'load()'.
    Mocks the EntryPoint class used by importlib.metadata.
    """
    module = 'test_plugin_framework'
    dist = _Dist('0.1')
    group = 'groupfoo'

    def __init__(self, name, load_obj, *,
                 error_on_import=None):
        self.name = name
        self.load_obj = load_obj
        self.error_on_import = error_on_import

    def load(self):
        if self.error_on_import is not None:
            raise self.error_on_import("as requested")
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
        plugins._plugin_load_from = [FakeEntryPoint('test_load', self.get_some_sigs),
                                     FakeEntryPoint('test_load', self.get_some_sigs, error_on_import=ModuleNotFoundError)]

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
        plugins._plugin_save_to = [FakeEntryPoint('test_save', FakeSaveClass),
                                   FakeEntryPoint('test_save', FakeSaveClass, error_on_import=ModuleNotFoundError)]

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


#
# Test basic features of the save_to plugin hook.
#

class FakeCommandClass(plugins.CommandLinePlugin):
    """
    A fake CLI class.
    """
    command = 'nifty'
    description = "do somethin' nifty"

    def __init__(self, parser):
        super().__init__(parser)
        parser.add_argument('arg1')
        parser.add_argument('--other', action='store_true')
        parser.add_argument('--do-fail', action='store_true')

    def main(self, args):
        super().main(args)
        print(f"hello, world! argument is: {args.arg1}")
        print(f"other is {args.other}")

        if args.do_fail:
            return 1
        return 0


class Test_EntryPointBasics_Command:
    # test the basics
    def setup_method(self):
        _ = plugins.get_cli_script_plugins()
        self.saved_plugins = plugins._plugin_cli
        plugins._plugin_cli_once = False
        plugins._plugin_cli = [FakeEntryPoint('test_command',
                                              FakeCommandClass)]

    def teardown_method(self):
        plugins._plugin_cli = self.saved_plugins

    def test_empty(self, runtmp):
        # empty out script plugins...
        plugins._plugin_cli = []

        with pytest.raises(utils.SourmashCommandFailed):
            runtmp.sourmash('scripts')
        out = runtmp.last_result.out
        err = runtmp.last_result.err
        print(out)
        print(err)
        assert '(No script plugins detected!)' in out

    def test_cmd_0(self, runtmp):
        # test default output with some plugins
        with pytest.raises(utils.SourmashCommandFailed):
            runtmp.sourmash('scripts')

        out = runtmp.last_result.out
        err = runtmp.last_result.err
        print(out)
        print(err)
        assert "do somethin' nifty" in out
        assert "sourmash scripts nifty" in out

    def test_cmd_1(self):
        # test descriptions
        ps = list(plugins.get_cli_scripts_descriptions())
        print(ps)
        assert len(ps) == 1

        descr0 = ps[0]
        assert "do somethin' nifty" in descr0
        assert "sourmash scripts nifty" in descr0

    def test_cmd_2(self):
        # test get_cli_script_plugins function
        ps = list(plugins.get_cli_script_plugins())
        print(ps)
        assert len(ps) == 1

    def test_cmd_3(self, runtmp):
        # test ability to run 'nifty' ;)
        with pytest.raises(utils.SourmashCommandFailed):
            runtmp.sourmash('scripts', 'nifty')

        out = runtmp.last_result.out
        err = runtmp.last_result.err
        print(out)
        print(err)

        assert 'nifty: error: the following arguments are required: arg1' in err
        assert 'usage:  nifty [-h] [-q] [-d] [--other] [--do-fail] arg1' in err

    def test_cmd_4(self, runtmp):
        # test basic argument parsing etc
        runtmp.sourmash('scripts', 'nifty', '--other', 'some arg')

        out = runtmp.last_result.out
        err = runtmp.last_result.err
        print(out)
        print(err)

        assert 'other is True' in out
        assert 'hello, world! argument is: some arg' in out

    def test_cmd_5(self, runtmp):
        # test exit code passthru
        with pytest.raises(utils.SourmashCommandFailed):
            runtmp.sourmash('scripts', 'nifty', '--do-fail', 'some arg')

        status = runtmp.last_result.status
        out = runtmp.last_result.out
        err = runtmp.last_result.err
        print(out)
        print(err)
        print(status)

        assert 'other is False' in out
        assert 'hello, world! argument is: some arg' in out


class FakeCommandClass_Second(plugins.CommandLinePlugin):
    """
    A fake CLI class.
    """
    command = 'more_nifty'
    description = "do somethin' else nifty"

    def __init__(self, parser):
        super().__init__(parser)
        parser.add_argument('arg1')
        parser.add_argument('--other', action='store_true')
        parser.add_argument('--do-fail', action='store_true')

    def main(self, args):
        super().main(args)
        print(f"hello, world! argument is: {args.arg1}")
        print(f"other is {args.other}")

        if args.do_fail:
            return 1
        return 0


class FakeCommandClass_Broken_1:
    """
    A fake CLI class.
    """
    # command = 'more_nifty' # no command

    def __init__(self, parser):
        assert 0

    def main(self, args):
        assert 0


class FakeCommandClass_Broken_2:
    """
    A fake CLI class.
    """
    command = 'broken'
    # no description

    def __init__(self, parser):
        pass

    def main(self, args):
        return 0


class Test_EntryPointBasics_TwoCommands:
    # test a second command
    def setup_method(self):
        _ = plugins.get_cli_script_plugins()
        self.saved_plugins = plugins._plugin_cli
        plugins._plugin_cli_once = False
        plugins._plugin_cli = [FakeEntryPoint('test_command',
                                              FakeCommandClass),
                               FakeEntryPoint('test_command2',
                                              FakeCommandClass_Second),
                               FakeEntryPoint('test_command3',
                                              FakeCommandClass_Broken_1),
                               FakeEntryPoint('test_command4',
                                              FakeCommandClass_Broken_2),
                               FakeEntryPoint('error-on-import',
                                              FakeCommandClass,
                                           error_on_import=ModuleNotFoundError)
                               ]

    def teardown_method(self):
        plugins._plugin_cli = self.saved_plugins

    def test_cmd_0(self, runtmp):
        # test default output for a few plugins
        with pytest.raises(utils.SourmashCommandFailed):
            runtmp.sourmash('scripts')

        out = runtmp.last_result.out
        err = runtmp.last_result.err
        print(out)
        print(err)
        assert "do somethin' nifty" in out
        assert "sourmash scripts nifty" in out

        assert "do somethin' else nifty" in out
        assert "sourmash scripts more_nifty" in out

    def test_cmd_1(self, runtmp):
        # test 'nifty'
        runtmp.sourmash('scripts', 'nifty', 'some arg')

        status = runtmp.last_result.status
        out = runtmp.last_result.out
        err = runtmp.last_result.err
        print(out)
        print(err)
        print(status)

        assert 'other is False' in out
        assert 'hello, world! argument is: some arg' in out

    def test_cmd_2(self, runtmp):
        # test 'more_nifty'
        runtmp.sourmash('scripts', 'more_nifty', 'some arg')

        status = runtmp.last_result.status
        out = runtmp.last_result.out
        err = runtmp.last_result.err
        print(out)
        print(err)
        print(status)

        assert 'other is False' in out
        assert 'hello, world! argument is: some arg' in out

    def test_sourmash_info(self, runtmp):
        # test 'sourmash info -v' => shows the plugins
        runtmp.sourmash('info', '-v')

        out = runtmp.last_result.out
        err = runtmp.last_result.err
        print(out)
        print(err)

        expected = """
groupfoo             test_plugin_framework          0.1   test_command
groupfoo             test_plugin_framework          0.1   test_command2
groupfoo             test_plugin_framework          0.1   test_command3
groupfoo             test_plugin_framework          0.1   test_command4
""".splitlines()
        for line in expected:
            assert line in err


def test_cli_scripts_getattr_fail():
    # test scripts.__getattr__ w/fail
    from sourmash.cli import scripts

    with pytest.raises(AttributeError):
        scripts.ThisAttrDoesNotExist


def test_cli_scripts_getattr_succ():
    # test scripts.__getattr__ w/success
    from sourmash.cli import scripts

    scripts.subparser
