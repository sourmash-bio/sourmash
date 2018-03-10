"Various utilities used by sourmash tests."

from __future__ import print_function
import sys
import os
import tempfile
import shutil
import subprocess

import pkg_resources
from pkg_resources import Requirement, resource_filename, ResolutionError
import traceback
from io import open  # pylint: disable=redefined-builtin
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO


SIG_FILES = [os.path.join('demo', f) for f in (
  "SRR2060939_1.sig", "SRR2060939_2.sig", "SRR2241509_1.sig",
  "SRR2255622_1.sig", "SRR453566_1.sig", "SRR453569_1.sig", "SRR453570_1.sig")
]


def scriptpath(scriptname='sourmash'):
    """Return the path to the scripts, in both dev and install situations."""
    # note - it doesn't matter what the scriptname is here, as long as
    # it's some script present in this version of sourmash.

    path = os.path.join(os.path.dirname(__file__), "../")
    if os.path.exists(os.path.join(path, scriptname)):
        return path

    path = os.path.join(os.path.dirname(__file__), "../../EGG-INFO/")
    if os.path.exists(os.path.join(path, scriptname)):
        return path

    for path in os.environ['PATH'].split(':'):
        if os.path.exists(os.path.join(path, scriptname)):
            return path


def _runscript(scriptname):
    """Find & run a script with exec (i.e. not via os.system or subprocess)."""
    namespace = {"__name__": "__main__"}
    namespace['sys'] = globals()['sys']

    try:
        pkg_resources.load_entry_point("sourmash", 'console_scripts', scriptname)()
        return 0
    except pkg_resources.ResolutionError:
        pass

    path = scriptpath()

    scriptfile = os.path.join(path, scriptname)
    if os.path.isfile(scriptfile):
        if os.path.isfile(scriptfile):
            exec(  # pylint: disable=exec-used
                compile(open(scriptfile).read(), scriptfile, 'exec'),
                namespace)
            return 0

    return -1


def runscript(scriptname, args, in_directory=None,
              fail_ok=False):
    """Run a Python script using exec().

    Run the given Python script, with the given args, in the given directory,
    using 'exec'.  Mimic proper shell functionality with argv, and capture
    stdout and stderr.

    When using :attr:`fail_ok`=False in tests, specify the expected error.
    """
    __tracebackhide__ = True
    sysargs = [scriptname]
    sysargs.extend(args)

    cwd = os.getcwd()

    try:
        status = -1
        oldargs = sys.argv
        sys.argv = sysargs

        oldout, olderr = sys.stdout, sys.stderr
        sys.stdout = StringIO()
        sys.stdout.name = "StringIO"
        sys.stderr = StringIO()

        if in_directory:
            os.chdir(in_directory)
        else:
            in_directory = cwd

        try:
            print('running:', scriptname, 'in:', in_directory, file=oldout)
            print('arguments', sysargs, file=oldout)

            status = _runscript(scriptname)
        except SystemExit as err:
            status = err.code
        except:  # pylint: disable=bare-except
            traceback.print_exc(file=sys.stderr)
            status = -1
    finally:
        sys.argv = oldargs
        out, err = sys.stdout.getvalue(), sys.stderr.getvalue()
        sys.stdout, sys.stderr = oldout, olderr

        os.chdir(cwd)

    if status != 0 and not fail_ok:
        print(out)
        print(err)
        assert False, (status, out, err)

    return status, out, err


def get_test_data(filename):
    filepath = None
    try:
        filepath = resource_filename(
            Requirement.parse("sourmash"), "sourmash/sourmash/test-data/"\
                + filename)
    except ResolutionError:
        pass
    if not filepath or not os.path.isfile(filepath):
        filepath = os.path.join(os.path.dirname(__file__), 'test-data',
                                filename)
    return filepath


class TempDirectory(object):
    def __init__(self):
        self.tempdir = tempfile.mkdtemp(prefix='sourmashtest_')

    def __enter__(self):
        return self.tempdir

    def __exit__(self, exc_type, exc_value, traceback):
        try:
            shutil.rmtree(self.tempdir, ignore_errors=True)
        except OSError:
            pass

        if exc_type:
            return False


def run_shell_cmd(cmd, fail_ok=False, in_directory=None):
    cwd = os.getcwd()
    if in_directory:
        os.chdir(in_directory)

    print('running: ', cmd)
    try:
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        (out, err) = proc.communicate()

        out = out.decode('utf-8')
        err = err.decode('utf-8')

        if proc.returncode != 0 and not fail_ok:
            print('out:', out)
            print('err:', err)
            raise AssertionError("exit code is non zero: %d" % proc.returncode)

        return (proc.returncode, out, err)
    finally:
        os.chdir(cwd)
