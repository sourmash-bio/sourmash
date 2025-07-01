import os
import sys

from hypothesis import settings, Verbosity
import pytest

import matplotlib.pyplot as plt

plt.rcParams.update({"figure.max_open_warning": 0})

from sourmash_tst_utils import TempDirectory, RunnerContext

sys.stdout = sys.stderr


# behavior is default behavior, present in both sourmash v4 and sourmash v5.
@pytest.fixture(params=["v4", "v5", "(default)"])
def cli_v4_and_v5(request):
    return request.param


# behavior is default behavior in v4, and maybe will be invoked by --v4 in v5.
@pytest.fixture(params=["v4", "(default)"])
def cli_v4_only(request):
    return request.param


# behavior is available with --v5 and will be default in sourmash v5.
@pytest.fixture(params=["v5"])
def cli_v5_only(request):
    return request.param


@pytest.fixture
def runtmp():
    with TempDirectory() as location:
        yield RunnerContext(location)


@pytest.fixture(scope="session")
def runtmp_session():
    with TempDirectory() as location:
        yield RunnerContext(location)


@pytest.fixture
def run():
    yield RunnerContext(os.getcwd())


@pytest.fixture(params=[True, False])
def track_abundance(request):
    return request.param


@pytest.fixture(params=[True, False])
def dayhoff(request):
    return request.param


@pytest.fixture(params=[True, False])
def hp(request):
    return request.param


@pytest.fixture(params=[True, False])
def keep_identifiers(request):
    return request.param


@pytest.fixture(params=[True, False])
def keep_versions(request):
    return request.param


@pytest.fixture(params=[2, 5, 10])
def n_children(request):
    return request.param


@pytest.fixture(params=["--linear", "--no-linear"])
def linear_gather(request):
    return request.param


@pytest.fixture(params=["--prefetch", "--no-prefetch"])
def prefetch_gather(request):
    return request.param


@pytest.fixture(params=["SBT", "rocksdb", "zip"])
def disk_index_type(request):
    return request.param


@pytest.fixture(params=[True, False])
def use_manifest(request):
    return request.param


@pytest.fixture(params=["json", "sql"])
def lca_db_format(request):
    return request.param


@pytest.fixture(params=["csv", "sql"])
def manifest_db_format(request):
    return request.param


@pytest.fixture(params=["sig", "sig.gz", "zip", ".d/", ".sqldb"])
def sig_save_extension(request):
    return request.param


@pytest.fixture(params=["sig", "sig.gz", "zip", ".d/"])
def sig_save_extension_abund(request):
    return request.param


# these should both always succeed for 'sig check' and 'sig collect' output
# manifests.
@pytest.fixture(params=["--abspath", "--relpath"])
def abspath_or_relpath(request):
    return request.param


# this will fail if subdirs used; see #3008. but this ensures v4 behavior of
# sig collect/sig check works, where manifest paths are interpreted relative
# to cwd.
@pytest.fixture(params=["--no-abspath", "--abspath", "--relpath"])
def abspath_relpath_v4(request):
    return request.param


# --- BEGIN - Only run tests using a particular fixture --- #
# Cribbed from: http://pythontesting.net/framework/pytest/pytest-run-tests-using-particular-fixture/
def pytest_collection_modifyitems(items, config):
    fixture_name = config.option.usesfixture
    if fixture_name is not None:
        selected_items = []
        deselected_items = []

        for item in items:
            if fixture_name in getattr(item, "fixturenames", ()):
                selected_items.append(item)
            else:
                deselected_items.append(item)
        config.hook.pytest_deselected(items=deselected_items)
        items[:] = selected_items


# --- END - Only run tests using a particular fixture --- #


def pytest_addoption(parser):
    parser.addoption(
        "--usesfixture",
        action="store",
        default=None,
        help="just run tests that use a particular fixture",
    )

    parser.addoption(
        "--run-hypothesis", action="store_true", help="run hypothesis tests"
    )


def pytest_runtest_setup(item):
    if item.config.getoption("--run-hypothesis"):
        if not any(mark for mark in item.iter_markers(name="hypothesis")):
            pytest.skip("--run-hypothesis option set, running only hypothesis tests")


settings.register_profile("ci", max_examples=1000)
settings.register_profile("dev", max_examples=10)
settings.register_profile("debug", max_examples=10, verbosity=Verbosity.verbose)
settings.load_profile(os.getenv("HYPOTHESIS_PROFILE", "default"))
