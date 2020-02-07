import os

from hypothesis import settings, Verbosity
import pytest


@pytest.fixture(params=[True, False])
def track_abundance(request):
    return request.param


@pytest.fixture(params=[True, False])
def dayhoff(request):
    return request.param


@pytest.fixture(params=[True, False])
def hp(request):
    return request.param


@pytest.fixture(params=[2, 5, 10])
def n_children(request):
    return request.param


# --- BEGIN - Only run tests using a particular fixture --- #
# Cribbed from: http://pythontesting.net/framework/pytest/pytest-run-tests-using-particular-fixture/
def pytest_collection_modifyitems(items, config):
    fixture_name = config.option.usesfixture
    if fixture_name is not None:
        selected_items = []
        deselected_items = []

        for item in items:
            if fixture_name in getattr(item, 'fixturenames', ()):
                selected_items.append(item)
            else:
                deselected_items.append(item)
        config.hook.pytest_deselected(items=deselected_items)
        items[:] = selected_items
# --- END - Only run tests using a particular fixture --- #

def pytest_addoption(parser):
    parser.addoption("--usesfixture",
                     action="store",
                     default=None,
                     help="just run tests that use a particular fixture")

    parser.addoption("--run-hypothesis", action="store_true",
                     help="run hypothesis tests")

def pytest_runtest_setup(item):
    hyp = any(mark for mark in item.iter_markers(name="hypothesis"))
    if hyp:
        if not item.config.getoption("--run-hypothesis"):
            pytest.skip("set --run-hypothesis option to run hypothesis tests")

settings.register_profile("ci", max_examples=1000)
settings.register_profile("dev", max_examples=10)
settings.register_profile("debug", max_examples=10, verbosity=Verbosity.verbose)
settings.load_profile(os.getenv(u'HYPOTHESIS_PROFILE', 'default'))
