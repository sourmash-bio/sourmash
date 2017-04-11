import pytest


@pytest.fixture(params=[True, False])
def track_abundance(request):
    return request.param


@pytest.fixture(params=[2, 5, 10])
def n_children(request):
    return request.param
