import pytest


@pytest.fixture(params=[True, False])
def track_abundance(request):
    return request.param
