import pytest
import sourmash_tst_utils as utils


def test_failed_sourmash_exception(runtmp):
    # this is current behavior:
    with pytest.raises(ValueError):
        runtmp.sourmash('')

    # ...this is desired behavior:
    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('')
