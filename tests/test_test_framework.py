import pytest
import sourmash_tst_utils as utils


def test_failed_sourmash_exception(runtmp):
    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('')
