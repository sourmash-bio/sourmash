import pytest
import sourmash_tst_utils as utils
from sourmash_tst_utils import SourmashCommandFailed


def test_failed_sourmash_exception(runtmp):
    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('')
