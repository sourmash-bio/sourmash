import json

from jsonschema import validate
import pytest

from sourmash.schema import SBT_SCHEMAS, MINHASH_SCHEMAS

from . import sourmash_tst_utils as utils


@pytest.mark.parametrize(
    "test_data_file,schema", [
    ('v5.sbt.json', SBT_SCHEMAS['5']),
    ('v4.sbt.json', SBT_SCHEMAS['5']),
    ('47.fa.sig', MINHASH_SCHEMAS['0.4']),
])
def test_validate_data(test_data_file, schema):
    with open(utils.get_test_data(test_data_file), 'r') as f:
        instance = json.load(f)

    # TODO: we can figure out version from each file and choose
    # schema here, instead of parametrizing?

    validate(instance=instance, schema=schema)
