import json

from jsonschema import validate
import pytest

from sourmash.schema import SBT_SCHEMAS, MINHASH_SCHEMAS

from . import sourmash_tst_utils as utils


@pytest.mark.parametrize("test_data_file", [
    ('v5.sbt.json'),
    ('v4.sbt.json'),
    ('47.fa.sig')
])
def test_validate_data(test_data_file):
    with open(utils.get_test_data(test_data_file), 'r') as f:
        instance = json.load(f)

    if test_data_file.endswith(".sig"):
        validate(instance=instance, schema=MINHASH_SCHEMAS["0.4"])
    elif test_data_file.endswith(".sbt.json"):
        validate(instance=instance, schema=SBT_SCHEMAS[instance["version"]])
    else:
        raise Exception("Couldn't decide schema to use")
