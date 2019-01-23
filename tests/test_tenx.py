from . import sourmash_tst_utils as utils

from sourmash.tenx import read_10x_folder, read_single_column, \
    _pass_alignment_qc, _parse_barcode_renamer, barcode_iterator


def test_read_single_column():
    filename = utils.get_test_data('10x-example/barcodes.tsv')
    barcodes = read_single_column(filename)
    assert len(barcodes) == 625


def test_read_10x_folder():
    utils.get_test_data('10x-example')