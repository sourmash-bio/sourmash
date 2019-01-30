from . import sourmash_tst_utils as utils

from sourmash.tenx import read_10x_folder, read_single_column, \
    _pass_alignment_qc, _parse_barcode_renamer, barcode_iterator


def test_read_single_column():
    filename = utils.get_test_data('10x-example/barcodes.tsv')
    barcodes = read_single_column(filename)
    assert len(barcodes) == 10


def test_read_10x_folder():
    import bamnostic as bs

    tenx_folder = utils.get_test_data('10x-example')

    barcodes, bam_file = read_10x_folder(tenx_folder)

    assert len(barcodes) == 10
    assert isinstance(bam_file, bs.AlignmentFile)

    total_alignments = sum(1 for _ in bam_file)
    assert total_alignments == 1714


def test__pass_alignment_qc():
    tenx_folder = utils.get_test_data('10x-example')

    barcodes, bam_file = read_10x_folder(tenx_folder)

    total_pass = sum(1 for alignment in bam_file if
                     _pass_alignment_qc(alignment, barcodes))
    assert total_pass == 1610


def test__parse_barcode_renamer():
    pass


def test_barcode_iterator():
    pass
