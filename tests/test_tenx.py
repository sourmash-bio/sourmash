from . import sourmash_tst_utils as utils
import sourmash.tenx as sourmash_tenx

import pysam as bs
import os


def test_read_barcodes_file():
    filename = utils.get_test_data('10x-example/barcodes.tsv')
    barcodes = sourmash_tenx.read_barcodes_file(filename)
    assert len(barcodes) == 10


def test_read_bam_file():
    filename = utils.get_test_data('10x-example/possorted_genome_bam.bam')
    bam_file = sourmash_tenx.read_bam_file(filename)
    assert isinstance(bam_file, bs.AlignmentFile)
    total_alignments = sum(1 for _ in bam_file)
    assert total_alignments == 1714


def test_shard_bam_file():
    filename = utils.get_test_data('10x-example/possorted_genome_bam.bam')
    length_alignment = 1714
    num_shards = 2
    bam_shard_files = sourmash_tenx.shard_bam_file(filename, length_alignment // num_shards)
    assert len(bam_shard_files) == 2
    for bam_file in bam_shard_files:
        if os.path.exists(bam_file):
            os.unlink(bam_file)
    bam_shard_files = sourmash_tenx.shard_bam_file(filename, length_alignment)
    assert len(bam_shard_files) == 1
    for bam_file in bam_shard_files:
        if os.path.exists(bam_file):
            os.unlink(bam_file)


def test_pass_alignment_qc():
    barcodes = sourmash_tenx.read_barcodes_file(
        utils.get_test_data('10x-example/barcodes.tsv'))
    bam = sourmash_tenx.read_bam_file(
        utils.get_test_data('10x-example/possorted_genome_bam.bam'))

    total_pass = sum(1 for alignment in bam if
                     sourmash_tenx.pass_alignment_qc(alignment, barcodes))
    assert total_pass == 439


def test_pass_alignment_qc_filtered():
    bam = sourmash_tenx.read_bam_file(
        utils.get_test_data('10x-example/possorted_genome_bam_filtered.bam'))
    total_pass = sum(1 for alignment in bam if
                     sourmash_tenx.pass_alignment_qc(alignment, None))
    assert total_pass == 192


def test_parse_barcode_renamer():
    filename = utils.get_test_data('10x-example/barcodes.tsv')
    barcodes = sourmash_tenx.read_barcodes_file(filename)
    renamer = sourmash_tenx.parse_barcode_renamer(barcodes, None)
    for key, value in renamer.items():
        assert key == value
    assert len(renamer) == len(barcodes)


def test_bam_to_cell_sequences():
    bam_file = utils.get_test_data('10x-example/possorted_genome_bam.bam')
    barcodes = sourmash_tenx.read_barcodes_file(
        utils.get_test_data('10x-example/barcodes.tsv'))
    cell_sequences = sourmash_tenx.bam_to_cell_sequences(
        barcodes, barcode_renamer=None, delimiter='X', bam_file=bam_file)
    assert len(cell_sequences) == 8


def test_filtered_bam_to_cell_sequences():
    bam_file = utils.get_test_data('10x-example/possorted_genome_bam_filtered.bam')
    cell_sequences = sourmash_tenx.bam_to_cell_sequences(
        barcodes=None, barcode_renamer=None, delimiter='X', bam_file=bam_file)
    assert len(cell_sequences) == 156
