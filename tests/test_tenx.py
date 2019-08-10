from . import sourmash_tst_utils as utils
import sourmash.tenx as sourmash_tenx

import pysam as bs
import tempfile
import os


def test_read_barcodes_file():
    filename = utils.get_test_data('10x-example/barcodes.tsv')
    barcodes = sourmash_tenx.read_barcodes_file(filename)
    assert len(barcodes) == 10


def test_read_bam_file():
    filename = utils.get_test_data('10x-example/possorted_genome_bam.bam')
    bam_file = sourmash_tenx.read_bam_file(filename)
    assert isinstance(bam_file, bs.AlignmentFile)


def test_write_bam_file():
    filename = utils.get_test_data('10x-example/possorted_genome_bam.bam')
    read_bam_file = sourmash_tenx.read_bam_file(filename)
    temp_folder = tempfile.mkdtemp()
    temp_file_name = os.path.join(temp_folder, 'temp.bam')
    write_bam_file = sourmash_tenx.write_bam_file(read_bam_file, temp_file_name)
    assert isinstance(read_bam_file, bs.AlignmentFile)
    assert isinstance(write_bam_file, bs.AlignmentFile)
    assert read_bam_file.header.to_dict() == write_bam_file.header.to_dict()
    for alignment1, alignment2 in zip(read_bam_file, write_bam_file):
        assert alignment1 == alignment2


def test_read_10x_folder():

    tenx_folder = utils.get_test_data('10x-example')

    barcodes, bam_file = sourmash_tenx.read_10x_folder(tenx_folder)

    assert len(barcodes) == 10
    assert isinstance(bam_file, bs.AlignmentFile)

    total_alignments = sum(1 for _ in bam_file)
    assert total_alignments == 1714


def test_tile():
    filename = utils.get_test_data('10x-example/possorted_genome_bam.bam')
    bam_tile_files = sourmash_tenx.tile(filename, 1714 // 2)
    assert len(bam_tile_files) == 2
