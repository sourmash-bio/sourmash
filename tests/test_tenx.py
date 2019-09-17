import pytest
pytest.importorskip('pysam')

from . import sourmash_tst_utils as utils
import sourmash.tenx as sourmash_tenx
import pysam as bs


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
    bam_file = sourmash_tenx.read_bam_file(filename)
    assert isinstance(bam_file, bs.AlignmentFile)

    expected_alignments = sum(1 for _ in bam_file)
    with utils.TempDirectory() as location:
        bam_shard_files = sourmash_tenx.shard_bam_file(
            filename, expected_alignments, location)
        assert len(bam_shard_files) == 1

    num_shards = 2
    with utils.TempDirectory() as location:
        bam_shard_files = sourmash_tenx.shard_bam_file(
            filename, expected_alignments // num_shards, location)
        assert len(bam_shard_files) == 2

        total_alignments = 0
        for bam_file in bam_shard_files:
            total_alignments += sum(1 for _ in sourmash_tenx.read_bam_file(bam_file))
        assert total_alignments == expected_alignments

        whole_bam_file = sourmash_tenx.read_bam_file(filename)
        for bam_file in bam_shard_files:
            for line in sourmash_tenx.read_bam_file(bam_file):
                assert line == next(whole_bam_file)


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
    total_alignments = sum(1 for _ in bam)
    bam = sourmash_tenx.read_bam_file(
        utils.get_test_data('10x-example/possorted_genome_bam_filtered.bam'))
    assert total_alignments == 1500
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

    renamer = sourmash_tenx.parse_barcode_renamer(
        barcodes, utils.get_test_data('10x-example/barcodes_renamer.tsv'))
    for key, value in renamer.items():
        assert key in value
        assert "epithelial_cell" in value
    assert len(renamer) == len(barcodes)


def test_bam_to_fasta():
    filename = utils.get_test_data('10x-example/barcodes.tsv')
    bam_file = utils.get_test_data('10x-example/possorted_genome_bam.bam')
    barcodes = sourmash_tenx.read_barcodes_file(filename)

    fastas = sourmash_tenx.bam_to_fasta(
        barcodes=barcodes, barcode_renamer=None, delimiter="X", umi_filter=False, bam_file=bam_file)
    assert len(list(fastas)) == 8


def test_filtered_bam_to_fasta():
    bam_file = utils.get_test_data('10x-example/possorted_genome_bam_filtered.bam')
    fastas = sourmash_tenx.bam_to_fasta(
        barcodes=None, barcode_renamer=None, delimiter='X', umi_filter=False, bam_file=bam_file)
    assert len(list(fastas)) == 32


def test_filtered_bam_to_umi_fasta():
    bam_file = utils.get_test_data('10x-example/possorted_genome_bam_filtered.bam')
    fastas = sourmash_tenx.bam_to_fasta(
        barcodes=None, barcode_renamer=None, delimiter='X', umi_filter=True, bam_file=bam_file)
    assert len(list(fastas)) == 32


def test_write_sequences():
    cell_sequences = {'AAATGCCCAAACTGCT-1': "atgc", 'AAATGCCCAAAGTGCT-1': "gtga"}
    fastas = list(sourmash_tenx.write_cell_sequences(cell_sequences))
    assert len(fastas) == len(cell_sequences)
    for fasta in fastas:
        assert fasta.endswith(".fasta")


def test_write_sequences_umi():
    cell_sequences = {'AAATGCCCAXAACTGCT-1': "atgc", 'AAATGCCXCAAAGTGCT-1': "gtga", 'AAATGCCXCAAAGTGCT-2': "gtgc"}
    fastas = list(sourmash_tenx.write_cell_sequences(cell_sequences, True))
    assert len(fastas) == len(cell_sequences)
    for fasta in fastas:
        assert fasta.endswith(".fasta")

