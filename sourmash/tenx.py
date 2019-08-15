from collections import defaultdict
import os
from .logging import notify
import time
import tempfile
import warnings
import pysam

CELL_BARCODES = ['CB', 'XC']
UMIS = ['UB', 'XB']
BAM_FILENAME = 'possorted_genome_bam.bam'
BARCODES_TSV = 'barcodes.tsv'


def pass_alignment_qc(alignment, barcodes):
    """Assert high quality mapping, QC-passing barcode and UMI of alignment"""
    high_quality_mapping = alignment.mapq == 255
    if barcodes is not None:
        good_cell_barcode = any([alignment.has_tag(cb) and alignment.get_tag(cb) in barcodes
                                 for cb in CELL_BARCODES])
    else:
        good_cell_barcode = any([alignment.has_tag(cb) for cb in CELL_BARCODES])
    good_molecular_barcode = any([alignment.has_tag(umi) for umi in UMIS])
    not_duplicate = not alignment.is_duplicate

    pass_qc = high_quality_mapping and good_cell_barcode and \
              good_molecular_barcode and not_duplicate
    return pass_qc


def parse_barcode_renamer(barcodes, barcode_renamer):
    """
    Return a dictionary with cell barcode and the renamed barcode
    :param barcodes: 10x barcodes list
    :param barcode_renamer: Tab-separated file mapping
        10x barcode name to new name, e.g. with channel or cell "
        "annotation label"
    :return: barcode renamer dictionary
    """
    if barcode_renamer is not None:
        renamer = {}

        with open(barcode_renamer) as f:
            for line in f.readlines():
                barcode, renamed = line.split()
                assert barcode in barcodes
                renamer[barcode] = renamed
    else:
        renamer = dict(zip(barcodes, barcodes))
    return renamer


def read_barcodes_file(barcode_path):
    """Read single-column barcodes.tsv and genes.tsv files from 10x
    Parameters
    ----------
    barcode_path : str
        Name of a 10x 'barcodes.tsv' files
    Returns
    -------
    barcodes : list
        List of QC-passing barcodes from 'barcodes.tsv'
    """
    with open(barcode_path) as f:
        lines = set(line.strip() for line in f)
    return lines


def read_bam_file(bam_path):
    """Read from a QC-pass bam file
    Parameters
    ----------
    bam_path : str
        Name of a 10x cellranger bam file
        'possorted_genome_bam.bam'
    Returns
    -------
    bam_file : pysam.AlignmentFile
        Iterator over possorted_genome_bam.bam file
    """
    return pysam.AlignmentFile(bam_path, mode='rb')


def read_10x_folder(tenx_folder):
    """Get QC-pass barcodes, genes, and bam file from a 10x folder
    Parameters
    ----------
    folder : str
        Name of a 10x cellranger output folder containing
        'possorted_genome_bam.bam', 'possorted_genome_bam.bam.bai'
        and 'barcodes.tsv' files
    Returns
    -------
    barcodes : list
        List of QC-passing barcodes from 'barcodes.tsv'
    bam_file : pysam.AlignmentFile
        Iterator over possorted_genome_bam.bam file
    """
    barcodes = read_barcodes_file(os.path.join(tenx_folder, BARCODES_TSV))
    if len(barcodes) > 1e5:
        warnings.warn(
            "This barcode file contains {} total "
            "number of barcodes, which is far greater than "
            "typical single-cell experiments as of 2019. Counting "
            "min-hashes on this file will take >2TB of memory. "
            "Is this barcode list filtered by gene, read, or UMI "
            "count?".format(len(barcodes)))

    bam_file = read_bam_file(os.path.join(tenx_folder, BAM_FILENAME))

    return barcodes, bam_file


def write_bam_file(bam_file, bam_write_path):
    """Write QC-pass bam file to bam_write_path
    Parameters
    ----------
    bam_write_path : str
        Name of a 10x cellranger bam file to write to
    Returns
    -------
    bam_file : pysam.AlignmentFile
        Iterator over possorted_genome_bam.bam file
    """
    with pysam.AlignmentFile(
            bam_write_path,
            "wb",
            template=bam_file, header=bam_file.header) as outf:
        for index, alignment in enumerate(bam_file):
            outf.write(alignment)
    return read_bam_file(bam_write_path)


def shard_bam_file(bam_file_path, chunked_file_line_count):
    """Shard QC-pass bam file with the given line count and save them to tmp dir
    Parameters
    ----------
    bam_file_path : str
        Bam file to shard
    chunked_file_line_count: int
        number of lines/alignment reads in each sharded bam file
    Returns
    -------
    list of sharded bam files
    """
    notify("Sharding a bam file")
    startt = time.time()
    file_names = []
    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
        line_count = 0
        file_count = 0
        header = bam_file.header
        for alignment in bam_file:
            if line_count == 0:
                tmpfilename = tempfile.NamedTemporaryFile(suffix=".bam", delete=False)
                file_names.append(tmpfilename.name)
                outf = pysam.AlignmentFile(tmpfilename, "wb", header=header)
            if line_count == chunked_file_line_count:
                file_count = file_count + 1
                line_count = 0
                outf.close()
                tmpfilename.close()
                notify("===== Sharding bam file ====== {}".format(file_count), end="\r")
            else:
                outf.write(alignment)
                line_count = line_count + 1

    notify("time taken to tile the large bam file is {:.5f} seconds".format(time.time() - startt))
    return file_names


def bam_to_fasta(barcodes, barcode_renamer, delimiter, one_file_per_cell, bam_file):
    """Convert 10x bam to one-record-per-cell fasta
    Parameters
    ----------
    bam : bamnostic.AlignmentFile
    barcodes : list of str
        QC-passing barcodes
    barcode_renamer : str or None
        Tab-separated filename mapping a barcode to a new name, e.g.
        AAATGCCCAAACTGCT-1    lung_epithelial_cell|AAATGCCCAAACTGCT-1
    delimiter : str, default "X"
        Non-DNA or protein alphabet character to be ignored, e.g. if a cell
        has two sequences 'AAAAAAAAA' and 'CCCCCCCC', they would be
        concatenated as 'AAAAAAAAAXCCCCCCCC'.
    Returns
    -------
    filenames: list
        list of fasta file for each cell sequence
    """
    bam = read_bam_file(bam_file)

    # If bam file is filtered already by barcodes, barcodes is None
    if barcodes is not None:
        bam_filtered = (x for x in bam if pass_alignment_qc(x, barcodes))
        renamer = parse_barcode_renamer(barcodes, barcode_renamer)
    else:
        bam_filtered = bam
        renamer = None

    cell_sequences = defaultdict(str)
    for alignment in bam_filtered:
        # Get barcode of alignment, looks like "AAATGCCCAAACTGCT-1"
        # a bam file might have good cell barcode as any of the tags in CELL_BARCODES
        for cb in CELL_BARCODES:
            if alignment.has_tag(cb):
                barcode = alignment.get_tag(cb)
        renamed = renamer[barcode] if renamer is not None else barcode

        # Make a long string of all the cell sequences, separated
        # by a non-alphabet letter
        cell_sequences[renamed] += alignment.seq + delimiter

    filenames = write_cell_sequences(cell_sequences)
    return filenames


def write_cell_sequences(cell_sequences):
    """Write each cell sequence to a fasta file
    Parameters
    ----------
    cell_sequences : str
        Name of a 10x cellranger bam file
        'possorted_genome_bam.bam'
    Returns
    -------
    fasta_file : pysam.AlignmentFile
        Iterator over possorted_genome_bam.bam file
    """
    temp_folder = tempfile.mkdtemp()

    for cell, seq in cell_sequences.items():
        filename = os.path.join(temp_folder, cell + '.fasta')
        with open(filename, "w") as f:
            f.write(">{}\n{}".format(cell, seq))
        yield filename
