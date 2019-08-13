from collections import defaultdict
import os
from .logging import notify
import time
import tempfile
import warnings
import pysam

CELL_BARCODE = 'CB'
UMI = 'UB'
BAM_FILENAME = 'possorted_genome_bam.bam'
BARCODES_TSV = 'barcodes.tsv'


def pass_alignment_qc(alignment, barcodes):
    """Assert high quality mapping, QC-passing barcode and UMI of alignment"""
    high_quality_mapping = alignment.mapq == 255
    good_cell_barcode = alignment.has_tag(CELL_BARCODE) and \
                   alignment.get_tag(CELL_BARCODE) in barcodes
    good_molecular_barcode = alignment.has_tag(UMI)
    not_duplicate = not alignment.is_duplicate

    pass_qc = high_quality_mapping and good_cell_barcode and \
              good_molecular_barcode and not_duplicate
    return pass_qc


def parse_barcode_renamer(barcodes, barcode_renamer):
    """
    :param barcodes:
    :param barcode_renamer:
    :return:
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


def write_bam_file(bam_file, bam_write_path, line_count=None):
    """Write a part of QC-pass bam file to bam_write_path
    Parameters
    ----------
    bam_write_path : str
        Name of a 10x cellranger bam file
        'possorted_genome_bam.bam'
    Returns
    -------
    bam_file : pysam.AlignmentFile
        Iterator over possorted_genome_bam.bam file
    """
    line_count = 0
    with pysam.AlignmentFile(bam_write_path, "wb", template=bam_file, header=bam_file.header) as outf:
        for index, alignment in enumerate(bam_file):
            if line_count is not None:
                if index == line_count:
                    outf.write(alignment)
            else:
                outf.write(alignment)
    return read_bam_file(bam_write_path)


def tile(bam_file_path, chunked_file_line_count):
    notify("Tiling a bam file")
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
                notify("===== Tiling bam file ====== {}".format(file_count), end="\r")
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
    """
    bam = read_bam_file(bam_file)
    bam_filtered = (x for x in bam if pass_alignment_qc(x, barcodes))

    renamer = parse_barcode_renamer(barcodes, barcode_renamer)

    cell_sequences = defaultdict(str)
    for alignment in bam_filtered:
        # Get barcode of alignment, looks like "AAATGCCCAAACTGCT-1"
        barcode = alignment.get_tag(CELL_BARCODE)
        renamed = renamer[barcode]

        # Make a long string of all the cell sequences, separated
        # by a non-alphabet letter
        cell_sequences[renamed] += alignment.seq + delimiter

    filenames = write_cell_sequences(cell_sequences)
    return filenames


def write_cell_sequences(cell_sequences):
    """Write each cell's sequences to an individual file"""
    temp_folder = tempfile.mkdtemp()

    for cell, seq in cell_sequences.items():
        filename = os.path.join(temp_folder, cell + '.fasta')
        with open(filename, "w") as f:
            f.write(">{}\n{}".format(cell, seq))
        yield filename
