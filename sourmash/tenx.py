import os
from .logging import notify
import time
import tempfile
import warnings
import pysam


BAM_FILENAME = 'possorted_genome_bam.bam'
BARCODES_TSV = 'barcodes.tsv'


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
