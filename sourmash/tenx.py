"""
10x-sequencing specific utility functions.
"""

import os
from .logging import notify
from collections import defaultdict
import tempfile
import time
import numpy as np

CELL_BARCODES = ['CB', 'XC']
UMIS = ['UB', 'XM']


def pass_alignment_qc(alignment, barcodes):
    """
    Check high quality mapping, QC-passing barcode and UMI of alignment.

    alignment :
        aligned bam segment
    barcodes : list
        List of cellular barcode strings
    Returns
    -------
    pass_qc : boolean
        true if a high quality, QC passing barcode with a UMI, false otherwise
    """
    high_quality_mapping = alignment.mapq == 255
    if barcodes is not None:
        good_cell_barcode = any(
            [alignment.has_tag(cb) and alignment.get_tag(cb) in barcodes
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
    Return a dictionary with cell barcode and the renamed barcode.

    barcodes : list
        List of cellular barcode strings
    barcode_renamer : str
        Path to tab-separated file mapping barcodes to their new name
        e.g. with channel or cell annotation label,
        e.g. AAATGCCCAAACTGCT-1    lung_epithelial_cell|AAATGCCCAAACTGCT-1
    Returns
    -------
    barcode_renamer : dict
        A (str, str) mapping of the original barcode to its new name
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
    """Read single-column barcodes.tsv and genes.tsv files from 10x.

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
        barcodes = np.unique([line.strip() for line in f])
    return barcodes


def read_bam_file(bam_path):
    """Read from a QC-pass bam file.

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
    import pysam

    return pysam.AlignmentFile(bam_path, mode='rb')


def shard_bam_file(bam_file_path, chunked_file_line_count, shards_folder):
    """Shard QC-pass bam file with the given line count and save to shards_folder

    Parameters
    ----------
    bam_file_path : str
        Bam file to shard
    chunked_file_line_count: int
        number of lines/alignment reads in each sharded bam file
    shards_folder: str
        absolute path to save the sharded bam files to
    Returns
    -------

    shards : list
        list of sharded bam filenames
    """
    import pysam

    notify("Sharding the bam file")
    startt = time.time()
    file_names = []

    with read_bam_file(bam_file_path) as bam_file:
        line_count = 0
        file_count = 0
        header = bam_file.header
        for index, alignment in enumerate(bam_file):
            if line_count == 0:
                file_name = os.path.join(
                    shards_folder,
                    "temp_bam_shard_{}.bam".format(file_count))
                file_names.append(file_name)
                outf = pysam.AlignmentFile(file_name, "wb", header=header)
            if line_count == chunked_file_line_count:
                file_count = file_count + 1
                line_count = 0
                outf.write(alignment)
                outf.close()
                notify("===== Sharding bam file ====== {}", file_count,
                       end="\r")
            else:
                outf.write(alignment)
                line_count = line_count + 1
        outf.close()

    notify("time taken to shard the bam file into {} shards is {:.5f} seconds",
           file_count, time.time() - startt)
    return file_names


def bam_to_fasta(barcodes, barcode_renamer, delimiter, umi_filter, bam_file):
    """Convert 10x bam to one-record-per-cell fasta.

    Parameters
    ----------
    bam : bamnostic.AlignmentFile
    barcodes : list of str
        QC-passing barcodes
    barcode_renamer : str or None
        Tab-separated filename mapping a barcode to a new name, e.g.
        AAATGCCCAAACTGCT-1    lung_epithelial_cell|AAATGCCCAAACTGCT-1
    delimiter : str
        Non-DNA or protein alphabet character to be ignored, e.g. if a cell
        has two sequences 'AAAAAAAAA' and 'CCCCCCCC', they would be
        concatenated as 'AAAAAAAAAXCCCCCCCC'.
    umi_filter : boolean
        If umi filter is True, then umi is written in place of fasta
            record name
        barcode is the fasta file name in the output fastas.
    Returns
    -------
    filenames: list
        one temp fasta filename for one cell's high-quality, non-duplicate
        reads

    """
    bam = read_bam_file(bam_file)

    # Filter out high quality alignments and/or alignments with selected
    # barcoddes
    bam_filtered = (x for x in bam if pass_alignment_qc(x, barcodes))
    if barcode_renamer is not None and barcodes is not None:
        renamer = parse_barcode_renamer(barcodes, barcode_renamer)
    else:
        renamer = None
    cell_sequences = defaultdict(str)

    for alignment in bam_filtered:
        # Get barcode of alignment, looks like "AAATGCCCAAACTGCT-1"
        # a bam file might have good cell barcode as any of the tags in
        # CELL_BARCODES
        for cb in CELL_BARCODES:
            if alignment.has_tag(cb):
                barcode = alignment.get_tag(cb)

        # If umi_filter is True, add the umi to the cell barcode separated by X
        # "AAATGCCCAAACTGCT-1" is the CELL_BARCODE and AGCT is the UMI
        # then the key in cell_sequences would be AAATGCCCAAACTGCT-1XAGCT
        renamed = renamer[barcode] if renamer is not None else barcode
        if umi_filter:
            for umi_tag in UMIS:
                if alignment.has_tag(umi_tag):
                    umi = alignment.get_tag(umi_tag)
            renamed = renamed + delimiter + umi

        # Make a long string of all the cell sequences, separated
        # by a non-alphabet letter
        cell_sequences[renamed] += alignment.seq + delimiter

    filenames = list(set(write_cell_sequences(cell_sequences, umi_filter, delimiter)))
    notify("bam_to_fasta conversion completed on {}", bam_file, end='\r', flush=True)

    bam.close()

    return filenames


def write_cell_sequences(cell_sequences, umi_filter=False, delimiter="X"):
    """
    Write each cell's sequences to an individual file
        Parameters
    ----------
    cell_sequences: dictionary with a cell and corresponding sequence
    if umi_filter is True, the cell key is expected to contain umi as well
    separated by the delimiter.
    e.g cell_key dictionary if umi_filter False is {AAAAAAAAA: AGCTACACTA},
    else {AAAAAAAAAXACTAG: AGCTACACTA} - In this case AAAAAAAAA would be cell
    barcode and ACTAG would be umi. The umi will be further used by downstream
    processing functions appropriately. The barcode is safely returned as the
    fasta filename and the umi is saved as record.name/sequence id in the
    fasta file
    umi_filter : boolean, default False
        If umi filter is True, then umi is written in the fasta file as
        record.name, cell barcode is the fasta file name in the output fastas.
    delimiter : str, default X
        Used to separate barcode and umi in the cell sequences dict.

    Returns
    -------
    filenames: generator
        one temp fasta filename for one cell/cell_umi with  sequence
    """
    temp_folder = tempfile.mkdtemp()

    for cell, seq in cell_sequences.items():
        if umi_filter:
            barcode, umi = cell.split(delimiter)
            filename = os.path.join(temp_folder, barcode + '.fasta')

            # Append to an existing barcode file with a different umi
            with open(filename, "a") as f:
                if os.stat(filename).st_size != 0:
                    f.write(">{}\n{}\n".format(umi, seq))
                else:
                    f.write(">{}\n{}".format(umi, seq))
            yield filename
        else:
            filename = os.path.join(temp_folder, cell + '.fasta')
            with open(filename, "w") as f:
                f.write(">{}\n{}".format(cell, seq))
            yield filename
