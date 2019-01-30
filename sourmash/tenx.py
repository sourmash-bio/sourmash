from collections import defaultdict
import os

CELL_BARCODE = 'CB'
UMI = 'UB'

BAM_FILENAME = 'possorted_genome_bam.bam'
BARCODES_TSV = 'barcodes.tsv'


def read_single_column(filename):
    """Read single-column barcodes.tsv and genes.tsv files from 10x"""
    with open(filename) as f:
        lines = set(line.strip() for line in f)
    return lines


def read_10x_folder(folder):
    """Get QC-pass barcodes, genes, and bam file from a 10x folder

    Parameters
    ----------
    folder : str
        Name of a 10x cellranger output folder containing
        'possorted_genome_bam.bam' and 'barcodes.tsv' files

    Returns
    -------
    barcodes : list
        List of QC-passing barcodes from 'barcodes.tsv'
    bam_file : bamnostic.AlignmentFile
        Iterator over possorted_genome_bam.bam file
    """
    import bamnostic as bs

    barcodes = read_single_column(os.path.join(folder, BARCODES_TSV))

    bam_file = bs.AlignmentFile(os.path.join(folder, BAM_FILENAME), mode='rb')

    return barcodes, bam_file


def _pass_alignment_qc(alignment, barcodes):
    """Assert high quality mapping, QC-passing barcode and UMI of alignment"""
    high_quality_mapping = alignment.mapq == 255
    good_barcode = CELL_BARCODE in alignment.tags and \
                   alignment.get_tag(CELL_BARCODE) in barcodes
    good_umi = UMI in alignment.tags

    pass_qc = high_quality_mapping and good_barcode and \
              good_umi
    return pass_qc


def _parse_barcode_renamer(barcodes, barcode_renamer):
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


def barcode_iterator(bam, barcodes, barcode_renamer):
    """Yield a (barcode, list of str) pair for each QC-pass barcode"""
    bam_filtered = (x for x in bam if _pass_alignment_qc(x, barcodes))

    renamer = _parse_barcode_renamer(barcodes, barcode_renamer)

    # alignments only have a CELL_BARCODE tag if they past QC
    bam_sort_by_barcode = sorted(bam_filtered,
                                 key=lambda x: x.get_tag(CELL_BARCODE))

    previous_barcode = None
    barcode_alignments = []
    for alignment in bam_sort_by_barcode:
        # Get barcode of alignment, looks like "AAATGCCCAAACTGCT-1"
        barcode = alignment.get_tag(CELL_BARCODE)

        # If this is a new non-null barcode, return all previous sequences
        if previous_barcode is not None and barcode != previous_barcode:
            yield renamer[previous_barcode], barcode_alignments

            # Reset the barcode alignments
            barcode_alignments = []

        # Add only the aligned sequence to this list of barcode alignments
        barcode_alignments.append(alignment.seq)

        # Set this current barcode as the previous one
        previous_barcode = barcode

    # Yield the final one
    yield renamer[previous_barcode], barcode_alignments


def bam_to_fasta(bam, barcodes, barcode_renamer, delimiter="X",
                 one_file_per_cell=False):
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
        Non-DNA or protein alphabet character to be ignored

    Returns
    -------

    """

    bam_filtered = (x for x in bam if _pass_alignment_qc(x, barcodes))

    renamer = _parse_barcode_renamer(barcodes, barcode_renamer)

    cell_sequences = defaultdict(str)

    for alignment in bam_filtered:

        if _pass_alignment_qc(alignment, barcodes):
            # Get barcode of alignment, looks like "AAATGCCCAAACTGCT-1"
            barcode = alignment.get_tag(CELL_BARCODE)
            renamed = renamer[barcode]

            # Make a long string of all the cell sequences, separated
            # by a non-alphabet letter
            cell_sequences[renamed] += alignment.seq + delimiter

    write_cell_sequences(cell_sequences, output_folder, one_file_per_cell)


def _write_one_file_per_cell(cell_sequences, output_folder):
    filename = os.path.join(output_folder,
                            "per_cell_aligned_sequences.fasta")

    with open(filename, "w") as f:
        for cell, seq in cell_sequences.items():
            f.write(f">{cell}\n{seq}")
    return filename


def _write_all_cells_in_one_file(cell_sequences, output_folder):
    output_folder = os.path.join(output_folder,
                                 "per_cell_aligned_sequences")
    os.makedirs(output_folder, exist_ok=True)

    filenames = []

    for cell, seq in cell_sequences.items():
        filename = os.path.join(output_folder, cell + '.fasta')
        with open(filename, "w"):
            f.write(f">{cell}\n{seq}")
        filenames.append(filename)
    return filenames


def write_cell_sequences(cell_sequences, output_folder, one_file_per_cell=False):
    if one_file_per_cell:
        filenames = _write_one_file_per_cell(cell_sequences, output_folder)
    else:
        filename = _write_all_cells_in_one_file(cell_sequences, output_folder)
        filenames = [filename]
    return filenames

