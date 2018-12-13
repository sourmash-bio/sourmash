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
    """Get QC-pass barcodes, genes, and bam file from a 10x folder"""
    import bamnostic as bs

    barcodes = read_single_column(os.path.join(folder, BARCODES_TSV))

    bam_file = bs.AlignmentFile(os.path.join(folder, BAM_FILENAME), mode='rb')

    return barcodes, bam_file


def _pass_alignment_qc(alignment, barcodes):
    high_quality_mapping = alignment.mapq == 255
    good_barcode = CELL_BARCODE in alignment.tags and \
                   alignment.get_tag(CELL_BARCODE) in barcodes
    good_umi = UMI in alignment.tags

    pass_qc = high_quality_mapping and good_barcode and \
              good_umi
    return pass_qc


def _parse_barcode_renamer(barcodes, barcode_renamer):
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
