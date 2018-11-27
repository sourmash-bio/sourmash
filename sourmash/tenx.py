import os

CELL_BARCODE = 'CB'


def read_single_column(filename):
    """Read single-column barcodes.tsv and genes.tsv files from 10x"""
    with open(filename) as f:
        lines = set(line.strip() for line in f)
    return lines


def read_10x_folder(folder):
    """Get QC-pass barcodes, genes, and bam file from a 10x folder"""
    import bamnostic as bs

    barcodes = read_single_column(os.path.join(folder, 'barcodes.tsv'))

    bam_file = bs.AlignmentFile(
        os.path.join(folder, 'possorted_genome_bam.bam'), mode='rb')

    return barcodes, bam_file

def _pass_alignment_qc(alignment, barcodes):
    high_quality_mapping = alignment.mapq == 255
    good_barcode = 'CB' in alignment.tags and \
                   alignment.get_tag('CB') in barcodes
    good_umi = 'UB' in alignment.tags

    pass_qc = high_quality_mapping and good_barcode and \
              good_umi
    return pass_qc


def barcode_iterator(bam):
    bam_filtered = (x for x in bam if _pass_alignment_qc(x))

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
            yield previous_barcode, barcode_alignments

            # Reset the barcode alignments
            barcode_alignments = []

        # Add only the aligned sequence to this list of barcode alignments
        barcode_alignments.append(alignment.seq)

        # Set this current barcode as the previous one
        previous_barcode = barcode

    # Yield the final one
    yield previous_barcode, barcode_alignments
