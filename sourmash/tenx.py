import os


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
