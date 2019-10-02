"""
Functions implementing the 'compute' command and related functions.
"""
from __future__ import print_function, division, absolute_import

import argparse
import glob
import os
import os.path
import sys
import random
import itertools
import time
from collections import defaultdict, OrderedDict
from functools import partial
import numpy as np

import screed
from .sourmash_args import SourmashArgumentParser
from . import DEFAULT_SEED, MinHash
from . import signature as sig
from . import signature_json
from . import sourmash_args
from . import np_utils
from .logging import notify, error, set_quiet, debug

from .sourmash_args import DEFAULT_N
DEFAULT_COMPUTE_K = '21,31,51'
DEFAULT_LINE_COUNT = 1500


def compute(args):
    """Compute the signature for one or more files.

    Use cases:
        sourmash compute multiseq.fa              => multiseq.fa.sig, etc.
        sourmash compute genome.fa --singleton    => genome.fa.sig
        sourmash compute file1.fa file2.fa -o file.sig
            => creates one output file file.sig, with one signature for each
               input file.
        sourmash compute file1.fa file2.fa --merge merged -o file.sig
            => creates one output file file.sig, with all sequences from
               file1.fa and file2.fa combined into one signature.
    """
    parser = SourmashArgumentParser()
    parser.add_argument('filenames', nargs='+',
                        help='file(s) of sequences')

    sourmash_args.add_construct_moltype_args(parser)

    parser.add_argument('-q', '--quiet', action='store_true',
                        help='suppress non-error output')
    parser.add_argument('--input-is-protein', action='store_true',
                        help='Consume protein sequences - no translation needed.')
    parser.add_argument('-k', '--ksizes',
                        default=DEFAULT_COMPUTE_K,
                        help='comma-separated list of k-mer sizes (default: %(default)s)')
    parser.add_argument('-n', '--num-hashes', type=int,
                        default=DEFAULT_N,
                        help='number of hashes to use in each sketch (default: %(default)i)')
    parser.add_argument('--check-sequence', action='store_true',
                        help='complain if input sequence is invalid (default: False)')
    parser.add_argument('-f', '--force', action='store_true',
                        help='recompute signatures even if the file exists (default: False)')
    parser.add_argument('-o', '--output', type=argparse.FileType('wt'),
                        help='output computed signatures to this file')
    parser.add_argument('--singleton', action='store_true',
                        help='compute a signature for each sequence record individually (default: False)')
    parser.add_argument('--merge', '--name', type=str, default='', metavar="MERGED",
                        help="merge all input files into one signature named this")
    parser.add_argument('--name-from-first', action='store_true',
                        help="name the signature generated from each file after the first record in the file (default: False)")
    parser.add_argument('--input-is-10x', action='store_true',
                        help="Input is 10x single cell output folder (default: False)")
    parser.add_argument('--count-valid-reads', default=0, type=int,
                        help="For 10x input only (i.e input-is-10x flag is True), "
                        "A barcode is only considered a valid barcode read "
                        "and its signature is written if number of umis are greater "
                        "than count-valid-reads. It is used to weed out cell barcodes "
                        "with few umis that might have been due to false rna enzyme reactions")
    parser.add_argument('--write-barcode-meta-csv', type=str,
                        help="For 10x input only (i.e input-is-10x flag is True), for each of the unique barcodes, "
                        "Write to a given path, number of reads and number of umis per barcode.")
    parser.add_argument('-p', '--processes', default=2, type=int,
                        help='For 10x input only (i.e input-is-10x flag is True, '
                        'Number of processes to use for reading 10x bam file')
    parser.add_argument('--save-fastas', type=str,
                        help='For 10x input only (i.e input-is-10x flag is True), '
                        'save merged fastas for all the unique barcodes to {CELL_BARCODE}.fasta '
                        'save by default in the same directory from which the command is run')
    parser.add_argument('--line-count', type=int,
                        help='For 10x input only (i.e input-is-10x flag is True), line count for each bam shard',
                        default=DEFAULT_LINE_COUNT)
    parser.add_argument('--track-abundance', action='store_true',
                        help='track k-mer abundances in the generated signature (default: False)')
    parser.add_argument('--scaled', type=float, default=0,
                        help='choose number of hashes as 1 in FRACTION of input k-mers')
    parser.add_argument('--seed', type=int,
                        help='seed used by MurmurHash (default: 42)',
                        default=DEFAULT_SEED)
    parser.add_argument('--randomize', action='store_true',
                        help='shuffle the list of input filenames randomly')
    parser.add_argument('--license', default='CC0', type=str,
                        help='signature license. Currently only CC0 is supported.')
    parser.add_argument('--rename-10x-barcodes', type=str,
                        help="Tab-separated file mapping 10x barcode name "
                        "to new name, e.g. with channel or cell "
                        "annotation label", required=False)
    parser.add_argument('--barcodes-file', type=str,
                        help="Barcodes file if the input is unfiltered 10x bam file", required=False)

    args = parser.parse_args(args)
    set_quiet(args.quiet)

    if args.license != 'CC0':
        error('error: sourmash only supports CC0-licensed signatures. sorry!')
        sys.exit(-1)

    if args.input_is_protein and args.dna:
        notify('WARNING: input is protein, turning off nucleotide hashing')
        args.dna = False
        args.protein = True

    if args.scaled:
        if args.scaled < 1:
            error('ERROR: --scaled value must be >= 1')
            sys.exit(-1)
        if args.scaled != round(args.scaled, 0):
            error('ERROR: --scaled value must be integer value')
            sys.exit(-1)
        if args.scaled >= 1e9:
            notify('WARNING: scaled value is nonsensical!? Continuing anyway.')

        if args.num_hashes != 0:
            notify('setting num_hashes to 0 because --scaled is set')
            args.num_hashes = 0
 
    notify('computing signatures for files: {}', ", ".join(args.filenames))

    if args.randomize:
        notify('randomizing file list because of --randomize')
        random.shuffle(args.filenames)

    # get list of k-mer sizes for which to compute sketches
    ksizes = args.ksizes
    if ',' in ksizes:
        ksizes = ksizes.split(',')
        ksizes = list(map(int, ksizes))
    else:
        ksizes = [int(ksizes)]

    notify('Computing signature for ksizes: {}', str(ksizes))
    num_sigs = 0
    if args.dna and args.protein:
        notify('Computing both nucleotide and protein signatures.')
        num_sigs = 2*len(ksizes)
    elif args.dna and args.dayhoff:
        notify('Computing both nucleotide and Dayhoff-encoded protein '
               'signatures.')
        num_sigs = 2*len(ksizes)
    elif args.dna:
        notify('Computing only nucleotide (and not protein) signatures.')
        num_sigs = len(ksizes)
    elif args.protein:
        notify('Computing only protein (and not nucleotide) signatures.')
        num_sigs = len(ksizes)
    elif args.dayhoff:
        notify('Computing only Dayhoff-encoded protein (and not nucleotide) '
               'signatures.')
        num_sigs = len(ksizes)

    if (args.protein or args.dayhoff) and not args.input_is_protein:
        bad_ksizes = [ str(k) for k in ksizes if k % 3 != 0 ]
        if bad_ksizes:
            error('protein ksizes must be divisible by 3, sorry!')
            error('bad ksizes: {}', ", ".join(bad_ksizes))
            sys.exit(-1)

    notify('Computing a total of {} signature(s).', num_sigs)

    if num_sigs == 0:
        error('...nothing to calculate!? Exiting!')
        sys.exit(-1)

    if args.merge and not args.output:
        error("must specify -o with --merge")
        sys.exit(-1)

    def make_minhashes():
        seed = args.seed

        # one minhash for each ksize
        Elist = []
        for k in ksizes:
            if args.protein:
                E = MinHash(ksize=k, n=args.num_hashes,
                            is_protein=True,
                            dayhoff=False,
                            track_abundance=args.track_abundance,
                            scaled=args.scaled,
                            seed=seed)
                Elist.append(E)
            if args.dayhoff:
                E = MinHash(ksize=k, n=args.num_hashes,
                            is_protein=True,
                            dayhoff=True,
                            track_abundance=args.track_abundance,
                            scaled=args.scaled,
                            seed=seed)
                Elist.append(E)
            if args.dna:
                E = MinHash(ksize=k, n=args.num_hashes,
                            is_protein=False,
                            dayhoff=False,
                            track_abundance=args.track_abundance,
                            scaled=args.scaled,
                            seed=seed)
                Elist.append(E)
        return Elist

    def add_seq(Elist, seq, input_is_protein, check_sequence):
        for E in Elist:
            if input_is_protein:
                E.add_protein(seq)
            else:
                E.add_sequence(seq, not check_sequence)

    def build_siglist(Elist, filename, name=None):
        return [ sig.SourmashSignature(E, filename=filename,
                                       name=name) for E in Elist ]

    def iter_split(string, sep=None):
        """Split a string by the given separator and
        return the results in a generator"""
        sep = sep or ' '
        groups = itertools.groupby(string, lambda s: s != sep)
        return (''.join(g) for k, g in groups if k)

    def fasta_to_sig_record(index):
        """Convert fasta to sig record"""
        if umi_filter:
            return filtered_umi_to_sig(index)
        else:
            return unfiltered_umi_to_sig(index)

    def unfiltered_umi_to_sig(index):
        """Returns signature records across fasta files for a unique barcode"""
        # Initializing time
        startt = time.time()

        # Getting all fastas for a given barcode
        # from different shards
        single_barcode_fastas = all_fastas_sorted[index]

        count = 0
        # Iterating through fasta files for single barcode from different
        # fastas
        for fasta in iter_split(single_barcode_fastas, ","):

            # Initializing the fasta file to write
            # all the sequences from all bam shards to
            if count == 0:
                unique_fasta_file = os.path.basename(fasta)
                if args.save_fastas:
                    f = open(unique_fasta_file, "w")

            # Add sequence
            for record in screed.open(fasta):
                sequence = record.sequence
                add_seq(Elist, sequence,
                        args.input_is_protein, args.check_sequence)

                # Updating the fasta file with each of the sequences
                if args.save_fastas:
                    f.write(">{}\n{}".format(filename, record.sequence))

            # Delete fasta file in tmp folder
            if os.path.exists(fasta):
                os.unlink(fasta)

            count += 1

        # Close the opened fasta file
        if args.save_fastas:
            f.close()

        # Build signature records
        barcode_name = unique_fasta_file.replace(".fasta", "")
        siglist = build_siglist(
            Elist,
            os.path.join(args.filenames[0], unique_fasta_file),
            name=barcode_name)
        records = signature_json.add_meta_save(siglist)

        notify(
            "time taken to build signature records for a barcode {} is {:.5f} seconds",
            unique_fasta_file, time.time() - startt, end='\r')
        return records

    def filtered_umi_to_sig(index):
        """Returns signature records for all the fasta files for a unique 
        barcode, only if it has more than count-valid-reads number of umis."""

        # Initializing time
        startt = time.time()

        # Getting all fastas for a given barcode
        # from different shards
        single_barcode_fastas = all_fastas_sorted[index]

        debug("calculating umi counts", end="\r", flush=True)
        # Tracking UMI Counts
        umis = defaultdict(int)
        # Iterating through fasta files for single barcode from different
        # fastas
        for fasta in iter_split(single_barcode_fastas, ","):
            # calculate unique umi, sequence counts
            for record in screed.open(fasta):
                umis[record.name] += record.sequence.count(delimiter)

        if args.write_barcode_meta_csv:
            unique_fasta_file = os.path.basename(fasta)
            unique_fasta_file = unique_fasta_file.replace(".fasta", "_meta.txt")
            with open(unique_fasta_file, "w") as f:
                f.write("{} {}".format(len(umis), sum(list(umis.values()))))

        debug("Completed tracking umi counts", end="\r", flush=True)
        if len(umis) < args.count_valid_reads:
            return []
        count = 0
        for fasta in iter_split(single_barcode_fastas, ","):

            # Initializing fasta file to save the sequence to
            if count == 0:
                unique_fasta_file = os.path.basename(fasta)
                if args.save_fastas:
                    f = open(unique_fasta_file, "w")

            # Add sequences of barcodes with more than count-valid-reads umis
            for record in screed.open(fasta):
                sequence = record.sequence
                add_seq(Elist, sequence,
                        args.input_is_protein, args.check_sequence)

                # Updating the fasta file with each of the sequences
                if args.save_fastas:
                    f.write(">{}\n{}".format(filename, record.sequence))

            # Delete fasta file in tmp folder
            if os.path.exists(fasta):
                os.unlink(fasta)
            count += 1

        debug("Added sequences of unique barcode,umi to Elist", end="\r",
              flush=True)
        # Close the opened fasta file
        if args.save_fastas:
            f.close()
        # Build signature records
        barcode_name = unique_fasta_file.replace(".fasta", "")
        siglist = build_siglist(
            Elist,
            os.path.join(args.filenames[0], unique_fasta_file),
            name=barcode_name)
        records = signature_json.add_meta_save(siglist)
        notify(
            "time taken to build signature records for a barcode {} is {:.5f} seconds",
            unique_fasta_file, time.time() - startt, end="\r", flush=True)
        return records

    def save_siglist(siglist, output_fp, filename=None):
        # save!
        if output_fp:
            sigfile_name = args.output.name
            sig.save_signatures(siglist, args.output)
        else:
            if filename is None:
                raise Exception("internal error, filename is None")
            with open(filename, 'w') as fp:
                sigfile_name = filename
                sig.save_signatures(siglist, fp)
        notify(
            'saved signature(s) to {}. Note: signature license is CC0.',
            sigfile_name,
            end="\r")

    def calculate_chunksize(length, num_jobs):
        chunksize, extra = divmod(length, num_jobs)
        if extra:
            chunksize += 1
        return chunksize

    if args.input_is_10x:
        all_fastas_sorted = []
        delimiter = "X"
        umi_filter = True if args.count_valid_reads != 0 else False
        Elist = make_minhashes()
        CELL_BARCODE = "CELL_BARCODE"
        UMI_COUNT = "UMI_COUNT"
        READ_COUNT = "READ_COUNT"

    if args.track_abundance:
        notify('Tracking abundance of input k-mers.')

    if not args.merge:
        if args.output:
            siglist = []

        for filename in args.filenames:
            sigfile = os.path.basename(filename) + '.sig'
            if not args.output and os.path.exists(sigfile) and not \
                args.force:
                notify('skipping {} - already done', filename)
                continue

            if args.singleton:
                siglist = []
                for n, record in enumerate(screed.open(filename)):
                    # make minhashes for each sequence
                    Elist = make_minhashes()
                    add_seq(Elist, record.sequence,
                            args.input_is_protein, args.check_sequence)

                    siglist += build_siglist(Elist, filename, name=record.name)

                notify('calculated {} signatures for {} sequences in {}',
                       len(siglist), n + 1, filename)
            elif args.input_is_10x:
                from .tenx import (read_barcodes_file, shard_bam_file,
                                   bam_to_fasta)
                from pathos import multiprocessing

                # Initializing time
                startt = time.time()

                # Setting barcodes file, some 10x files don't have a filtered
                # barcode file
                if args.barcodes_file is not None:
                    barcodes = read_barcodes_file(args.barcodes_file)
                else:
                    barcodes = None

                # Shard bam file to smaller bam file
                notify('... reading bam file from {}', filename)
                n_jobs = args.processes
                filenames, mmap_file = np_utils.to_memmap(np.array(
                    shard_bam_file(filename, args.line_count, os.getcwd())))

                # Create a per-cell fasta generator of sequences
                # If the reads should be filtered by barcodes and umis
                # umis are saved in fasta file as record name and name of
                # the fasta file is the barcode
                func = partial(
                    bam_to_fasta,
                    barcodes,
                    args.rename_10x_barcodes,
                    delimiter,
                    umi_filter)

                length_sharded_bam_files = len(filenames)
                chunksize = calculate_chunksize(length_sharded_bam_files,
                                                n_jobs)
                pool = multiprocessing.Pool(processes=n_jobs)
                notify(
                    "multiprocessing pool processes {} and chunksize {} calculated", n_jobs, chunksize)
                # All the fastas are stored in a string instead of a list
                # This saves memory per element of the list by 8 bits
                # If we have unique barcodes in the order of 10^6 before
                # filtering that would result in a huge list if each barcode
                # is saved as a separate element, hence the string
                all_fastas = "," .join(itertools.chain(*(
                    pool.imap(
                        lambda x: func(x.encode('utf-8')), filenames, chunksize=chunksize))))
                pool.close()
                pool.join()

                # clean up the memmap and sharded intermediary bam files
                [os.unlink(file) for file in filenames if os.path.exists(file)]
                del filenames
                os.unlink(mmap_file)
                notify("Deleted intermediary bam and memmap files")

                # Build a dictionary with each unique barcode as key and
                # their fasta files from different shards
                fasta_files_dict = OrderedDict()
                for fasta in iter_split(all_fastas, ","):
                    barcode = os.path.basename(fasta).replace(".fasta", "")
                    value = fasta_files_dict.get(barcode, "")
                    fasta_files_dict[barcode] = value + fasta + ","

                # Cleaning up to retrieve memory from unused large variables
                del all_fastas
                notify("Created fasta_files_dict")

                # Find unique barcodes
                all_fastas_sorted = list(fasta_files_dict.values())
                del fasta_files_dict
                unique_barcodes = len(all_fastas_sorted)
                notify("Found {} unique barcodes", unique_barcodes)

                # For each barcode, add all their sequences, find
                # minhash and convert them to signature The below
                # fasta_to_sig_record also takes into consideration if
                # umi_filter is True and acts accordingly to create
                # records for barcodes with considerable number of
                # umis as provided in count-valid-reads This allows us
                # to save unique barcodes with valid reads and
                # significant umis, otherwise signature files save
                # every barcode even with one umi/one read resulting
                # in passing every alignment in bam file without file
                # i.e resulting in a signature file in GB, whereas the
                # expected .sig file should be in the order of MB

                pool = multiprocessing.Pool(processes=n_jobs)
                chunksize = calculate_chunksize(unique_barcodes, n_jobs)
                notify("Pooled {} and chunksize {} mapped", n_jobs, chunksize)

                records = list(itertools.chain(*pool.imap(
                    lambda index: fasta_to_sig_record(index),
                    range(unique_barcodes),
                    chunksize=chunksize)))

                pool.close()
                pool.join()

                if args.write_barcode_meta_csv:

                    barcodes_meta_txts = glob.glob("*_meta.txt")

                    with open(args.write_barcode_meta_csv, "w") as fp:
                        fp.write("{},{},{}".format(CELL_BARCODE, UMI_COUNT,
                                                   READ_COUNT))
                        fp.write('\n')
                        for barcode_meta_txt in barcodes_meta_txts:
                            with open(barcode_meta_txt, 'r') as f:
                                umi_count, read_count = f.readline().split()
                                umi_count = int(umi_count)
                                read_count = int(read_count)

                                barcode_name = barcode_meta_txt.replace('_meta.txt', '')
                                fp.write("{},{},{}\n".format(barcode_name,
                                                             umi_count,
                                                             read_count))

                            os.unlink(barcode_meta_txt)

                filtered_barcode_signatures = len(records)
                notify("Signature records created for {}",
                       filtered_barcode_signatures)

                # Save the signature records in .sig file
                if args.output is not None:
                    signature_json.write_records_to_json(records, args.output)
                else:
                    signature_json.write_records_to_json(records, open(sigfile, "wt"))
                del records

                notify("time taken to save {} signature records for 10x folder is {:.5f} seconds",
                       filtered_barcode_signatures, time.time() - startt)
            else:
                # make minhashes for the whole file
                Elist = make_minhashes()

                # consume & calculate signatures
                notify('... reading sequences from {}', filename)
                name = None
                for n, record in enumerate(screed.open(filename)):
                    if n % 10000 == 0:
                        if n:
                            notify('\r...{} {}', filename, n, end='')
                        elif args.name_from_first:
                            name = record.name

                    add_seq(Elist, record.sequence,
                            args.input_is_protein, args.check_sequence)

                notify('...{} {} sequences', filename, n, end='')

                sigs = build_siglist(Elist, filename, name)
                if args.output:
                    siglist += sigs
                else:
                    siglist = sigs

                notify('calculated {} signatures for {} sequences in {}',
                       len(sigs), n + 1, filename)

            if not args.output and not args.input_is_10x:
                save_siglist(siglist, args.output, sigfile)

        if args.output and not args.input_is_10x:
            save_siglist(siglist, args.output, sigfile)
    else:                             # single name specified - combine all
        # make minhashes for the whole file
        Elist = make_minhashes()

        total_seq = 0
        for filename in args.filenames:
            # consume & calculate signatures
            notify('... reading sequences from {}', filename)

            for n, record in enumerate(screed.open(filename)):
                if n % 10000 == 0 and n:
                    notify('\r... {} {}', filename, n, end='')

                add_seq(Elist, record.sequence,
                        args.input_is_protein, args.check_sequence)
            notify('... {} {} sequences', filename, n + 1)

            total_seq += n + 1

        siglist = build_siglist(Elist, filename, name=args.merge)
        notify('calculated {} signatures for {} sequences taken from {} files',
               len(siglist), total_seq, len(args.filenames))

        # at end, save!
        save_siglist(siglist, args.output)
