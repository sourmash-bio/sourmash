from __future__ import print_function, division

import argparse
import csv
import os
import os.path
import sys
from collections import namedtuple
import random

import screed
import sourmash_lib
from . import signature as sig
from . import sourmash_args
from .logging import notify, error, print_results, set_quiet

from .sourmash_args import DEFAULT_LOAD_K
DEFAULT_COMPUTE_K = '21,31,51'

DEFAULT_N = 500
WATERMARK_SIZE = 10000


def info(args):
    "Report sourmash version + version of installed dependencies."
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='report versions of khmer and screed')
    args = parser.parse_args(args)

    from . import VERSION
    notify('sourmash version {}', VERSION)
    notify('- loaded from path: {}', os.path.dirname(__file__))
    notify('')

    if args.verbose:
        import khmer
        notify('khmer version {}', khmer.__version__)
        notify('- loaded from path: {}', os.path.dirname(khmer.__file__))
        notify('')

        import screed
        notify('screed version {}', screed.__version__)
        notify('- loaded from path: {}', os.path.dirname(screed.__file__))


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
    parser = argparse.ArgumentParser()
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
    parser.add_argument('--track-abundance', action='store_true',
                        help='track k-mer abundances in the generated signature (default: False)')
    parser.add_argument('--scaled', type=float, default=0,
                        help='choose number of hashes as 1 in FRACTION of input k-mers')
    parser.add_argument('--seed', type=int,
                        help='seed used by MurmurHash (default: 42)',
                        default=sourmash_lib.DEFAULT_SEED)
    parser.add_argument('--randomize', action='store_true',
                        help='shuffle the list of input filenames randomly')
    parser.add_argument('--license', default='CC0', type=str,
                        help='signature license. Currently only CC0 is supported.')

    args = parser.parse_args(args)
    set_quiet(args.quiet)

    if args.license != 'CC0':
        error('error: sourmash only supports CC0-licensed signatures. sorry!')
        sys.exit(-1)

    if args.input_is_protein and args.dna:
        notify('WARNING: input is protein, turning off DNA hashing')
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
        notify('Computing both DNA and protein signatures.')
        num_sigs = 2*len(ksizes)
    elif args.dna:
        notify('Computing only DNA (and not protein) signatures.')
        num_sigs = len(ksizes)
    elif args.protein:
        notify('Computing only protein (and not DNA) signatures.')
        num_sigs = len(ksizes)

    if args.protein:
        bad_ksizes = [ str(k) for k in ksizes if k % 3 != 0 ]
        if bad_ksizes:
            error('protein ksizes must be divisible by 3, sorry!')
            error('bad ksizes: {}', ", ".join(bad_ksizes))
            sys.exit(-1)

    notify('Computing a total of {} signatures.', num_sigs)

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
                E = sourmash_lib.MinHash(ksize=k, n=args.num_hashes,
                                            is_protein=True,
                                    track_abundance=args.track_abundance,
                                            scaled=args.scaled,
                                            seed=seed)
                Elist.append(E)
            if args.dna:
                E = sourmash_lib.MinHash(ksize=k, n=args.num_hashes,
                                            is_protein=False,
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

    def save_siglist(siglist, output_fp, filename=None):
        # save!
        if output_fp:
            sig.save_signatures(siglist, args.output)
        else:
            if filename is None:
                raise Exception("internal error, filename is None")
            with open(filename, 'w') as fp:
                sig.save_signatures(siglist, fp)
        notify('saved {} signature(s). Note: signature license is CC0.'.format(len(siglist)))

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

                notify('calculated {} signatures for {} sequences in {}'.\
                          format(len(siglist), n + 1, filename))
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

                    s = record.sequence
                    add_seq(Elist, record.sequence,
                            args.input_is_protein, args.check_sequence)
                notify('')

                sigs = build_siglist(Elist, filename, name)
                if args.output:
                    siglist += sigs
                else:
                    siglist = sigs

                notify('calculated {} signatures for {} sequences in {}'.\
                          format(len(siglist), n + 1, filename))

            if not args.output:
                save_siglist(siglist, args.output, sigfile)

        if args.output:
            save_siglist(siglist, args.output, sigfile)
    else:                             # single name specified - combine all
        # make minhashes for the whole file
        Elist = make_minhashes()

        for filename in args.filenames:
            # consume & calculate signatures
            notify('... reading sequences from {}', filename)
            for n, record in enumerate(screed.open(filename)):
                if n % 10000 == 0 and n:
                    notify('\r... {} {}', filename, n, end='')

                add_seq(Elist, record.sequence,
                        args.input_is_protein, args.check_sequence)

        siglist = build_siglist(Elist, filename, name=args.merge)
        notify('calculated {} signatures for {} sequences taken from {}'.\
               format(len(siglist), n + 1, " ".join(args.filenames)))
        # at end, save!
        save_siglist(siglist, args.output)


def compare(args):
    "Compare multiple signature files and create a distance matrix."
    import numpy

    parser = argparse.ArgumentParser()
    parser.add_argument('signatures', nargs='+', help='list of signatures')
    parser.add_argument('-o', '--output')
    parser.add_argument('--ignore-abundance', action='store_true',
                        help='do NOT use k-mer abundances if present')
    sourmash_args.add_ksize_arg(parser, DEFAULT_LOAD_K)
    sourmash_args.add_moltype_args(parser)
    parser.add_argument('--csv', type=argparse.FileType('w'),
                        help='save matrix in CSV format (with column headers)')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='suppress non-error output')
    args = parser.parse_args(args)
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    # load in the various signatures
    siglist = []
    ksizes = set()
    moltypes = set()
    for filename in args.signatures:
        notify('loading {}', filename, end='\r')
        loaded = sig.load_signatures(filename, ksize=args.ksize,
                                     select_moltype=moltype)
        loaded = list(loaded)
        if not loaded:
            notify('\nwarning: no signatures loaded at given ksize/molecule type from {}', filename)
        siglist.extend(loaded)

        # track ksizes/moltypes
        for s in loaded:
            ksizes.add(s.minhash.ksize)
            moltypes.add(sourmash_args.get_moltype(s))

        # error out while loading if we have more than one ksize/moltype
        if len(ksizes) > 1 or len(moltypes) > 1:
            break

    # check ksizes and type
    if len(ksizes) > 1:
        error('multiple k-mer sizes loaded; please specify one with -k.')
        ksizes = sorted(ksizes)
        error('(saw k-mer sizes {})'.format(', '.join(map(str, ksizes))))
        sys.exit(-1)

    if len(moltypes) > 1:
        error('multiple molecule types loaded; please specify --dna, --protein')
        sys.exit(-1)

    notify(' '*79, end='\r')
    notify('loaded {} signatures total.'.format(len(siglist)))

    # check to make sure they're potentially compatible - either using
    # max_hash/scaled, or not.
    scaled_sigs = [s.minhash.max_hash for s in siglist]
    is_scaled = all(scaled_sigs)
    is_scaled_2 = any(scaled_sigs)

    # complain if it's not all one or the other
    if is_scaled != is_scaled_2:
        error('cannot mix scaled signatures with bounded signatures')
        sys.exit(-1)

    # if using --scaled, downsample appropriately
    if is_scaled:
        max_scaled = max(s.minhash.scaled for s in siglist)
        notify('downsampling to scaled value of {}'.format(max_scaled))
        for s in siglist:
            s.minhash = s.minhash.downsample_scaled(max_scaled)

    if len(siglist) == 0:
        error('no signatures!')
        sys.exit(-1)

    notify('')

    # build the distance matrix
    D = numpy.zeros([len(siglist), len(siglist)])
    numpy.set_printoptions(precision=3, suppress=True)

    # do all-by-all calculation
    labeltext = []
    for i, E in enumerate(siglist):
        for j, E2 in enumerate(siglist):
            if i < j:
                continue
            similarity = E.similarity(E2, args.ignore_abundance)
            D[i][j] = similarity
            D[j][i] = similarity

        if len(siglist) < 30:
            # for small matrices, pretty-print some output
            name_num = '{}-{}'.format(i, E.name())
            if len(name_num) > 20:
                name_num = name_num[:17] + '...'
            print_results('{:20s}\t{}'.format(name_num, D[i, :, ],))

        labeltext.append(E.name())

    print_results('min similarity in matrix: {:.3f}', numpy.min(D))

    # shall we output a matrix?
    if args.output:
        labeloutname = args.output + '.labels.txt'
        notify('saving labels to: {}', labeloutname)
        with open(labeloutname, 'w') as fp:
            fp.write("\n".join(labeltext))

        notify('saving distance matrix to: {}', args.output)
        with open(args.output, 'wb') as fp:
            numpy.save(fp, D)

    # output CSV?
    if args.csv:
        w = csv.writer(args.csv)
        w.writerow(labeltext)

        for i in range(len(labeltext)):
            y = []
            for j in range(len(labeltext)):
                y.append('{}'.format(D[i][j]))
            args.csv.write(','.join(y) + '\n')


def plot(args):
    "Produce a clustering and plot."
    import matplotlib as mpl
    mpl.use('Agg')
    import numpy
    import scipy
    import pylab
    import scipy.cluster.hierarchy as sch
    from . import fig as sourmash_fig

    # set up cmd line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('distances', help="output from 'sourmash compare'")
    parser.add_argument('--pdf', action='store_true',
                        help='output PDF, not PNG.')
    parser.add_argument('--labels', action='store_true',
                        help='show sample labels on dendrogram/matrix')
    parser.add_argument('--indices', action='store_false',
                        help='show sample indices but not labels')
    parser.add_argument('--vmax', default=1.0, type=float,
                        help='upper limit of heatmap scale; (default: %(default)f)')
    parser.add_argument('--vmin', default=0.0, type=float,
                        help='lower limit of heatmap scale; (default: %(default)f)')
    parser.add_argument("--subsample", type=int,
                        help="randomly downsample to this many samples, max.")
    parser.add_argument("--subsample-seed", type=int, default=1,
                        help="random seed for --subsample; default=1")
    args = parser.parse_args(args)

    # load files
    D_filename = args.distances
    labelfilename = D_filename + '.labels.txt'

    D = numpy.load(open(D_filename, 'rb'))
    labeltext = [ x.strip() for x in open(labelfilename) ]

    # build filenames, decide on PDF/PNG output
    dendrogram_out = os.path.basename(D_filename) + '.dendro'
    if args.pdf:
        dendrogram_out += '.pdf'
    else:
        dendrogram_out += '.png'

    matrix_out = os.path.basename(D_filename) + '.matrix'
    if args.pdf:
        matrix_out += '.pdf'
    else:
        matrix_out += '.png'

    ### make the dendrogram:
    fig = pylab.figure(figsize=(8,5))
    ax1 = fig.add_axes([0.1, 0.1, 0.7, 0.8])
    ax1.set_xticks([])
    ax1.set_yticks([])

    # subsample?
    if args.subsample:
        numpy.random.seed(args.subsample_seed)

        sample_idx = list(range(len(labeltext)))
        numpy.random.shuffle(sample_idx)
        sample_idx = sample_idx[:args.subsample]

        np_idx = numpy.array(sample_idx)
        D = D[numpy.ix_(np_idx, np_idx)]
        labeltext = [ labeltext[idx] for idx in sample_idx ]

    ### do clustering
    Y = sch.linkage(D, method='single')
    Z1 = sch.dendrogram(Y, orientation='right', labels=labeltext)
    fig.savefig(dendrogram_out)
    notify('wrote dendrogram to: {}', dendrogram_out)

    ### make the dendrogram+matrix:
    fig = sourmash_fig.plot_composite_matrix(D, labeltext,
                                             show_labels=args.labels,
                                             show_indices=args.indices,
                                             vmin=args.vmin,
                                             vmax=args.vmax)
    fig.savefig(matrix_out)
    notify('wrote numpy distance matrix to: {}', matrix_out)

    # print out sample numbering for FYI.
    for i, name in enumerate(labeltext):
        print_results('{}\t{}', i, name)


def import_csv(args):
    "Import a CSV file full of signatures/hashes."
    p = argparse.ArgumentParser()
    p.add_argument('mash_csvfile')
    p.add_argument('-o', '--output', type=argparse.FileType('wt'),
                   default=sys.stdout, help='(default: stdout)')
    args = p.parse_args(args)

    with open(args.mash_csvfile, 'r') as fp:
        reader = csv.reader(fp)
        siglist = []
        for row in reader:
            hashfn = row[0]
            hashseed = int(row[1])

            # only support a limited import type, for now ;)
            assert hashfn == 'murmur64'
            assert hashseed == 42

            _, _, ksize, name, hashes = row
            ksize = int(ksize)

            hashes = hashes.strip()
            hashes = list(map(int, hashes.split(' ' )))

            e = sourmash_lib.MinHash(len(hashes), ksize)
            e.add_many(hashes)
            s = sig.SourmashSignature(e, filename=name)
            siglist.append(s)
            notify('loaded signature: {} {}', name, s.md5sum()[:8])

        notify('saving {} signatures to JSON', len(siglist))
        sig.save_signatures(siglist, args.output)


def dump(args):
    parser = argparse.ArgumentParser()
    parser.add_argument('filenames', nargs='+')
    parser.add_argument('-k', '--ksize', type=int, default=DEFAULT_LOAD_K, help='k-mer size (default: %(default)i)')
    args = parser.parse_args(args)

    for filename in args.filenames:
        notify('loading {}', filename)
        siglist = sig.load_signatures(filename, ksize=args.ksize)
        siglist = list(siglist)
        assert len(siglist) == 1

        s = siglist[0]

        fp = open(filename + '.dump.txt', 'w')
        fp.write(" ".join((map(str, s.minhash.get_hashes()))))
        fp.close()


def sbt_combine(args):
    parser = argparse.ArgumentParser()
    parser.add_argument('sbt_name', help='name to save SBT into')
    parser.add_argument('sbts', nargs='+',
                        help='SBTs to combine to a new SBT')
    parser.add_argument('-x', '--bf-size', type=float, default=1e5)

    sourmash_args.add_moltype_args(parser)

    args = parser.parse_args(args)
    moltype = sourmash_args.calculate_moltype(args)

    inp_files = list(args.sbts)
    notify('combining {} SBTs', len(inp_files))

    tree = sourmash_lib.load_sbt_index(inp_files.pop(0))

    for f in inp_files:
        new_tree = sourmash_lib.load_sbt_index(f)
        # TODO: check if parameters are the same for both trees!
        tree.combine(new_tree)

    notify('saving SBT under "{}".', args.sbt_name)
    tree.save(args.sbt_name)


def index(args):
    import sourmash_lib.sbtmh

    parser = argparse.ArgumentParser()
    parser.add_argument('sbt_name', help='name to save SBT into')
    parser.add_argument('signatures', nargs='+',
                        help='signatures to load into SBT')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='suppress non-error output')
    parser.add_argument('-k', '--ksize', type=int, default=None,
                        help='k-mer size for which to build the SBT.')
    parser.add_argument('-d', '--n_children', type=int, default=2,
                        help='Number of children for internal nodes')
    parser.add_argument('--traverse-directory', action='store_true',
                        help='load all signatures underneath this directory.')
    parser.add_argument('--append', action='store_true', default=False,
                        help='add signatures to an existing SBT.')
    parser.add_argument('-x', '--bf-size', type=float, default=1e5,
                        help='Bloom filter size used for internal nodes.')
    parser.add_argument('-s', '--sparseness', type=float, default=.0,
                        help='What percentage of internal nodes will not be saved. '
                             'Ranges from 0.0 (save all nodes) to 1.0 (no nodes saved)')

    sourmash_args.add_moltype_args(parser)

    args = parser.parse_args(args)
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    if args.append:
        tree = sourmash_lib.load_sbt_index(args.sbt_name)
    else:
        tree = sourmash_lib.create_sbt_index(args.bf_size,
                                             n_children=args.n_children)

    if args.traverse_directory:
        inp_files = list(sourmash_args.traverse_find_sigs(args.signatures))
    else:
        inp_files = list(args.signatures)

    if args.sparseness < 0 or args.sparseness > 1.0:
        error('sparseness must be in range [0.0, 1.0].')

    notify('loading {} files into SBT', len(inp_files))

    n = 0
    ksizes = set()
    moltypes = set()
    for f in inp_files:
        siglist = sig.load_signatures(f, ksize=args.ksize,
                                      select_moltype=moltype)

        # load all matching signatures in this file
        for ss in siglist:
            ksizes.add(ss.minhash.ksize)
            moltypes.add(sourmash_args.get_moltype(ss))

            leaf = sourmash_lib.sbtmh.SigLeaf(ss.md5sum(), ss)
            tree.add_node(leaf)
            n += 1

        # check to make sure we aren't loading incompatible signatures
        if len(ksizes) > 1 or len(moltypes) > 1:
            error('multiple k-mer sizes or molecule types present; fail.')
            error('specify --dna/--protein and --ksize as necessary')
            error('ksizes: {}; moltypes: {}',
                  ", ".join(map(str, ksizes)), ", ".join(moltypes))
            sys.exit(-1)

    # did we load any!?
    if n == 0:
        error('no signatures found to load into tree!? failing.')
        sys.exit(-1)

    notify('loaded {} sigs; saving SBT under "{}"', n, args.sbt_name)
    tree.save(args.sbt_name, sparseness=args.sparseness)


def search(args):
    from .search import search_databases

    parser = argparse.ArgumentParser()
    parser.add_argument('query', help='query signature')
    parser.add_argument('databases', help='signatures/SBTs to search',
                        nargs='+')
    parser.add_argument('--traverse-directory', action='store_true',
                        help='search all signatures underneath directories.')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='suppress non-error output')
    parser.add_argument('--threshold', default=0.08, type=float,
                        help='minimum threshold for reporting matches')
    parser.add_argument('--save-matches', type=argparse.FileType('wt'),
                        help='output matching signatures to this file.')
    parser.add_argument('--best-only', action='store_true',
                        help='report only the best match (with greater speed).')
    parser.add_argument('-n', '--num-results', default=3, type=int,
                        help='number of results to report')
    parser.add_argument('--containment', action='store_true',
                        help='evaluate containment rather than similarity')
    parser.add_argument('--scaled', type=float, default=0,
                        help='downsample query to this scaled factor (yields greater speed)')
    parser.add_argument('-o', '--output', type=argparse.FileType('wt'),
                        help='output CSV containing matches to this file')

    sourmash_args.add_ksize_arg(parser, DEFAULT_LOAD_K)
    sourmash_args.add_moltype_args(parser)

    args = parser.parse_args(args)
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    # set up the query.
    query = sourmash_args.load_query_signature(args.query,
                                               ksize=args.ksize,
                                               select_moltype=moltype)
    query_moltype = sourmash_args.get_moltype(query)
    query_ksize = query.minhash.ksize
    notify('loaded query: {}... (k={}, {})', query.name()[:30],
                                             query_ksize,
                                             query_moltype)

    # downsample if requested
    if args.scaled:
        if query.minhash.max_hash == 0:
            error('cannot downsample a signature not created with --scaled')
            sys.exit(-1)

        notify('downsampling query from scaled={} to {}',
               query.minhash.scaled, int(args.scaled))
        query.minhash = query.minhash.downsample_scaled(args.scaled)

    # set up the search databases
    databases = sourmash_args.load_sbts_and_sigs(args.databases,
                                                 query_ksize, query_moltype,
                                                 args.traverse_directory)

    if not len(databases):
        error('Nothing found to search!')
        sys.exit(-1)

    # do the actual search
    results = search_databases(query, databases,
                               args.threshold, args.containment,
                               args.best_only)

    n_matches = len(results)
    if args.best_only:
        args.num_results = 1

    if n_matches <= args.num_results:
        print_results('{} matches:'.format(len(results)))
    else:
        print_results('{} matches; showing first {}:',
               len(results), args.num_results)
        n_matches = args.num_results

    # output!
    print_results("similarity   match")
    print_results("----------   -----")
    for sr in results[:n_matches]:
        pct = '{:.1f}%'.format(sr.similarity*100)
        name = sr.match_sig._display_name(60)
        print_results('{:>6}       {}', pct, name)

    if args.best_only:
        notify("** reporting only one match because --best-only was set")

    if args.output:
        fieldnames = ['similarity', 'name', 'filename', 'md5']
        w = csv.DictWriter(args.output, fieldnames=fieldnames)

        w.writeheader()
        for sr in results:
            d = dict(sr._asdict())
            del d['match_sig']
            w.writerow(d)

    # save matching signatures upon request
    if args.save_matches:
        outname = args.save_matches.name
        notify('saving all matched signatures to "{}"', outname)
        sig.save_signatures([ sr.match_sig for sr in results ],
                            args.save_matches)


def categorize(args):
    parser = argparse.ArgumentParser()
    parser.add_argument('sbt_name', help='name of SBT to load')
    parser.add_argument('queries', nargs='+',
                        help='list of signatures to categorize')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='suppress non-error output')
    parser.add_argument('-k', '--ksize', type=int, default=None)
    parser.add_argument('--threshold', default=0.08, type=float)
    parser.add_argument('--traverse-directory', action="store_true")

    sourmash_args.add_moltype_args(parser)

    parser.add_argument('--csv', type=argparse.FileType('at'))
    parser.add_argument('--load-csv', default=None)

    args = parser.parse_args(args)
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    already_names = set()
    if args.load_csv:
        with open(args.load_csv, 'rt') as fp:
            r = csv.reader(fp)
            for row in r:
                already_names.add(row[0])

    tree = sourmash_lib.load_sbt_index(args.sbt_name)

    if args.traverse_directory:
        inp_files = set(sourmash_args.traverse_find_sigs(args.queries))
    else:
        inp_files = set(args.queries) - already_names

    inp_files = set(inp_files) - already_names

    notify('found {} files to query', len(inp_files))

    loader = sourmash_args.LoadSingleSignatures(inp_files,
                                                args.ksize, moltype)

    for queryfile, query, query_moltype, query_ksize in loader:
        notify('loaded query: {}... (k={}, {})', query.name()[:30],
               query_ksize, query_moltype)

        results = []
        search_fn = sourmash_lib.sbtmh.SearchMinHashesFindBest().search

        for leaf in tree.find(search_fn, query, args.threshold):
            if leaf.data.md5sum() != query.md5sum(): # ignore self.
                results.append((query.similarity(leaf.data), leaf.data))

        best_hit_sim = 0.0
        best_hit_query_name = ""
        if results:
            results.sort(key=lambda x: -x[0])   # reverse sort on similarity
            best_hit_sim, best_hit_query = results[0]
            notify('for {}, found: {:.2f} {}', query.name(),
                                               best_hit_sim,
                                               best_hit_query.name())
            best_hit_query_name = best_hit_query.name()
        else:
            notify('for {}, no match found', query.name())

        if args.csv:
            w = csv.writer(args.csv)
            w.writerow([queryfile, best_hit_query_name, best_hit_sim])

    if loader.skipped_ignore:
        notify('skipped/ignore: {}', loader.skipped_ignore)
    if loader.skipped_nosig:
        notify('skipped/nosig: {}', loader.skipped_nosig)


def gather(args):
    from .search import gather_databases

    parser = argparse.ArgumentParser()
    parser.add_argument('query', help='query signature')
    parser.add_argument('databases', help='signatures/SBTs to search',
                        nargs='+')
    parser.add_argument('--traverse-directory', action='store_true',
                        help='search all signatures underneath directories.')
    parser.add_argument('-o', '--output', type=argparse.FileType('wt'),
                        help='output CSV containing matches to this file')
    parser.add_argument('--save-matches', type=argparse.FileType('wt'),
                        help='save the matched signatures from the database to this file.')
    parser.add_argument('--threshold-bp', type=float, default=5e4,
                        help='threshold (in bp) for reporting results')
    parser.add_argument('--output-unassigned', type=argparse.FileType('wt'),
                        help='output unassigned portions of the query as a signature to this file')
    parser.add_argument('--scaled', type=float, default=0,
                        help='downsample query to this scaled factor')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='suppress non-error output')

    sourmash_args.add_ksize_arg(parser, DEFAULT_LOAD_K)
    sourmash_args.add_moltype_args(parser)

    args = parser.parse_args(args)
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    # load the query signature & figure out all the things
    query = sourmash_args.load_query_signature(args.query,
                                               ksize=args.ksize,
                                               select_moltype=moltype)
    query_moltype = sourmash_args.get_moltype(query)
    query_ksize = query.minhash.ksize
    notify('loaded query: {}... (k={}, {})', query.name()[:30],
                                             query_ksize,
                                             query_moltype)

    # verify signature was computed right.
    if query.minhash.max_hash == 0:
        error('query signature needs to be created with --scaled')
        sys.exit(-1)

    # downsample if requested
    if args.scaled:
        notify('downsampling query from scaled={} to {}',
               query.minhash.scaled, int(args.scaled))
        query.minhash = query.minhash.downsample_scaled(args.scaled)

    # empty?
    if not query.minhash.get_mins():
        error('no query hashes!? exiting.')
        sys.exit(-1)

    # set up the search databases
    databases = sourmash_args.load_sbts_and_sigs(args.databases,
                                                 query_ksize, query_moltype,
                                                 args.traverse_directory)

    if not len(databases):
        error('Nothing found to search!')
        sys.exit(-1)

    # pretty-printing code.
    def format_bp(bp):
        bp = float(bp)
        if bp < 500:
            return '{:.0f} bp '.format(bp)
        elif bp <= 500e3:
            return '{:.1f} kbp'.format(round(bp / 1e3, 1))
        elif bp < 500e6:
            return '{:.1f} Mbp'.format(round(bp / 1e6, 1))
        elif bp < 500e9:
            return '{:.1f} Gbp'.format(round(bp / 1e9, 1))
        return '???'

    found = []
    sum_found = 0
    for result, n_intersect_mins, new_max_hash, next_query in gather_databases(query, databases,
                                                     args.threshold_bp):
        # print interim result & save in a list for later use
        pct_query = '{:.1f}%'.format(result.f_orig_query*100)
        pct_genome = '{:.1f}%'.format(result.f_match*100)

        name = result.leaf._display_name(40)

        if not len(found):                # first result? print header.
            print_results("")
            print_results("overlap     p_query p_match ")
            print_results("---------   ------- --------")

        # print interim result & save in a list for later use
        pct_query = '{:.1f}%'.format(result.f_orig_query*100)
        pct_genome = '{:.1f}%'.format(result.f_match*100)

        name = result.leaf._display_name(40)

        print_results('{:9}   {:>6}  {:>6}      {}',
                      format_bp(result.intersect_bp), pct_query, pct_genome,
                      name)
        sum_found += n_intersect_mins
        found.append(result)


    # basic reporting
    print_results('\nfound {} matches total;', len(found))

    sum_found /= len(query.minhash.get_hashes())
    print_results('the recovered matches hit {:.1f}% of the query',
           sum_found * 100)
    print_results('')

    if not found:
        sys.exit(0)

    if args.output:
        fieldnames = ['intersect_bp', 'f_orig_query', 'f_match',
                      'f_unique_to_query', 'name', 'filename', 'md5']
        w = csv.DictWriter(args.output, fieldnames=fieldnames)
        w.writeheader()
        for result in found:
            d = dict(result._asdict())
            del d['leaf']                 # actual signature not in CSV.
            w.writerow(d)

    if args.save_matches:
        outname = args.save_matches.name
        notify('saving all matches to "{}"', outname)
        sig.save_signatures([ r.leaf for r in found ], args.save_matches)

    if args.output_unassigned:
        if not found:
            notify('nothing found - entire query signature unassigned.')
        if not query.minhash.get_mins():
            notify('no unassigned hashes! not saving.')
        else:
            outname = args.output_unassigned.name
            notify('saving unassigned hashes to "{}"', outname)

            e = sourmash_lib.MinHash(ksize=query_ksize, n=0,
                                     max_hash=new_max_hash)
            e.add_many(query.minhash.get_mins())
            sig.save_signatures([ sig.SourmashSignature(e) ],
                                args.output_unassigned)


def watch(args):
    "Build a signature from raw FASTA/FASTQ coming in on stdin, search."
    from sourmash_lib.sbtmh import search_minhashes, SearchMinHashesFindBest

    parser = argparse.ArgumentParser()
    parser.add_argument('sbt_name', help='name of SBT to search')
    parser.add_argument('inp_file', nargs='?', default='/dev/stdin')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='suppress non-error output')
    parser.add_argument('-o', '--output', type=argparse.FileType('wt'),
                        help='save signature generated from data here')
    parser.add_argument('--threshold', default=0.05, type=float,
                        help='minimum threshold for matches')
    parser.add_argument('--input-is-protein', action='store_true',
                        help='Consume protein sequences - no translation needed')
    sourmash_args.add_construct_moltype_args(parser)
    parser.add_argument('-n', '--num-hashes', type=int,
                        default=DEFAULT_N,
                        help='number of hashes to use in each sketch (default: %(default)i)')
    parser.add_argument('--name', type=str, default='stdin',
                        help='name to use for generated signature')
    sourmash_args.add_ksize_arg(parser, DEFAULT_LOAD_K)
    args = parser.parse_args(args)
    set_quiet(args.quiet)

    if args.input_is_protein and args.dna:
        notify('WARNING: input is protein, turning off DNA hashing.')
        args.dna = False
        args.protein = True

    if args.dna and args.protein:
        notify('ERROR: cannot use "watch" with both DNA and protein.')

    if args.dna:
        moltype = 'DNA'
        is_protein = False
    else:
        moltype = 'protein'
        is_protein = True

    tree = sourmash_lib.load_sbt_index(args.sbt_name)

    def get_ksize(tree):
        """Walk nodes in `tree` to find out ksize"""
        for node in tree.nodes.values():
            if isinstance(node, sourmash_lib.sbtmh.SigLeaf):
                return node.data.minhash.ksize

    # deduce ksize from the SBT we are loading
    ksize = args.ksize
    if ksize is None:
        ksize = get_ksize(tree)

    E = sourmash_lib.MinHash(ksize=ksize, n=args.num_hashes,
                             is_protein=is_protein)
    streamsig = sig.SourmashSignature(E, filename='stdin',
                                      name=args.name)

    notify('Computing signature for k={}, {} from stdin',
           ksize, moltype)

    def do_search():
        search_fn = SearchMinHashesFindBest().search

        results = []
        for leaf in tree.find(search_fn, streamsig, args.threshold):
            results.append((streamsig.similarity(leaf.data),
                            leaf.data))

        return results

    notify('reading sequences from stdin')
    screed_iter = screed.open(args.inp_file)
    watermark = WATERMARK_SIZE

    # iterate over input records
    n = 0
    for n, record in enumerate(screed_iter):
        # at each watermark, print status & check cardinality
        if n >= watermark:
            notify('\r... read {} sequences', n, end='')
            watermark += WATERMARK_SIZE

            if do_search():
                break

        if args.input_is_protein:
            E.add_protein(record.sequence)
        else:
            E.add_sequence(record.sequence, False)

    results = do_search()
    if not results:
        notify('... read {} sequences, no matches found.', n)
    else:
        results.sort(key=lambda x: -x[0])   # take best
        similarity, found_sig = results[0]
        print_results('FOUND: {}, at {:.3f}', found_sig.name(),
               similarity)

    if args.output:
        notify('saving signature to {}', args.output.name)
        sig.save_signatures([streamsig], args.output)


def storage(args):
    from .sbt import convert_cmd

    parser = argparse.ArgumentParser()
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='suppress non-error output')

    subparsers = parser.add_subparsers()
    convert_parser = subparsers.add_parser('convert')
    convert_parser.add_argument('sbt', help='SBT to convert')
    convert_parser.add_argument('-b', "--backend", type=str,
                                help='Backend to convert to')
    convert_parser.set_defaults(command='convert')

    args = parser.parse_args(args)
    set_quiet(args.quiet)
    if args.command == 'convert':
        convert_cmd(args.sbt, args.backend)
