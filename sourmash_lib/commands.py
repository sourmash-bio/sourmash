import argparse
import csv
import os
import os.path
import sys

import screed
import sourmash_lib
from . import signature as sig
from . import sourmash_args
from .logging import notify, error

DEFAULT_K = 31
DEFAULT_N = 500
WATERMARK_SIZE = 10000


def search(args):
    "Search a query sig against one or more signatures; report up to the 3 top matches."
    parser = argparse.ArgumentParser()
    parser.add_argument('query', help='query signature')
    parser.add_argument('against', nargs='+', help='list of signatures')
    parser.add_argument('--threshold', default=0.08, type=float)
    parser.add_argument('-n', '--num-results', default=3, type=int)
    parser.add_argument('-k', '--ksize', default=DEFAULT_K, type=int)
    parser.add_argument('-f', '--force', action='store_true')
    parser.add_argument('--save-matches', type=argparse.FileType('wt'))

    sourmash_args.add_moltype_args(parser)

    args = parser.parse_args(args)
    moltype = sourmash_args.calculate_moltype(args)

    # get the query signature
    query = sourmash_args.load_query_signature(args.query,
                                               select_ksize=args.ksize,
                                               select_moltype=moltype)
    query_moltype = sourmash_args.get_moltype(query)
    query_ksize = query.estimator.ksize
    notify('loaded query: {}... (k={}, {})', query.name()[:30],
                                             query_ksize,
                                             query_moltype)

    # get the signatures to query
    notify('loading db of signatures from {} files', len(args.against))
    against = []
    for filename in args.against:
        if filename == args.query and not args.force:
            notify('excluding query from database (file {})', filename)
            continue

        sl = sig.load_signatures(filename,
                                 select_ksize=args.ksize,
                                 select_moltype=moltype)

        for x in sl:
            against.append((x, filename))
    notify('loaded {} signatures total.', len(against))

    # compute query x db
    distances = []
    for (x, filename) in against:
        distance = query.similarity(x)
        if distance >= args.threshold:
            distances.append((distance, x, filename))

    # any matches? sort, show.
    if distances:
        distances.sort(reverse=True, key = lambda x: x[0])
        notify('{} matches; showing {}:', len(distances), args.num_results)
        for distance, match, filename in distances[:args.num_results]:

            print('\t', match.name(), '\t', "%.3f" % distance,
                  '\t', filename)

        if args.save_matches:
            outname = args.save_matches.name
            notify('saving all matches to "{}"', outname)
            sig.save_signatures([ m for (d, m, f) in distances ],
                                args.save_matches)
    else:
        notify('** no matches in {} signatures', len(against))


def compute(args):
    """Compute the signature for one or more files.

    Use cases:
        sourmash compute multiseq.fa              => multiseq.fa.sig, etc.
        sourmash compute genome.fa --singleton    => genome.fa.sig
        sourmash compute file1.fa file2.fa --name => specify w/-o
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('filenames', nargs='+',
                        help='file(s) of sequences')

    sourmash_args.add_moltype_args(parser, default_dna=True)

    parser.add_argument('--input-is-protein', action='store_true')
    parser.add_argument('-k', '--ksizes',
                        default=str(DEFAULT_K),
                        help='comma-separated list of k-mer sizes (default: %(default)s)')
    parser.add_argument('-n', '--num-hashes', type=int,
                        default=DEFAULT_N,
                        help='number of hashes to use in each sketch (default: %(default)i)')
    parser.add_argument('--check-sequence', action='store_true')
    parser.add_argument('-f', '--force', action='store_true')
    parser.add_argument('-o', '--output', type=argparse.FileType('wt'))
    parser.add_argument('--email', type=str, default='')
    parser.add_argument('--singleton', action='store_true')
    parser.add_argument('--name', type=str, default='')
    parser.add_argument('--name-from-first', action='store_true')
    parser.add_argument('--with-cardinality', action='store_true',
                        help='calculate # of k-mers in input sequences')
    parser.add_argument('--track-abundance', action='store_true',
                        help='track k-mer abundances')
    parser.add_argument('--scaled', type=float,
                        help='set resolution of signature size')
    parser.add_argument('--seed', type=int, help='hash seed',
                        default=sourmash_lib.DEFAULT_SEED)
    args = parser.parse_args(args)

    if args.input_is_protein and args.dna:
        notify('WARNING: input is protein, turning off DNA hashing')
        args.dna = False
        args.protein = True

    if not args.dna and args.protein:
        if args.with_cardinality:
            error('Cannot compute cardinality for protein sequences.')
            sys.exit(-1)

    if args.scaled:
        if args.num_hashes != 0:
            notify('setting num_hashes to 0 because --scaled is set')
            args.num_hashes = 0

    notify('computing signatures for files: {}', ", ".join(args.filenames))

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

    if args.name and not args.output:
        error("must specify -o with --name")
        sys.exit(-1)

    def make_estimators():
        seed = args.seed
        max_hash = 0
        if args.scaled:
            max_hash = 2**64 / float(args.scaled)

        # one estimator for each ksize
        Elist = []
        for k in ksizes:
            if args.protein:
                E = sourmash_lib.Estimators(ksize=k, n=args.num_hashes,
                                            is_protein=True,
                                    track_abundance=args.track_abundance,
                                            max_hash=max_hash,
                                            seed=seed)
                Elist.append(E)
            if args.dna:
                E = sourmash_lib.Estimators(ksize=k, n=args.num_hashes,
                                            is_protein=False,
                                    with_cardinality=args.with_cardinality,
                                    track_abundance=args.track_abundance,
                                            max_hash=max_hash,
                                            seed=seed)
                Elist.append(E)
        return Elist

    def add_seq(Elist, seq, input_is_protein, check_sequence):
        for E in Elist:
            if input_is_protein:
                E.mh.add_protein(seq)
            else:
                E.add_sequence(seq, not check_sequence)

    def build_siglist(email, Elist, filename, name=None):
        return [ sig.SourmashSignature(email, E, filename=filename,
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

    notify('Computing signature for ksizes: {}', str(ksizes))

    if args.with_cardinality:
        print('Calculating k-mer cardinality of input sequences.',
              file=sys.stderr)

    if args.track_abundance:
        print('Tracking abundance of input k-mers.',
              file=sys.stderr)

    if not args.name:
        for filename in args.filenames:
            sigfile = os.path.basename(filename) + '.sig'
            if not args.output and os.path.exists(sigfile) and not \
                args.force:
                print('skipping', filename, '- already done',
                      file=sys.stderr)
                continue

            if args.singleton:
                siglist = []
                for n, record in enumerate(screed.open(filename)):
                    # make estimators for each sequence
                    Elist = make_estimators()
                    add_seq(Elist, record.sequence,
                            args.input_is_protein, args.check_sequence)

                    siglist += build_siglist(args.email, Elist, filename,
                                             name=record.name)
                print('calculated {} signatures for {} sequences in {}'.\
                          format(len(siglist), n + 1, filename))
            else:
                # make estimators for the whole file
                Elist = make_estimators()

                # consume & calculate signatures
                print('... reading sequences from', filename,
                      file=sys.stderr)
                name = None
                for n, record in enumerate(screed.open(filename)):
                    if n % 10000 == 0:
                        if n:
                            print('...', filename, n, file=sys.stderr)
                        elif args.name_from_first:
                            name = record.name

                    s = record.sequence
                    add_seq(Elist, record.sequence,
                            args.input_is_protein, args.check_sequence)

                siglist = build_siglist(args.email, Elist, filename, name)
                print('calculated {} signatures for {} sequences in {}'.\
                          format(len(siglist), n + 1, filename))
            # at end, save!
            save_siglist(siglist, args.output, sigfile)
    else:                             # single name specified - combine all
        # make estimators for the whole file
        Elist = make_estimators()

        for filename in args.filenames:
            # consume & calculate signatures
            print('... reading sequences from', filename,
                  file=sys.stderr)
            for n, record in enumerate(screed.open(filename)):
                if n % 10000 == 0 and n:
                    print('...', filename, n, file=sys.stderr)

                add_seq(Elist, record.sequence,
                        args.input_is_protein, args.check_sequence)

        siglist = build_siglist(args.email, Elist, filename,
                                name=args.name)
        print('calculated {} signatures for {} sequences taken from {}'.\
               format(len(siglist), n + 1, " ".join(args.filenames)))
        # at end, save!
        save_siglist(siglist, args.output)


def compare(args):
    "Compare multiple signature files and create a distance matrix."
    import numpy

    parser = argparse.ArgumentParser()
    parser.add_argument('signatures', nargs='+', help='list of signatures')
    parser.add_argument('-k', '--ksize', type=int, default=DEFAULT_K, help='k-mer size (default: %(default)s)')
    parser.add_argument('-o', '--output')
    parser.add_argument('--ignore-abundance', action='store_true',
                        help='do NOT use k-mer abundances if present')
    args = parser.parse_args(args)

    # load in the various signatures
    siglist = []
    for filename in args.signatures:
        print('loading', filename, file=sys.stderr)
        loaded = sig.load_signatures(filename, select_ksize=args.ksize)
        loaded = list(loaded)
        if not loaded:
            notify('warning: no signatures loaded at given ksize from {}',
                   filename)
        siglist.extend(loaded)

    if len(siglist) == 0:
        error('no signatures!')
        sys.exit(-1)

    # build the distance matrix
    D = numpy.zeros([len(siglist), len(siglist)])
    numpy.set_printoptions(precision=3, suppress=True)

    # do all-by-all calculation
    labeltext = []
    for i, E in enumerate(siglist):
        for j, E2 in enumerate(siglist):
            D[i][j] = E.similarity(E2, args.ignore_abundance)

        print('%d-%20s\t%s' % (i, E.name(), D[i, :, ],))
        labeltext.append(E.name())

    notify('min similarity in matrix: {}', numpy.min(D))

    # shall we output a matrix?
    if args.output:
        labeloutname = args.output + '.labels.txt'
        notify('saving labels to: {}', labeloutname)
        with open(labeloutname, 'w') as fp:
            fp.write("\n".join(labeltext))

        notify('saving distance matrix to: {}', args.output)
        with open(args.output, 'wb') as fp:
            numpy.save(fp, D)


def plot(args):
    "Produce a clustering and plot."
    import numpy
    import scipy
    import pylab
    import scipy.cluster.hierarchy as sch
    from . import fig as sourmash_fig

    # set up cmd line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('distances', help="output from 'sourmash compare'")
    parser.add_argument('--pdf', action='store_true')
    parser.add_argument('--labels', action='store_true')
    parser.add_argument('--indices', action='store_false')
    parser.add_argument('--vmax', default=1.0, type=float,
                        help='(default: %(default)f)')
    parser.add_argument('--vmin', default=0.0, type=float,
                        help='(default: %(default)f)')
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

    Y = sch.linkage(D, method='single') # cluster!
    Z1 = sch.dendrogram(Y, orientation='right', labels=labeltext)
    fig.savefig(dendrogram_out)
    notify('wrote {}', dendrogram_out)

    ### make the dendrogram+matrix:
    fig = sourmash_fig.plot_composite_matrix(D, labeltext,
                                             show_labels=args.labels,
                                             show_indices=args.indices,
                                             vmin=args.vmin,
                                             vmax=args.vmax)
    fig.savefig(matrix_out)
    notify('wrote {}', matrix_out)

    # print out sample numbering for FYI.
    for i, name in enumerate(labeltext):
        print(i, '\t', name)


def import_csv(args):
    "Import a CSV file full of signatures/hashes."
    p = argparse.ArgumentParser()
    p.add_argument('mash_csvfile')
    p.add_argument('-o', '--output', type=argparse.FileType('wt'),
                   default=sys.stdout, help='(default: stdout)')
    p.add_argument('--email', type=str, default='', help='(default: %(default)s)')
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

            e = sourmash_lib.Estimators(len(hashes), ksize)
            e.add_many(hashes)
            s = sig.SourmashSignature(args.email, e, filename=name)
            siglist.append(s)
            notify('loaded signature: {} {}', name, s.md5sum()[:8])

        notify('saving {} signatures to JSON', len(siglist))
        sig.save_signatures(siglist, args.output)


def dump(args):
    parser = argparse.ArgumentParser()
    parser.add_argument('filenames', nargs='+')
    parser.add_argument('-k', '--ksize', type=int, default=DEFAULT_K, help='k-mer size (default: %(default)i)')
    args = parser.parse_args(args)

    for filename in args.filenames:
        notify('loading {}', filename)
        siglist = sig.load_signatures(filename, select_ksize=args.ksize)
        siglist = list(siglist)
        assert len(siglist) == 1

        s = siglist[0]

        fp = open(filename + '.dump.txt', 'w')
        fp.write(" ".join((map(str, s.estimator.get_hashes()))))
        fp.close()


def sbt_index(args):
    from sourmash_lib.sbt import SBT, GraphFactory
    from sourmash_lib.sbtmh import search_minhashes, SigLeaf

    parser = argparse.ArgumentParser()
    parser.add_argument('sbt_name', help='name to save SBT into')
    parser.add_argument('signatures', nargs='+',
                        help='signatures to load into SBT')
    parser.add_argument('-k', '--ksize', type=int, default=None)
    parser.add_argument('--traverse-directory', action='store_true')
    parser.add_argument('-x', '--bf-size', type=float, default=1e5)

    sourmash_args.add_moltype_args(parser)

    args = parser.parse_args(args)
    moltype = sourmash_args.calculate_moltype(args)

    factory = GraphFactory(1, args.bf_size, 4)
    tree = SBT(factory)

    if args.traverse_directory:
        inp_files = list(sourmash_args.traverse_find_sigs(args.signatures))
    else:
        inp_files = list(args.signatures)


    notify('loading {} files into SBT', len(inp_files))

    n = 0
    ksizes = set()
    moltypes = set()
    for f in inp_files:
        siglist = sig.load_signatures(f, select_ksize=args.ksize,
                                      select_moltype=moltype)

        # load all matching signatures in this file
        for ss in siglist:
            ksizes.add(ss.estimator.ksize)
            moltypes.add(sourmash_args.get_moltype(ss))

            leaf = SigLeaf(ss.md5sum(), ss)
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
    tree.save(args.sbt_name)


def sbt_search(args):
    from sourmash_lib.sbt import SBT, GraphFactory
    from sourmash_lib.sbtmh import search_minhashes, SigLeaf
    from sourmash_lib.sbtmh import SearchMinHashesFindBest

    parser = argparse.ArgumentParser()
    parser.add_argument('sbt_name', help='name of SBT to load')
    parser.add_argument('query', help='signature to query')
    parser.add_argument('-k', '--ksize', type=int, default=DEFAULT_K)
    parser.add_argument('--threshold', default=0.08, type=float)
    parser.add_argument('--save-matches', type=argparse.FileType('wt'))
    parser.add_argument('--best-only', action='store_true')

    sourmash_args.add_moltype_args(parser)
    args = parser.parse_args(args)
    moltype = sourmash_args.calculate_moltype(args)

    search_fn = search_minhashes
    if args.best_only:
        search_fn = SearchMinHashesFindBest().search

    tree = SBT.load(args.sbt_name, leaf_loader=SigLeaf.load)
    query = sourmash_args.load_query_signature(args.query,
                                               select_ksize=args.ksize,
                                               select_moltype=moltype)
    query_moltype = sourmash_args.get_moltype(query)
    query_ksize = query.estimator.ksize
    notify('loaded query: {}... (k={}, {})', query.name()[:30],
                                             query_ksize,
                                             query_moltype)

    results = []
    for leaf in tree.find(search_fn, query, args.threshold):
        results.append((query.similarity(leaf.data), leaf.data))
        #results.append((leaf.data.similarity(ss), leaf.data))

    results.sort(key=lambda x: -x[0])   # reverse sort on similarity
    for (similarity, query) in results:
        print('{:.2f} {}'.format(similarity, query.name()))

    if args.save_matches:
        outname = args.save_matches.name
        notify('saving all matches to "{}"', outname)
        sig.save_signatures([ m for (sim, m) in results ],
                            args.save_matches)


def categorize(args):
    from sourmash_lib.sbt import SBT, GraphFactory
    from sourmash_lib.sbtmh import search_minhashes, SigLeaf
    from sourmash_lib.sbtmh import SearchMinHashesFindBest

    parser = argparse.ArgumentParser()
    parser.add_argument('sbt_name', help='name of SBT to load')
    parser.add_argument('queries', nargs='+',
                        help='list of signatures to categorize')
    parser.add_argument('-k', '--ksize', type=int, default=None)
    parser.add_argument('--threshold', default=0.08, type=float)
    parser.add_argument('--traverse-directory', action="store_true")

    sourmash_args.add_moltype_args(parser)

    parser.add_argument('--csv', type=argparse.FileType('at'))
    parser.add_argument('--load-csv', default=None)

    args = parser.parse_args(args)
    moltype = sourmash_args.calculate_moltype(args)

    already_names = set()
    if args.load_csv:
        with open(args.load_csv, 'rt') as fp:
            r = csv.reader(fp)
            for row in r:
                already_names.add(row[0])

    tree = SBT.load(args.sbt_name, leaf_loader=SigLeaf.load)

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
        search_fn = SearchMinHashesFindBest().search

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

def sbt_gather(args):
    from sourmash_lib.sbt import SBT, GraphFactory
    from sourmash_lib.sbtmh import search_minhashes, SigLeaf
    from sourmash_lib.sbtmh import SearchMinHashesFindBestIgnoreMaxHash

    parser = argparse.ArgumentParser()
    parser.add_argument('sbt_name', help='name of SBT to search')
    parser.add_argument('query', help='query signature')
    parser.add_argument('-k', '--ksize', type=int, default=DEFAULT_K)
    parser.add_argument('--threshold', default=0.05, type=float)
    parser.add_argument('-o', '--output', type=argparse.FileType('wt'))
    parser.add_argument('--csv', type=argparse.FileType('wt'))
    parser.add_argument('--save-matches', type=argparse.FileType('wt'))

    sourmash_args.add_moltype_args(parser)

    args = parser.parse_args(args)
    moltype = sourmash_args.calculate_moltype(args)

    tree = SBT.load(args.sbt_name, leaf_loader=SigLeaf.load)
    query = sourmash_args.load_query_signature(args.query,
                                               select_ksize=args.ksize,
                                               select_moltype=moltype)
    query_moltype = sourmash_args.get_moltype(query)
    query_ksize = query.estimator.ksize
    notify('loaded query: {}... (k={}, {})', query.name()[:30],
                                             query_ksize,
                                             query_moltype)

    if query.estimator.max_hash == 0:
        error('query signature needs to be created with --scaled')
        error('or using --with-cardinality.')
        sys.exit(-1)

    notify('query signature has max_hash: {}', query.estimator.max_hash)
    orig_query = query

    R_metagenome = 2**64 / float(orig_query.estimator.max_hash)

    new_mins = query.estimator.get_hashes()
    e = sourmash_lib.Estimators(ksize=args.ksize, n=len(new_mins))
    e.update(query.estimator)
    query = sig.SourmashSignature('', e)

    sum_found = 0.
    found = []
    while 1:
        search_fn = SearchMinHashesFindBestIgnoreMaxHash().search

        results = []
        # use super low threshold for this part of the search
        for leaf in tree.find(search_fn, query, 0.00001):
            results.append((query.estimator.similarity_ignore_maxhash(leaf.data.estimator), leaf.data))

        if not len(results):          # no matches at all!
            break

        # take the best result
        results.sort(key=lambda x: -x[0])   # reverse sort on similarity
        best_sim, best_ss = results[0]

        # subtract found hashes from search hashes, construct new search
        new_mins = set(query.estimator.get_hashes())
        found_mins = best_ss.estimator.get_hashes()

        if best_ss.estimator.max_hash:
            R_genome = 2**64 / float(best_ss.estimator.max_hash)
        elif best_ss.estimator.hll:
            genome_size = best_ss.estimator.hll.estimate_cardinality()
            genome_max_hash = max(found_mins)
            R_genome = float(genome_size) / float(genome_max_hash)
        else:
            error('Best hash match in sbt_gather has no cardinality')
            error('Please prepare database of sequences with --scaled')
            error('...or with --with-cardinality')
            sys.exit(-1)

        R_comparison = max(R_metagenome, R_genome)
        new_max_hash = 2**64 / float(R_comparison)
        new_mins = set([ i for i in new_mins if i < new_max_hash ])
        found_mins = set([ i for i in found_mins if i < new_max_hash ])

        # intersection:
        intersect_mins = new_mins.intersection(found_mins)

        if len(intersect_mins) < 5:   # hard cutoff for now
            notify('found only {} hashes in common.', len(intersect_mins))
            notify('this is below a sane threshold => exiting.')
            break

        # first denominator - genome size
        genome_n_mins = len(found_mins)
        f_genome = len(intersect_mins) / float(genome_n_mins)

        # second denominator - metagenome size
        query_n_mins = len(orig_query.estimator.get_hashes())
        f_query = len(intersect_mins) / float(query_n_mins)

        # print interim & save
        notify('found: {:.2f} {:.2f} {}', f_genome, f_query,
               best_ss.name())
        found.append((f_genome, best_ss))

        new_mins -= set(found_mins)
        e = sourmash_lib.Estimators(ksize=args.ksize, n=len(new_mins))
        e.add_many(new_mins)
        query = sig.SourmashSignature('', e)

    notify('found {}, total fraction {:.3f}', len(found), sum_found)
    notify('')

    if not found:
        sys.exit(0)

    found.sort(key=lambda x: x[0])
    found.reverse()

    notify('Composition:')
    for (frac, leaf_sketch) in found:
        notify('{:.2f} {}', frac, leaf_sketch.name())

    if args.output:
        print('Composition:', file=args.output)
        for (frac, leaf_sketch) in found:
            print('{:.2f} {}'.format(frac, leaf_sketch.name()),
                  file=args.output)

    if args.csv:
        fieldnames = ['fraction', 'name', 'sketch_kmers']
        w = csv.DictWriter(args.csv, fieldnames=fieldnames)

        w.writeheader()
        for (frac, leaf_sketch) in found:
            cardinality = leaf_sketch.estimator.hll.estimate_cardinality()
            w.writerow(dict(fraction=frac, name=leaf_sketch.name(),
                            sketch_kmers=cardinality))
    if args.save_matches:
        outname = args.save_matches.name
        notify('saving all matches to "{}"', outname)
        sig.save_signatures([ ss for (f, ss) in found ],
                            args.save_matches)


def watch(args):
    "Build a signature from raw FASTA/FASTQ coming in on stdin, search."
    from sourmash_lib.sbt import SBT, GraphFactory
    from sourmash_lib.sbtmh import search_minhashes, SigLeaf
    from sourmash_lib.sbtmh import SearchMinHashesFindBest

    parser = argparse.ArgumentParser()
    parser.add_argument('sbt_name', help='name of SBT to search')
    parser.add_argument('inp_file', nargs='?', default='/dev/stdin')
    parser.add_argument('-o', '--output', type=argparse.FileType('wt'))
    parser.add_argument('-k', '--ksize', type=int, default=DEFAULT_K)
    parser.add_argument('--threshold', default=0.05, type=float)
    parser.add_argument('--input-is-protein', action='store_true')
    sourmash_args.add_moltype_args(parser, default_dna=True)
    parser.add_argument('-n', '--num-hashes', type=int,
                        default=DEFAULT_N,
                        help='number of hashes to use in each sketch (default: %(default)i)')
    parser.add_argument('--name', type=str, default='stdin')
    args = parser.parse_args(args)

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

    E = sourmash_lib.Estimators(ksize=args.ksize, n=args.num_hashes,
                                is_protein=is_protein)
    streamsig = sig.SourmashSignature('', E, filename='stdin',
                                      name=args.name)

    notify('Computing signature for k={}, {} from stdin',
           args.ksize, moltype)


    tree = SBT.load(args.sbt_name, leaf_loader=SigLeaf.load)

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
            notify('... read {} sequences', n)
            watermark += WATERMARK_SIZE

            if do_search():
                break

        if args.input_is_protein:
            E.mh.add_protein(record.sequence)
        else:
            E.add_sequence(record.sequence, False)

    results = do_search()
    if not results:
        notify('... read {} sequences, no matches found.', n)
    else:
        results.sort(key=lambda x: -x[0])   # take best
        similarity, found_sig = results[0]
        notify('FOUND: {}, at {:.3f}', found_sig.name(),
               similarity)

    if args.output:
        sig.save_signatures([streamsig], args.output)


def convert(args):
    "Convert YAML to JSON"

    import sourmash_lib.signature

    parser = argparse.ArgumentParser("""
Ensure that signature files in YAML (old format) are converted to signature files in JSON (new format).
The JSON-for-sure signature files are created in the same directory as the YAML-may-be files and are added
the extension ".json".
""")
    parser.add_argument('--minified', action="store_true",
                        help='Store the JSON minified (uses less space but is less human-readable) or not.')
    parser.add_argument('-f', '--force',
                        action="store_true",
                        help="Overwrite existing files.")
    parser.add_argument('path', nargs='*',
                        help='Path to YAML file with signatures, or directory containing signature files.')

    args = parser.parse_args(args)

    if args.minified:
        kwargs = {'indent': None,
                  'sort_keys': False}
    else:
        kwargs = {}


    filenames = list()

    while len(args.path) > 0:
        path = args.path.pop()

        if not os.path.exists(path):
            notify("The path name %s does not exist" % path)
            sys.exit(1)

        if os.path.isdir(path):
            notify("The path %s is a directory (and we will search signatures in it)." % path)
            args.path.extend(os.path.join(path, x) for x in os.listdir(path))
            continue

        # "path" is a file (not a directory) past this point
        if not path.endswith(".sig"):
            notify("The file name %s does not end with '.sig'. Skipping." % path)
            continue

        # fail early if output already existing
        out_fn = path + ".json"
        if not args.force and os.path.exists(out_fn):
            notify("The output file %s is already present. Use --force to force overwriting." % out_fn)
            sys.exit(1)

        # path is a file and should be converted
        filenames.append(path)


    for i, path in enumerate(filenames, 1):
        notify("\rConverting file %i/%i" % (i, len(filenames)), end="", flush=True)
        with open(path) as fh:
            signatures = tuple(sourmash_lib.signature.load_signatures(fh))

        out_fn = path + ".json"
        with open(out_fn, 'w') as fh:
            sourmash_lib.signature.save_signatures(signatures, fp=fh, **kwargs)
    notify("\rConverting file %i/%i" % (i, len(filenames)))
