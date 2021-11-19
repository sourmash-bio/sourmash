"""
Functions implementing the main command-line subcommands.
"""
import csv
import os
import os.path
import sys

import screed
from .compare import (compare_all_pairs, compare_serial_containment,
                      compare_serial_max_containment)
from . import MinHash
from .sbtmh import load_sbt_index, create_sbt_index
from . import signature as sig
from . import sourmash_args
from .logging import notify, error, print_results, set_quiet
from .sourmash_args import (FileOutput, FileOutputCSV,
                            SaveSignaturesToLocation)
from .search import prefetch_database
from .index import LazyLinearIndex

WATERMARK_SIZE = 10000


def compare(args):
    "Compare multiple signature files and create a distance matrix."
    import numpy

    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)
    picklist = sourmash_args.load_picklist(args)

    inp_files = list(args.signatures)
    if args.from_file:
        more_files = sourmash_args.load_pathlist_from_file(args.from_file)
        inp_files.extend(more_files)

    progress = sourmash_args.SignatureLoadingProgress()

    # load in the various signatures
    siglist = []
    ksizes = set()
    moltypes = set()
    for filename in inp_files:
        notify(f"loading '{filename}'", end='\r')
        loaded = sourmash_args.load_file_as_signatures(filename,
                                                       ksize=args.ksize,
                                                       select_moltype=moltype,
                                                       picklist=picklist,
                                                       yield_all_files=args.force,
                                                       progress=progress)
        loaded = list(loaded)
        if not loaded:
            notify(f'\nwarning: no signatures loaded at given ksize/molecule type/picklist from {filename}')
        siglist.extend(loaded)

        # track ksizes/moltypes
        for s in loaded:
            ksizes.add(s.minhash.ksize)
            moltypes.add(sourmash_args.get_moltype(s))

        # error out while loading if we have more than one ksize/moltype
        if len(ksizes) > 1 or len(moltypes) > 1:
            break

    if not siglist:
        error('no signatures found! exiting.')
        sys.exit(-1)

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
    notify(f'loaded {format(len(siglist))} signatures total.')

    if picklist:
        sourmash_args.report_picklist(args, picklist)

    # check to make sure they're potentially compatible - either using
    # scaled, or not.
    scaled_sigs = [s.minhash.scaled for s in siglist]
    is_scaled = all(scaled_sigs)
    is_scaled_2 = any(scaled_sigs)

    # complain if it's not all one or the other
    if is_scaled != is_scaled_2:
        error('cannot mix scaled signatures with bounded signatures')
        sys.exit(-1)

    is_containment = False
    if args.containment or args.max_containment:
        is_containment = True

        if args.containment and args.max_containment:
            notify("ERROR: cannot specify both --containment and --max-containment!")
            sys.exit(-1)

    # complain if --containment and not is_scaled
    if is_containment and not is_scaled:
        error('must use scaled signatures with --containment and --max-containment')
        sys.exit(-1)

    # notify about implicit --ignore-abundance:
    if is_containment:
        track_abundances = any(( s.minhash.track_abundance for s in siglist ))
        if track_abundances:
            notify('NOTE: --containment and --max-containment ignore signature abundances.')

    # if using --scaled, downsample appropriately
    printed_scaled_msg = False
    if is_scaled:
        max_scaled = max(s.minhash.scaled for s in siglist)
        for s in siglist:
            if s.minhash.scaled != max_scaled:
                if not printed_scaled_msg:
                    notify(f'downsampling to scaled value of {format(max_scaled)}')
                    printed_scaled_msg = True
                s.minhash = s.minhash.downsample(scaled=max_scaled)

    if len(siglist) == 0:
        error('no signatures!')
        sys.exit(-1)

    notify('')

    # build the distance matrix
    numpy.set_printoptions(precision=3, suppress=True)

    # do all-by-all calculation

    labeltext = [str(item) for item in siglist]
    if args.containment:
        similarity = compare_serial_containment(siglist)
    elif args.max_containment:
        similarity = compare_serial_max_containment(siglist)
    else:
        similarity = compare_all_pairs(siglist, args.ignore_abundance,
                                       n_jobs=args.processes)

    if len(siglist) < 30:
        for i, E in enumerate(siglist):
            # for small matrices, pretty-print some output
            name_num = '{}-{}'.format(i, str(E))
            if len(name_num) > 20:
                name_num = name_num[:17] + '...'
            print_results('{:20s}\t{}'.format(name_num, similarity[i, :, ],))

    print_results('min similarity in matrix: {:.3f}', numpy.min(similarity))
    # shall we output a matrix?
    if args.output:
        labeloutname = args.output + '.labels.txt'
        notify(f'saving labels to: {labeloutname}')
        with open(labeloutname, 'w') as fp:
            fp.write("\n".join(labeltext))

        notify(f'saving comparison matrix to: {args.output}')
        with open(args.output, 'wb') as fp:
            numpy.save(fp, similarity)

    # output CSV?
    if args.csv:
        with FileOutputCSV(args.csv) as csv_fp:
            w = csv.writer(csv_fp)
            w.writerow(labeltext)

            for i in range(len(labeltext)):
                y = []
                for j in range(len(labeltext)):
                    y.append('{}'.format(similarity[i][j]))
                w.writerow(y)


def plot(args):
    "Produce a clustering matrix and plot."
    import matplotlib as mpl
    mpl.use('Agg')
    import numpy
    import pylab
    import scipy.cluster.hierarchy as sch
    from . import fig as sourmash_fig

    # load files
    D_filename = args.distances
    labelfilename = D_filename + '.labels.txt'

    notify(f'loading comparison matrix from {D_filename}...')
    D = numpy.load(open(D_filename, 'rb'))
    # not sure how to change this to use f-strings
    notify('...got {} x {} matrix.', *D.shape)

    if args.labeltext:
        labelfilename = args.labeltext
    notify(f'loading labels from {labelfilename}')
    labeltext = [ x.strip() for x in open(labelfilename) ]
    if len(labeltext) != D.shape[0]:
        error('{} labels != matrix size, exiting', len(labeltext))
        sys.exit(-1)

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

    hist_out = os.path.basename(D_filename) + '.hist'
    if args.pdf:
        hist_out += '.pdf'
    else:
        hist_out += '.png'

    # output to a different directory?
    if args.output_dir:
        if not os.path.isdir(args.output_dir):
            os.mkdir(args.output_dir)
        dendrogram_out = os.path.join(args.output_dir, dendrogram_out)
        matrix_out = os.path.join(args.output_dir, matrix_out)
        hist_out = os.path.join(args.output_dir, hist_out)

    # make the histogram
    notify(f'saving histogram of matrix values => {hist_out}')
    fig = pylab.figure(figsize=(8,5))
    pylab.hist(numpy.array(D.flat), bins=100)
    fig.savefig(hist_out)

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
    sch.dendrogram(Y, orientation='right', labels=labeltext)
    fig.savefig(dendrogram_out)
    notify(f'wrote dendrogram to: {dendrogram_out}')

    ### make the dendrogram+matrix:
    (fig, rlabels, rmat) = sourmash_fig.plot_composite_matrix(D, labeltext,
                                             show_labels=args.labels,
                                             show_indices=args.indices,
                                             vmin=args.vmin,
                                             vmax=args.vmax,
                                             force=args.force)
    fig.savefig(matrix_out)
    notify(f'wrote numpy distance matrix to: {matrix_out}')

    if len(labeltext) < 30:
        # for small matrices, print out sample numbering for FYI.
        for i, name in enumerate(labeltext):
            print_results('{}\t{}', i, name)

    # write out re-ordered matrix and labels
    if args.csv:
        with FileOutputCSV(args.csv) as csv_fp:
            w = csv.writer(csv_fp)
            w.writerow(rlabels)

            for i in range(len(rlabels)):
                y = []
                for j in range(len(rlabels)):
                    y.append('{}'.format(rmat[i][j]))
                w.writerow(y)
        notify(f'Wrote clustered matrix and labels out to {args.csv}')


def import_csv(args):
    "Import a CSV file full of signatures/hashes."

    with open(args.mash_csvfile, newline='') as fp:
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

            e = MinHash(len(hashes), ksize)
            e.add_many(hashes)
            s = sig.SourmashSignature(e, filename=name)
            siglist.append(s)
            notify(f'loaded signature: {name} {s.md5sum()[:8]}')

        notify(f'saving {len(siglist)} signatures to JSON')
        with FileOutput(args.output, 'wt') as outfp:
            sig.save_signatures(siglist, outfp)


def sbt_combine(args):
    inp_files = list(args.sbts)
    notify(f'combining {len(inp_files)} SBTs')

    tree = load_sbt_index(inp_files.pop(0))

    for f in inp_files:
        new_tree = load_sbt_index(f)
        # TODO: check if parameters are the same for both trees!
        tree.combine(new_tree)

    notify(f'saving SBT under "{args.sbt_name}".')
    tree.save(args.sbt_name)


def index(args):
    """
    Build a Sequence Bloom Tree index of the given signatures.
    """
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)
    picklist = sourmash_args.load_picklist(args)

    if args.append:
        tree = load_sbt_index(args.sbt_name)
    else:
        tree = create_sbt_index(args.bf_size, n_children=args.n_children)

    if args.sparseness < 0 or args.sparseness > 1.0:
        error('sparseness must be in range [0.0, 1.0].')

    if args.scaled:
        args.scaled = int(args.scaled)
        notify(f'downsampling signatures to scaled={args.scaled}')

    inp_files = list(args.signatures)
    if args.from_file:
        more_files = sourmash_args.load_pathlist_from_file(args.from_file)
        inp_files.extend(more_files)

    if not inp_files:
        error("ERROR: no files to index!? Supply on command line or use --from-file")
        sys.exit(-1)

    notify(f'loading {len(inp_files)} files into SBT')

    progress = sourmash_args.SignatureLoadingProgress()

    n = 0
    ksizes = set()
    moltypes = set()
    nums = set()
    scaleds = set()
    for f in inp_files:
        siglist = sourmash_args.load_file_as_signatures(f,
                                                        ksize=args.ksize,
                                                        select_moltype=moltype,
                                                        yield_all_files=args.force,
                                                        picklist=picklist,
                                                        progress=progress)

        # load all matching signatures in this file
        ss = None
        for ss in siglist:
            ksizes.add(ss.minhash.ksize)
            moltypes.add(sourmash_args.get_moltype(ss))
            nums.add(ss.minhash.num)

            if args.scaled:
                ss.minhash = ss.minhash.downsample(scaled=args.scaled)

            if ss.minhash.track_abundance:
                ss.minhash = ss.minhash.flatten()
            scaleds.add(ss.minhash.scaled)

            tree.insert(ss)
            n += 1

        if not ss:
            continue

        # check to make sure we aren't loading incompatible signatures
        if len(ksizes) > 1 or len(moltypes) > 1:
            error('multiple k-mer sizes or molecule types present; fail.')
            error('specify --dna/--protein and --ksize as necessary')
            error('ksizes: {}; moltypes: {}',
                  ", ".join(map(str, ksizes)), ", ".join(moltypes))
            sys.exit(-1)

        if nums == { 0 } and len(scaleds) == 1:
            pass # good
        elif scaleds == { 0 } and len(nums) == 1:
            pass # also good
        else:
            error('trying to build an SBT with incompatible signatures.')
            error('nums = {}; scaleds = {}', repr(nums), repr(scaleds))
            sys.exit(-1)

    notify('')

    # did we load any!?
    if n == 0:
        error('no signatures found to load into tree!? failing.')
        sys.exit(-1)

    if picklist:
        sourmash_args.report_picklist(args, picklist)

    notify(f'loaded {n} sigs; saving SBT under "{args.sbt_name}"')
    tree.save(args.sbt_name, sparseness=args.sparseness)
    if tree.storage:
        tree.storage.close()


def search(args):
    from .search import (search_databases_with_flat_query,
                         search_databases_with_abund_query)

    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)
    picklist = sourmash_args.load_picklist(args)

    # set up the query.
    query = sourmash_args.load_query_signature(args.query,
                                               ksize=args.ksize,
                                               select_moltype=moltype,
                                               select_md5=args.md5)
    notify(f'loaded query: {str(query)[:30]}... (k={query.minhash.ksize}, {sourmash_args.get_moltype(query)})')

    if args.scaled:
        if not query.minhash.scaled:
            error('cannot downsample a signature not created with --scaled')
            sys.exit(-1)
        if args.scaled != query.minhash.scaled:
            notify(f'downsampling query from scaled={query.minhash.scaled} to {int(args.scaled)}')
        query.minhash = query.minhash.downsample(scaled=args.scaled)

    # set up the search databases
    is_containment = args.containment or args.max_containment
    if is_containment:
        if args.containment and args.max_containment:
            notify("ERROR: cannot specify both --containment and --max-containment!")
            sys.exit(-1)

    databases = sourmash_args.load_dbs_and_sigs(args.databases, query,
                                                not is_containment,
                                                picklist=picklist)

    if not len(databases):
        error('Nothing found to search!')
        sys.exit(-1)

    # forcibly ignore abundances if query has no abundances
    if not query.minhash.track_abundance:
        args.ignore_abundance = True
    else:
        if args.ignore_abundance:
            query.minhash = query.minhash.flatten()

    # do the actual search
    if query.minhash.track_abundance:
        results = search_databases_with_abund_query(query, databases,
                                   threshold=args.threshold,
                                   do_containment=args.containment,
                                   do_max_containment=args.max_containment,
                                   best_only=args.best_only,
                                   unload_data=True)
    else:
        results = search_databases_with_flat_query(query, databases,
                                   threshold=args.threshold,
                                   do_containment=args.containment,
                                   do_max_containment=args.max_containment,
                                   best_only=args.best_only,
                                   unload_data=True)

    n_matches = len(results)
    if args.best_only:
        args.num_results = 1

    if not args.num_results or n_matches <= args.num_results:
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
        name = sr.match._display_name(60)
        print_results('{:>6}       {}', pct, name)

    if args.best_only:
        notify("** reporting only one match because --best-only was set")

    if args.output:
        fieldnames = ['similarity', 'name', 'filename', 'md5',
                      'query_filename', 'query_name', 'query_md5']

        with FileOutputCSV(args.output) as fp:
            w = csv.DictWriter(fp, fieldnames=fieldnames)

            w.writeheader()
            for sr in results:
                d = dict(sr._asdict())
                del d['match']
                del d['query']
                w.writerow(d)

    # save matching signatures upon request
    if args.save_matches:
        notify(f'saving all matched signatures to "{args.save_matches}"')

        with SaveSignaturesToLocation(args.save_matches) as save_sig:
            for sr in results:
                save_sig.add(sr.match)

    if picklist:
        sourmash_args.report_picklist(args, picklist)


def categorize(args):
    "Use a database to find the best match to many signatures."
    from .index import MultiIndex
    from .search import make_jaccard_search_query

    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    # eliminate names we've already categorized
    already_names = set()
    if args.load_csv:
        with open(args.load_csv, newline='') as fp:
            r = csv.reader(fp)
            for row in r:
                already_names.add(row[0])

    # load search database
    db = sourmash_args.load_file_as_index(args.database)
    if args.ksize or moltype:
        db = db.select(ksize=args.ksize, moltype=moltype)

    # utility function to load & select relevant signatures.
    def _yield_all_sigs(queries, ksize, moltype):
        for filename in queries:
            mi = MultiIndex.load_from_path(filename, False)
            mi = mi.select(ksize=ksize, moltype=moltype)
            for ss, loc in mi.signatures_with_location():
                yield ss, loc

    csv_w = None
    csv_fp = None
    if args.csv:
        csv_fp = open(args.csv, 'w', newline='')
        csv_w = csv.writer(csv_fp)

    search_obj = make_jaccard_search_query(threshold=args.threshold)
    for orig_query, loc in _yield_all_sigs(args.queries, args.ksize, moltype):
        # skip if we've already done signatures from this file.
        if loc in already_names:
            continue

        notify(f'loaded query: {str(orig_query)[:30]}... (k={orig_query.minhash.ksize}, {orig_query.minhash.moltype})')

        if args.ignore_abundance:
            query = orig_query.copy()
            query.minhash = query.minhash.flatten()
        else:
            if orig_query.minhash.track_abundance:
                notify("ERROR: this search cannot be done on signatures calculated with abundance.")
                notify("ERROR: please specify --ignore-abundance.")
                sys.exit(-1)

            query = orig_query

        results = []
        for sr in db.find(search_obj, query):
            match = sr.signature
            if match.md5sum() != query.md5sum(): # ignore self.
                results.append((orig_query.similarity(match), match))

        if results:
            results.sort(key=lambda x: -x[0])   # reverse sort on similarity
            best_hit_sim, best_hit_query = results[0]
            notify(f'for {query}, found: {best_hit_sim:.2f} {best_hit_query}')
            best_hit_query_name = best_hit_query.name
            if csv_w:
                csv_w.writerow([loc, query, best_hit_query_name,
                               best_hit_sim])
        else:
            notify(f'for {query}, no match found')

    if csv_fp:
        csv_fp.close()


def gather(args):
    from .search import GatherDatabases, format_bp

    set_quiet(args.quiet, args.debug)
    moltype = sourmash_args.calculate_moltype(args)
    picklist = sourmash_args.load_picklist(args)

    # load the query signature & figure out all the things
    query = sourmash_args.load_query_signature(args.query,
                                               ksize=args.ksize,
                                               select_moltype=moltype,
                                               select_md5=args.md5)
    notify(f'loaded query: {str(query)[:30]}... (k={query.minhash.ksize}, {sourmash_args.get_moltype(query)})')

    # verify signature was computed right.
    if not query.minhash.scaled:
        error('query signature needs to be created with --scaled')
        sys.exit(-1)

    if args.scaled:
        notify(f'downsampling query from scaled={query.minhash.scaled} to {int(args.scaled)}')
        query.minhash = query.minhash.downsample(scaled=args.scaled)

    # empty?
    if not len(query.minhash):
        error('no query hashes!? exiting.')
        sys.exit(-1)

    # set up the search databases
    cache_size = args.cache_size
    if args.cache_size == 0:
        cache_size = None
    databases = sourmash_args.load_dbs_and_sigs(args.databases, query, False,
                                                cache_size=cache_size,
                                                picklist=picklist)

    if not len(databases):
        error('Nothing found to search!')
        sys.exit(-1)

    if args.linear:             # force linear traversal?
        databases = [ LazyLinearIndex(db) for db in databases ]

    if args.prefetch:           # note: on by default!
        notify("Starting prefetch sweep across databases.")
        prefetch_query = query.copy()
        prefetch_query.minhash = prefetch_query.minhash.flatten()
        noident_mh = prefetch_query.minhash.to_mutable()
        save_prefetch = SaveSignaturesToLocation(args.save_prefetch)
        save_prefetch.open()
        # set up prefetch CSV output, write headers, etc.
        prefetch_csvout_fp = None
        prefetch_csvout_w = None
        if args.save_prefetch_csv:
            fieldnames = ['intersect_bp', 'jaccard',
                          'max_containment', 'f_query_match', 'f_match_query',
                          'match_filename', 'match_name', 'match_md5', 'match_bp',
                          'query_filename', 'query_name', 'query_md5', 'query_bp']

            prefetch_csvout_fp = FileOutput(args.save_prefetch_csv, 'wt').open()
            prefetch_csvout_w = csv.DictWriter(prefetch_csvout_fp, fieldnames=fieldnames)
            prefetch_csvout_w.writeheader()

        counters = []
        for db in databases:
            counter = None
            try:
                counter, prefetch_info = db.counter_gather(prefetch_query, args.threshold_bp, args.save_prefetch_csv)
            except ValueError:
                if picklist:
                    # catch "no signatures to search" ValueError...
                    continue
                else:
                    raise       # re-raise other errors, if no picklist.

            save_prefetch.add_many(counter.siglist)

            if prefetch_csvout_fp:
                for prefetch_result in prefetch_info:
                    prefetch_csvout_w.writerow(prefetch_result)

            # subtract found hashes as we can.
            for found_sig in counter.siglist:
                noident_mh.remove_many(found_sig.minhash)
            counters.append(counter)

        notify(f"Found {len(save_prefetch)} signatures via prefetch; now doing gather.")
        save_prefetch.close()
        if prefetch_csvout_fp:
            prefetch_csvout_fp.close()
    else:
        counters = databases
        # we can't track unidentified hashes w/o prefetch
        noident_mh = None

    ## ok! now do gather -

    found = []
    weighted_missed = 1
    is_abundance = query.minhash.track_abundance and not args.ignore_abundance
    orig_query_mh = query.minhash
    gather_iter = GatherDatabases(query, counters,
                                  threshold_bp=args.threshold_bp,
                                  ignore_abundance=args.ignore_abundance,
                                  noident_mh=noident_mh)

    for result, weighted_missed in gather_iter:
        if not len(found):                # first result? print header.
            if is_abundance:
                print_results("")
                print_results("overlap     p_query p_match avg_abund")
                print_results("---------   ------- ------- ---------")
            else:
                print_results("")
                print_results("overlap     p_query p_match")
                print_results("---------   ------- -------")


        # print interim result & save in `found` list for later use
        pct_query = '{:.1f}%'.format(result.f_unique_weighted*100)
        pct_genome = '{:.1f}%'.format(result.f_match*100)
        name = result.match._display_name(40)

        if is_abundance:
            average_abund ='{:.1f}'.format(result.average_abund)
            print_results('{:9}   {:>7} {:>7} {:>9}    {}',
                      format_bp(result.intersect_bp), pct_query, pct_genome,
                      average_abund, name)
        else:
            print_results('{:9}   {:>7} {:>7}    {}',
                      format_bp(result.intersect_bp), pct_query, pct_genome,
                      name)
        found.append(result)

        if args.num_results and len(found) >= args.num_results:
            break


    # report on thresholding -
    if gather_iter.query:
        # if still a query, then we failed the threshold.
        notify(f'found less than {format_bp(args.threshold_bp)} in common. => exiting')

    # basic reporting:
    print_results(f'\nfound {len(found)} matches total;')
    if args.num_results and len(found) == args.num_results:
        print_results(f'(truncated gather because --num-results={args.num_results})')

    p_covered = (1 - weighted_missed) * 100
    print_results(f'the recovered matches hit {p_covered:.1f}% of the query')
    print_results('')
    if gather_iter.scaled != query.minhash.scaled:
        print_results(f'WARNING: final scaled was {gather_iter.scaled}, vs query scaled of {query.minhash.scaled}')

    # save CSV?
    if found and args.output:
        fieldnames = ['intersect_bp', 'f_orig_query', 'f_match',
                      'f_unique_to_query', 'f_unique_weighted',
                      'average_abund', 'median_abund', 'std_abund', 'name',
                      'filename', 'md5', 'f_match_orig', 'unique_intersect_bp',
                      'gather_result_rank', 'remaining_bp',
                      'query_filename', 'query_name', 'query_md5', 'query_bp']

        with FileOutputCSV(args.output) as fp:
            w = csv.DictWriter(fp, fieldnames=fieldnames)
            w.writeheader()
            for result in found:
                d = dict(result._asdict())
                del d['match']                 # actual signature not in CSV.
                w.writerow(d)

    # save matching signatures?
    if found and args.save_matches:
        notify(f"saving all matches to '{args.save_matches}'")
        with SaveSignaturesToLocation(args.save_matches) as save_sig:
            for sr in found:
                save_sig.add(sr.match)

    # save unassigned hashes?
    if args.output_unassigned:
        remaining_query = gather_iter.query
        if not (remaining_query.minhash or noident_mh):
            notify('no unassigned hashes to save with --output-unassigned!')
        else:
            notify(f"saving unassigned hashes to '{args.output_unassigned}'")

            if noident_mh:
                remaining_mh = remaining_query.minhash.to_mutable()
                remaining_mh += noident_mh
                remaining_query.minhash = remaining_mh

            if is_abundance:
                abund_query_mh = remaining_query.minhash.inflate(orig_query_mh)
                remaining_query.minhash = abund_query_mh

            with FileOutput(args.output_unassigned, 'wt') as fp:
                sig.save_signatures([ remaining_query ], fp)

    if picklist:
        sourmash_args.report_picklist(args, picklist)

    # DONE w/gather function.


def multigather(args):
    "Gather many signatures against multiple databases."
    from .search import GatherDatabases, format_bp

    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    if not args.db:
        error('Error! must specify at least one database with --db')
        sys.exit(-1)

    if not args.query and not args.query_from_file:
        error('Error! must specify at least one query signature with --query')
        sys.exit(-1)

    # flatten --db and --query
    args.db = [item for sublist in args.db for item in sublist]
    inp_files = [item for sublist in args.query for item in sublist]
    if args.query_from_file:
        more_files = sourmash_args.load_pathlist_from_file(args.query_from_file)
        inp_files.extend(more_files)

    # need a query to get ksize, moltype for db loading
    query = next(iter(sourmash_args.load_file_as_signatures(inp_files[0], ksize=args.ksize, select_moltype=moltype)))

    databases = sourmash_args.load_dbs_and_sigs(args.db, query, False)

    if not len(databases):
        error('Nothing found to search!')
        sys.exit(-1)

    # run gather on all the queries.
    n=0
    for queryfile in inp_files:
        # load the query signature(s) & figure out all the things
        for query in sourmash_args.load_file_as_signatures(queryfile,
                                                       ksize=args.ksize,
                                                       select_moltype=moltype):
            notify(f'loaded query: {str(query)[:30]}... (k={query.minhash.ksize}, {sourmash_args.get_moltype(query)})')

            # verify signature was computed right.
            if not query.minhash.scaled:
                error('query signature needs to be created with --scaled; skipping')
                continue

            if args.scaled:
                notify(f'downsampling query from scaled={query.minhash.scaled} to {int(args.scaled)}')
                query.minhash = query.minhash.downsample(scaled=args.scaled)

            # empty?
            if not len(query.minhash):
                error('no query hashes!? skipping to next..')
                continue

            counters = []
            prefetch_query = query.copy()
            prefetch_query.minhash = prefetch_query.minhash.flatten()
            noident_mh = prefetch_query.minhash.to_mutable()

            counters = []
            for db in databases:
                counter = db.counter_gather(prefetch_query, args.threshold_bp)
                for found_sig in counter.siglist:
                    noident_mh.remove_many(found_sig.minhash)
                counters.append(counter)

            found = []
            weighted_missed = 1
            is_abundance = query.minhash.track_abundance and not args.ignore_abundance
            gather_iter = GatherDatabases(query, counters,
                                          threshold_bp=args.threshold_bp,
                                          ignore_abundance=args.ignore_abundance,
                                          noident_mh=noident_mh)
            for result, weighted_missed in gather_iter:
                if not len(found):                # first result? print header.
                    if is_abundance:
                        print_results("")
                        print_results("overlap     p_query p_match avg_abund")
                        print_results("---------   ------- ------- ---------")
                    else:
                        print_results("")
                        print_results("overlap     p_query p_match")
                        print_results("---------   ------- -------")


                # print interim result & save in a list for later use
                pct_query = '{:.1f}%'.format(result.f_unique_weighted*100)
                pct_genome = '{:.1f}%'.format(result.f_match*100)
                name = result.match._display_name(40)

                if is_abundance:
                    average_abund ='{:.1f}'.format(result.average_abund)
                    print_results('{:9}   {:>7} {:>7} {:>9}    {}',
                              format_bp(result.intersect_bp), pct_query, pct_genome,
                              average_abund, name)
                else:
                    print_results('{:9}   {:>7} {:>7}    {}',
                              format_bp(result.intersect_bp), pct_query, pct_genome,
                              name)
                found.append(result)

            # report on thresholding -
            if gather_iter.query.minhash:
                # if still a query, then we failed the threshold.
                notify(f'found less than {format_bp(args.threshold_bp)} in common. => exiting')

            # basic reporting
            print_results('\nfound {} matches total;', len(found))

            print_results('the recovered matches hit {:.1f}% of the query',
                   (1 - weighted_missed) * 100)
            print_results('')

            if not found:
                notify('nothing found... skipping.')
                continue

            query_filename = query.filename
            if not query_filename:
                # use md5sum if query.filename not properly set
                query_filename = query.md5sum()

            output_base = os.path.basename(query_filename)
            output_csv = output_base + '.csv'

            fieldnames = ['intersect_bp', 'f_orig_query', 'f_match',
                          'f_unique_to_query', 'f_unique_weighted',
                          'average_abund', 'median_abund', 'std_abund', 'name',
                          'filename', 'md5', 'f_match_orig',
                          'unique_intersect_bp', 'gather_result_rank',
                          'remaining_bp', 'query_filename', 'query_name',
                          'query_md5', 'query_bp']
            with FileOutputCSV(output_csv) as fp:
                w = csv.DictWriter(fp, fieldnames=fieldnames)
                w.writeheader()
                for result in found:
                    d = dict(result._asdict())
                    del d['match']      # actual signature not output to CSV!
                    w.writerow(d)

            output_matches = output_base + '.matches.sig'
            with open(output_matches, 'wt') as fp:
                outname = output_matches
                notify(f'saving all matches to "{outname}"')
                sig.save_signatures([ r.match for r in found ], fp)

            output_unassigned = output_base + '.unassigned.sig'
            with open(output_unassigned, 'wt') as fp:
                remaining_query = gather_iter.query
                if noident_mh:
                    remaining_mh = remaining_query.minhash.to_mutable()
                    remaining_mh += noident_mh.downsample(scaled=remaining_mh.scaled)
                    remaining_query.minhash = remaining_mh

                if not found:
                    notify('nothing found - entire query signature unassigned.')
                elif not remaining_query:
                    notify('no unassigned hashes! not saving.')
                else:
                    notify(f'saving unassigned hashes to "{output_unassigned}"')

                    # CTB: note, multigather does not save abundances
                    sig.save_signatures([ remaining_query ], fp)
            n += 1

        # fini, next query!
    notify(f'\nconducted gather searches on {n} signatures')


def watch(args):
    "Build a signature from raw FASTA/FASTQ coming in on stdin, search."
    set_quiet(args.quiet)

    if args.input_is_protein and args.dna:
        notify('WARNING: input is protein, turning off nucleotide hashing.')
        args.dna = False
        args.protein = True

    if args.dna and args.protein:
        notify('ERROR: cannot use "watch" with both nucleotide and protein.')

    if args.dna:
        moltype = 'DNA'
        is_protein = False
        dayhoff = False
        hp = False
    elif args.protein:
        moltype = 'protein'
        is_protein = True
        dayhoff = False
        hp = False
    elif args.dayhoff:
        moltype = 'dayhoff'
        is_protein = True
        dayhoff = True
        hp = False
    else:
        moltype = 'hp'
        is_protein = True
        dayhoff = False
        hp = True

    tree = load_sbt_index(args.sbt_name)

    # check ksize from the SBT we are loading
    ksize = args.ksize
    if ksize is None:
        leaf = next(iter(tree.leaves()))
        tree_mh = leaf.data.minhash
        ksize = tree_mh.ksize

    E = MinHash(ksize=ksize, n=args.num_hashes, is_protein=is_protein, dayhoff=dayhoff, hp=hp)

    notify(f'Computing signature for k={ksize}, {moltype} from stdin')

    def do_search():
        results = []
        streamsig = sig.SourmashSignature(E, filename='stdin', name=args.name)
        for similarity, match, _ in tree.search(streamsig,
                                                threshold=args.threshold,
                                                best_only=True,
                                                ignore_abundance=True,
                                                do_containment=False):
            results.append((similarity, match))

        return results

    notify('reading sequences from stdin')
    screed_iter = screed.open(args.inp_file)
    watermark = WATERMARK_SIZE

    # iterate over input records
    n = 0
    for n, record in enumerate(screed_iter):
        # at each watermark, print status & check cardinality
        if n >= watermark:
            notify(f'\r... read {n} sequences', end='')
            watermark += WATERMARK_SIZE

            if do_search():
                break

        if args.input_is_protein:
            E.add_protein(record.sequence)
        else:
            E.add_sequence(record.sequence, False)

    results = do_search()
    if not results:
        notify(f'... read {n} sequences, no matches found.')
    else:
        results.sort(key=lambda x: -x[0])   # take best
        similarity, found_sig = results[0]
        print_results('FOUND: {}, at {:.3f}', found_sig,
               similarity)

    if args.output:
        notify(f'saving signature to {args.output}')
        with FileOutput(args.output, 'wt') as fp:
            streamsig = sig.SourmashSignature(E, filename='stdin',
                                              name=args.name)
        sig.save_signatures([streamsig], fp)


def migrate(args):
    "Migrate an SBT database to the latest version."
    tree = load_sbt_index(args.sbt_name, print_version_warning=False)

    notify(f'saving SBT under "{args.sbt_name}".')
    tree.save(args.sbt_name, structure_only=True)


def prefetch(args):
    "Output the 'raw' results of a containment/overlap search."

    # load databases from files, too.
    if args.db_from_file:
        more_db = sourmash_args.load_pathlist_from_file(args.db_from_file)
        args.databases.extend(more_db)

    if not args.databases:
        notify("ERROR: no databases or signatures to search!?")
        sys.exit(-1)

    if not (args.save_unmatched_hashes or args.save_matching_hashes or
            args.save_matches or args.output):
        notify("WARNING: no output(s) specified! Nothing will be saved from this prefetch!")

    # figure out what k-mer size and molecule type we're looking for here
    ksize = args.ksize
    moltype = sourmash_args.calculate_moltype(args)
    picklist = sourmash_args.load_picklist(args)

    # load the query signature & figure out all the things
    query = sourmash_args.load_query_signature(args.query,
                                               ksize=args.ksize,
                                               select_moltype=moltype,
                                               select_md5=args.md5)
    notify(f'loaded query: {str(query)[:30]}... (k={query.minhash.ksize}, {sourmash_args.get_moltype(query)})')

    # verify signature was computed with scaled.
    if not query.minhash.scaled:
        error('query signature needs to be created with --scaled')
        sys.exit(-1)

    # if with track_abund, flatten me
    query_mh = query.minhash
    if query_mh.track_abundance:
        query_mh = query_mh.flatten()

    if args.scaled:
        notify(f'downsampling query from scaled={query_mh.scaled} to {int(args.scaled)}')
        query_mh = query_mh.downsample(scaled=args.scaled)
    notify(f"all sketches will be downsampled to scaled={query_mh.scaled}")

    # empty?
    if not len(query_mh):
        error('no query hashes!? exiting.')
        sys.exit(-1)

    query.minhash = query_mh

    # set up CSV output, write headers, etc.
    csvout_fp = None
    csvout_w = None
    if args.output:
        fieldnames = ['intersect_bp', 'jaccard',
                      'max_containment', 'f_query_match', 'f_match_query',
                      'match_filename', 'match_name', 'match_md5', 'match_bp',
                      'query_filename', 'query_name', 'query_md5', 'query_bp']

        csvout_fp = FileOutput(args.output, 'wt').open()
        csvout_w = csv.DictWriter(csvout_fp, fieldnames=fieldnames)
        csvout_w.writeheader()

    # track & maybe save matches progressively
    matches_out = SaveSignaturesToLocation(args.save_matches)
    matches_out.open()
    if args.save_matches:
        notify(f"saving all matching database signatures to '{args.save_matches}'")

    # iterate over signatures in db one at a time, for each db;
    # find those with sufficient overlap
    ident_mh = query_mh.copy_and_clear()
    noident_mh = query_mh.to_mutable()

    did_a_search = False        # track whether we did _any_ search at all!
    for dbfilename in args.databases:
        notify(f"loading signatures from '{dbfilename}'")

        db = sourmash_args.load_file_as_index(dbfilename)

        # force linear traversal?
        if args.linear:
            db = LazyLinearIndex(db)

        db = db.select(ksize=ksize, moltype=moltype,
                       containment=True, scaled=True,
                       picklist=picklist)

        if not db:
            notify(f"...no compatible signatures in '{dbfilename}'; skipping")
            continue

        for result in prefetch_database(query, db, args.threshold_bp):
            match = result.match

            # track found & "untouched" hashes.
            match_mh = match.minhash.downsample(scaled=query.minhash.scaled)
            ident_mh += query.minhash & match_mh.flatten()
            noident_mh.remove_many(match.minhash)

            # output match info as we go
            if csvout_fp:
                d = dict(result._asdict())
                del d['match']                 # actual signatures not in CSV.
                del d['query']
                csvout_w.writerow(d)

            # output match signatures as we go (maybe)
            matches_out.add(match)

            if matches_out.count % 10 == 0:
                notify(f"total of {matches_out.count} matching signatures so far.",
                       end="\r")

        did_a_search = True

        # flush csvout so that things get saved progressively
        if csvout_fp:
            csvout_fp.flush()

        # delete db explicitly ('cause why not)
        del db

    if not did_a_search:
        notify("ERROR in prefetch: no compatible signatures in any databases?!")
        sys.exit(-1)

    notify(f"total of {matches_out.count} matching signatures.")
    matches_out.close()

    if csvout_fp:
        notify(f"saved {matches_out.count} matches to CSV file '{args.output}'")
        csvout_fp.close()

    assert len(query_mh) == len(ident_mh) + len(noident_mh)
    notify(f"of {len(query_mh)} distinct query hashes, {len(ident_mh)} were found in matches above threshold.")
    notify(f"a total of {len(noident_mh)} query hashes remain unmatched.")

    if args.save_matching_hashes:
        filename = args.save_matching_hashes
        notify(f"saving {len(ident_mh)} matched hashes to '{filename}'")

        sig_name = ''
        if query.name:
            sig_name = f"{query.name}-known"

        ss = sig.SourmashSignature(ident_mh, name=sig_name)
        with open(filename, "wt") as fp:
            sig.save_signatures([ss], fp)

    if args.save_unmatched_hashes:
        filename = args.save_unmatched_hashes

        sig_name = ''
        if query.name:
            sig_name = f"{query.name}-unknown"

        notify(f"saving {len(noident_mh)} unmatched hashes to '{filename}'")
        ss = sig.SourmashSignature(noident_mh, name=sig_name)
        with open(filename, "wt") as fp:
            sig.save_signatures([ss], fp)

    if picklist:
        sourmash_args.report_picklist(args, picklist)

    return 0
