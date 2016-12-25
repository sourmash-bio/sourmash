from __future__ import print_function
import sys
import os, os.path
import argparse
import csv

import screed
import sourmash_lib
from sourmash_lib import signature as sig
from sourmash_lib import fig as sourmash_fig
from . import sourmash_args

DEFAULT_K = 31
DEFAULT_N = 500


class SourmashCommands(object):

    def __init__(self):
        parser = argparse.ArgumentParser(description='work with RNAseq signatures',
                                         usage='''sourmash <command> [<args>]

Commands can be:

   compute <filenames>         Compute signatures for sequences in these files.
   compare <filenames.sig>     Compute distance matrix for given signatures.
   search <query> <against>    Search for matching signatures.
   plot <matrix>               Plot a distance matrix made by 'compare'.

   import_csv                  Import signatures from a CSV file.
.
''')
        parser.add_argument('command')
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            print('Unrecognized command')
            parser.print_help()
            sys.exit(1)

        cmd = getattr(self, args.command)
        print('# running sourmash subcommand: %s' % args.command,
              file=sys.stderr)
        cmd(sys.argv[2:])

    def search(self, args):
        "Search a query sig against one or more signatures; report up to the 3 top matches."
        parser = argparse.ArgumentParser()
        parser.add_argument('query')
        parser.add_argument('against', nargs='+')
        parser.add_argument('--threshold', default=0.08, type=float)
        parser.add_argument('-n', '--num-results', default=3, type=int)
        parser.add_argument('-k', '--ksize', default=DEFAULT_K, type=int)
        parser.add_argument('-f', '--force', action='store_true')
        parser.add_argument('--save-matches', type=argparse.FileType('wt'))

        sourmash_args.add_moltype_args(parser)

        args = parser.parse_args(args)

        if args.protein:
            if args.dna is True:
                raise Exception('cannot specify both --dna and --protein!')
            args.dna = False

        moltype = None
        if args.protein:
            moltype = 'protein'
        elif args.dna:
            moltype = 'dna'

        # get the query signature
        sl = sig.load_signatures(args.query,
                                 select_ksize=args.ksize,
                                 select_moltype=moltype)
        sl = list(sl)

        if len(sl) != 1:
            print('When loading query from "{}",'.format(args.query),
                  file=sys.stderr)
            print('{} query signatures matching ksize and molecule type; need exactly one.'.format(len(sl)))
            sys.exit(-1)
        query = sl[0]

        query_moltype = 'UNKNOWN'
        if query.estimator.is_molecule_type('dna'):
            query_moltype = 'DNA'
        elif query.estimator.is_molecule_type('protein'):
            query_moltype = 'protein'
        query_ksize = query.estimator.ksize
        print('loaded query: {}... (k={}, {})'.format(query.name()[:30],
                                                      query_ksize,
                                                      query_moltype))

        # get the signatures to query
        print('loading db of signatures from %d files' % len(args.against),
              file=sys.stderr)
        against = []
        for filename in args.against:
            if filename == args.query and not args.force:
                print('excluding query from database (file %s)' % filename,
                      file=sys.stderr)
                continue

            sl = sig.load_signatures(filename,
                                     select_ksize=args.ksize,
                                     select_moltype=moltype)

            for x in sl:
                against.append((x, filename))
        print('loaded {} signatures total.'.format(len(against)))

        # compute query x db
        distances = []
        for (x, filename) in against:
            distance = query.similarity(x)
            if distance >= args.threshold:
                distances.append((distance, x, filename))

        # any matches? sort, show.
        if distances:
            distances.sort(reverse=True, key = lambda x: x[0])
            print('{} matches; showing {}:'.format(len(distances),
                                                   args.num_results))
            for distance, match, filename in distances[:args.num_results]:

                print('\t', match.name(), '\t', "%.3f" % distance,
                      '\t', filename)

            if args.save_matches:
                outname = args.save_matches.name
                print('saving all matches to "{}"'.format(outname))
                sig.save_signatures([ m for (d, m, f) in distances ],
                                    args.save_matches)
        else:
            print('** no matches in %d signatures' % len(against),
                  file=sys.stderr)

    def compute(self, args):
        """Compute the signature for one or more files.

        Use cases:
            sourmash compute multiseq.fa              => multiseq.fa.sig, etc.
            sourmash compute genome.fa --singleton    => genome.fa.sig
            sourmash compute file1.fa file2.fa --name => specify w/-o
        """
        parser = argparse.ArgumentParser()
        parser.add_argument('filenames', nargs='+')

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
        parser.add_argument('--with-cardinality', action='store_true')
        args = parser.parse_args(args)

        if args.input_is_protein and args.dna:
            print('WARNING: input is protein, turning off DNA hash computing.',
                  file=sys.stderr)
            args.dna = False
            args.protein = True

        if not args.dna and args.protein:
            if args.with_cardinality:
                print('Cannot compute cardinality for protein sequences.',
                       file=sys.stderr)
                sys.exit(-1)

        print('computing signatures for files:', args.filenames,
              file=sys.stderr)

        # get list of k-mer sizes for which to compute sketches
        ksizes = args.ksizes
        if ',' in ksizes:
            ksizes = ksizes.split(',')
            ksizes = list(map(int, ksizes))
        else:
            ksizes = [int(ksizes)]

        print('Computing signature for ksizes: %s' % str(ksizes),
              file=sys.stderr)

        num_sigs = 0
        if args.dna and args.protein:
            print('Computing both DNA and protein signatures.',
                  file=sys.stderr)
            num_sigs = 2*len(ksizes)
        elif args.dna:
            print('Computing only DNA (and not protein) signatures.',
                  file=sys.stderr)
            num_sigs = len(ksizes)
        elif args.protein:
            print('Computing only protein (and not DNA) signatures.',
                  file=sys.stderr)
            num_sigs = len(ksizes)

        if args.protein:
            bad_ksizes = [ str(k) for k in ksizes if k % 3 != 0 ]
            if bad_ksizes:
                print('protein ksizes must be divisible by 3, sorry!',
                      file=sys.stderr)
                print('bad ksizes: {}'.format(", ".join(bad_ksizes)),
                      file=sys.stderr)
                sys.exit(-1)

        print('Computing a total of {} signatures.'.format(num_sigs),
              file=sys.stderr)

        if num_sigs == 0:
            print('...nothing to calculate!? Exiting!', file=sys.stderr)
            sys.exit(-1)

        if args.name and not args.output:
            print("must specify -o with --name", file=sys.stderr)
            sys.exit(-1)

        def make_estimators():
            # one estimator for each ksize
            Elist = []
            for k in ksizes:
                if args.protein:
                    E = sourmash_lib.Estimators(ksize=k, n=args.num_hashes,
                                                protein=True)
                    Elist.append(E)
                if args.dna:
                    E = sourmash_lib.Estimators(ksize=k, n=args.num_hashes,
                                                protein=False,
                                        with_cardinality=args.with_cardinality)
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

        print('Computing signature for ksizes: %s' % str(ksizes),
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

    def compare(self, args):
        "Compare multiple signature files and create a distance matrix."
        import numpy

        parser = argparse.ArgumentParser()
        parser.add_argument('signatures', nargs='+')
        parser.add_argument('-k', '--ksize', type=int, default=DEFAULT_K, help='k-mer size (default: %(default)s)')
        parser.add_argument('-o', '--output')
        args = parser.parse_args(args)

        # load in the various signatures
        siglist = []
        for filename in args.signatures:
            print('loading', filename, file=sys.stderr)
            loaded = sig.load_signatures(filename, select_ksize=args.ksize)
            loaded = list(loaded)
            if not loaded:
                print('warning: no signatures loaded at given ksize from %s' %
                          filename, file=sys.stderr)
            siglist.extend(loaded)

        if len(siglist) == 0:
            print('no signatures!', file=sys.stderr)
            sys.exit(-1)

        # build the distance matrix
        D = numpy.zeros([len(siglist), len(siglist)])
        numpy.set_printoptions(precision=3, suppress=True)

        # do all-by-all calculation
        i = 0
        labeltext = []
        for i, E in enumerate(siglist):
            for j, E2 in enumerate(siglist):
                D[i][j] = E.similarity(E2)

            print('%d-%20s\t%s' % (i, E.name(), D[i, :, ],))
            labeltext.append(E.name())
            i += 1

        print('min similarity in matrix:', numpy.min(D), file=sys.stderr)

        # shall we output a matrix?
        if args.output:
            labeloutname = args.output + '.labels.txt'
            print('saving labels to:', labeloutname, file=sys.stderr)
            with open(labeloutname, 'w') as fp:
                fp.write("\n".join(labeltext))

            print('saving distance matrix to:', args.output,
                  file=sys.stderr)
            with open(args.output, 'wb') as fp:
                numpy.save(fp, D)


    def plot(self, args):
        "Produce a clustering and plot."
        import numpy
        import scipy
        import pylab
        import scipy.cluster.hierarchy as sch

        # set up cmd line arguments
        parser = argparse.ArgumentParser()
        parser.add_argument('distances', help="output from 'sourmash compare'")
        parser.add_argument('--pdf', action='store_true')
        parser.add_argument('--labels', action='store_true')
        parser.add_argument('--indices', action='store_false')
        parser.add_argument('--vmax', default=1.0, type=float, help='(default: %(default)f)')
        parser.add_argument('--vmin', default=0.0, type=float, help='(default: %(default)f)')
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
        print('wrote', dendrogram_out)

        ### make the dendrogram+matrix:
        fig = sourmash_fig.plot_composite_matrix(D, labeltext,
                                                 show_labels=args.labels,
                                                 show_indices=args.indices,
                                                 vmin=args.vmin,
                                                 vmax=args.vmax)
        fig.savefig(matrix_out)
        print('wrote', matrix_out)

        # print out sample numbering for FYI.
        for i, name in enumerate(labeltext):
            print(i, '\t', name)

    def import_csv(self, args):
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
                for h in hashes:
                    e.mh.add_hash(h)
                s = sig.SourmashSignature(args.email, e, filename=name)
                siglist.append(s)
                print('loaded signature:', name,
                      s.md5sum()[:8], file=sys.stderr)

            print('saving %d signatures to YAML' % (len(siglist),),
                  file=sys.stderr)
            sig.save_signatures(siglist, args.output)

    def dump(self, args):
        parser = argparse.ArgumentParser()
        parser.add_argument('filenames', nargs='+')
        parser.add_argument('-k', '--ksize', type=int, default=DEFAULT_K, help='k-mer size (default: %(default)i)')
        args = parser.parse_args(args)

        for filename in args.filenames:
            print('loading', filename)
            siglist = sig.load_signatures(filename, select_ksize=args.ksize)
            siglist = list(siglist)
            assert len(siglist) == 1

            s = siglist[0]

            fp = open(filename + '.dump.txt', 'w')
            fp.write(" ".join((map(str, s.estimator.mh.get_mins()))))
            fp.close()

    def sbt_index(self, args):
        from sourmash_lib.sbt import SBT, GraphFactory
        from sourmash_lib.sbtmh import search_minhashes, SigLeaf

        parser = argparse.ArgumentParser()
        parser.add_argument('sbt_name')
        parser.add_argument('signatures', nargs='+')
        parser.add_argument('-k', '--ksize', type=int, default=DEFAULT_K)
        parser.add_argument('--traverse-directory', action='store_true')
        parser.add_argument('-x', '--bf-size', type=float, default=1e5)

        sourmash_args.add_moltype_args(parser)

        args = parser.parse_args(args)

        if args.protein:
            if args.dna is True:
                raise Exception('cannot specify both --dna and --protein!')
            args.dna = False
            moltype = 'protein'
        else:
            args.dna = True
            moltype = 'dna'

        factory = GraphFactory(1, args.bf_size, 4)
        tree = SBT(factory)

        inp_files = list(args.signatures)

        if args.traverse_directory:
            inp_files = []
            for dirname in args.signatures:
                for root, dirs, files in os.walk(dirname):
                    for name in files:
                        if name.endswith('.sig'):
                            fullname = os.path.join(root, name)
                            inp_files.append(fullname)

        print('loading {} files into SBT'.format(len(inp_files)))

        n = 0
        for f in inp_files:
            s = sig.load_signatures(f, select_ksize=args.ksize,
                                    select_moltype=moltype)

            for ss in s:
                leaf = SigLeaf(ss.md5sum(), ss)
                tree.add_node(leaf)
                n += 1

        print('loaded {} sigs; saving SBT under "{}".'.format(n,
                                                              args.sbt_name))
        tree.save(args.sbt_name)

    def sbt_search(self, args):
        from sourmash_lib.sbt import SBT, GraphFactory
        from sourmash_lib.sbtmh import search_minhashes, SigLeaf
        from sourmash_lib.sbtmh import SearchMinHashesFindBest

        parser = argparse.ArgumentParser()
        parser.add_argument('sbt_name')
        parser.add_argument('query')
        parser.add_argument('-k', '--ksize', type=int, default=DEFAULT_K)
        parser.add_argument('--threshold', default=0.08, type=float)
        parser.add_argument('--save-matches', type=argparse.FileType('wt'))
        parser.add_argument('--best-only', action='store_true')

        sourmash_args.add_moltype_args(parser)
        args = parser.parse_args(args)

        if args.protein:
            if args.dna is True:
                raise Exception('cannot specify both --dna and --protein!')
            args.dna = False

        moltype = None
        if args.protein:
            moltype = 'protein'
        elif args.dna:
            moltype = 'dna'

        search_fn = search_minhashes
        if args.best_only:
            search_fn = SearchMinHashesFindBest().search

        tree = SBT.load(args.sbt_name, leaf_loader=SigLeaf.load)
        sl = sig.load_signatures(args.query, select_ksize=args.ksize,
                                 select_moltype=moltype)
        sl = list(sl)
        if len(sl) != 1:
            print('When loading query from "{}",'.format(args.query),
                  file=sys.stderr)
            print('{} query signatures matching ksize and molecule type; need exactly one.'.format(len(sl)))
            sys.exit(-1)

        query = sl[0]

        query_moltype = 'UNKNOWN'
        if query.estimator.is_molecule_type('dna'):
            query_moltype = 'DNA'
        elif query.estimator.is_molecule_type('protein'):
            query_moltype = 'protein'
        query_ksize = query.estimator.ksize
        print('loaded query: {}... (k={}, {})'.format(query.name()[:30],
                                                      query_ksize,
                                                      query_moltype))

        results = []
        for leaf in tree.find(search_fn, query, args.threshold):
            results.append((query.similarity(leaf.data), leaf.data))
            #results.append((leaf.data.similarity(ss), leaf.data))

        results.sort(key=lambda x: -x[0])   # reverse sort on similarity
        for (similarity, query) in results:
            print('{:.2f} {}'.format(similarity, query.name()))

        if args.save_matches:
            outname = args.save_matches.name
            print('saving all matches to "{}"'.format(outname))
            sig.save_signatures([ m for (sim, m) in results ],
                                args.save_matches)


    def categorize(self, args):
        from sourmash_lib.sbt import SBT, GraphFactory
        from sourmash_lib.sbtmh import search_minhashes, SigLeaf
        from sourmash_lib.sbtmh import SearchMinHashesFindBest

        parser = argparse.ArgumentParser()
        parser.add_argument('sbt_name')
        parser.add_argument('queries', nargs='+')
        parser.add_argument('-k', '--ksize', type=int, default=DEFAULT_K)
        parser.add_argument('--threshold', default=0.08, type=float)
        parser.add_argument('--traverse-directory', action="store_true")

        sourmash_args.add_moltype_args(parser)

        parser.add_argument('--csv', type=argparse.FileType('at'))
        parser.add_argument('--load-csv', default=None)
        
        args = parser.parse_args(args)

        if args.protein:
            if args.dna is True:
                raise Exception('cannot specify both --dna and --protein!')
            args.dna = False

        moltype = None
        if args.protein:
            moltype = 'protein'
        elif args.dna:
            moltype = 'dna'

        already_names = set()
        if args.load_csv:
            with open(args.load_csv, 'rt') as fp:
                r = csv.reader(fp)
                for row in r:
                    already_names.add(row[0])

        search_fn = SearchMinHashesFindBest().search

        tree = SBT.load(args.sbt_name, leaf_loader=SigLeaf.load)

        if args.traverse_directory:
            inp_files = []
            for dirname in args.queries:
                for root, dirs, files in os.walk(dirname):
                    for name in files:
                        if name.endswith('.sig'):
                            fullname = os.path.join(root, name)
                            inp_files.append(fullname)
        else:
            inp_files = args.queries

        print('found {} files to query'.format(len(inp_files)))

        n_skipped = 0
        for queryfile in inp_files:
            if queryfile in already_names:
                n_skipped += 1
                continue

            sl = sig.load_signatures(queryfile, select_ksize=args.ksize,
                                     select_moltype=moltype)
            sl = list(sl)
            if len(sl) != 1:
                print('When loading query from "{}",'.format(queryfile),
                      file=sys.stderr)
                print('{} query signatures matching ksize and molecule type; need exactly one.'.format(len(sl)))
                continue

            query = sl[0]

            query_moltype = 'UNKNOWN'
            if query.estimator.is_molecule_type('dna'):
                query_moltype = 'DNA'
            elif query.estimator.is_molecule_type('protein'):
                query_moltype = 'protein'
            query_ksize = query.estimator.ksize
            print('loaded query: {}... (k={}, {})'.format(query.name()[:30],
                                                          query_ksize,
                                                          query_moltype))

            results = []
            for leaf in tree.find(search_fn, query, args.threshold):
                # ignore self
                if leaf.data.md5sum() != query.md5sum():
                    results.append((query.similarity(leaf.data), leaf.data))

            best_hit_sim = 0.0
            best_hit_query_name = ""
            if results:
                results.sort(key=lambda x: -x[0])   # reverse sort on similarity
                best_hit_sim, best_hit_query = results[0]
                print('for {}, found: {:.2f} {}'.format(query.name(),
                                                        best_hit_sim,
                                                        best_hit_query.name()))
                best_hit_query_name = best_hit_query.name()
            else:
                print('for {}, no match found'.format(query.name()))

            if args.csv:
                w = csv.writer(args.csv)
                w.writerow([queryfile, best_hit_query_name, best_hit_sim])

        print('note: skipped {}'.format(n_skipped))

    def sbt_gather(self, args):
        from sourmash_lib.sbt import SBT, GraphFactory
        from sourmash_lib.sbtmh import search_minhashes, SigLeaf
        from sourmash_lib.sbtmh import SearchMinHashesFindBest

        parser = argparse.ArgumentParser()
        parser.add_argument('sbt_name')
        parser.add_argument('query')
        parser.add_argument('-k', '--ksize', type=int, default=DEFAULT_K)
        parser.add_argument('--threshold', default=0.05, type=float)
        parser.add_argument('-o', '--output', type=argparse.FileType('wt'))
        parser.add_argument('--csv', type=argparse.FileType('wt'))

        sourmash_args.add_moltype_args(parser)

        args = parser.parse_args(args)

        if args.protein:
            if args.dna is True:
                raise Exception('cannot specify both --dna and --protein!')
            args.dna = False

        moltype = None
        if args.protein:
            moltype = 'protein'
        elif args.dna:
            moltype = 'dna'

        tree = SBT.load(args.sbt_name, leaf_loader=SigLeaf.load)
        sl = sig.load_signatures(args.query, select_ksize=args.ksize,
                                 select_moltype=moltype)
        sl = list(sl)
        if len(sl) != 1:
            print('When loading query from "{}",'.format(args.query),
                  file=sys.stderr)
            print('{} query signatures matching ksize and molecule type; need exactly one.'.format(len(sl)))
            sys.exit(-1)

        query = sl[0]

        query_moltype = 'UNKNOWN'
        if query.estimator.is_molecule_type('dna'):
            query_moltype = 'DNA'
        elif query.estimator.is_molecule_type('protein'):
            query_moltype = 'protein'
        query_ksize = query.estimator.ksize
        print('loaded query: {}... (k={}, {})'.format(query.name()[:30],
                                                      query_ksize,
                                                      query_moltype))

        tree = SBT.load(args.sbt_name, leaf_loader=SigLeaf.load)
        s = sig.load_signatures(args.query, select_ksize=args.ksize)
        orig_query = query

        sum_found = 0.
        found = []
        while 1:
            search_fn = SearchMinHashesFindBest().search

            results = []
            # use super low threshold for this part of the search
            for leaf in tree.find(search_fn, query, 0.00001):
                results.append((query.similarity(leaf.data), leaf.data))
                #results.append((leaf.data.similarity(ss), leaf.data))

            if not len(results):          # no matches at all!
                break

            # take the best result
            results.sort(key=lambda x: -x[0])   # reverse sort on similarity
            best_sim, best_ss = results[0]
            sim = best_ss.similarity(orig_query)

            # adjust by size of leaf (kmer cardinality of original genome)
            if best_ss.estimator.hll:
                leaf_kmers = best_ss.estimator.hll.estimate_cardinality()
                query_kmers = orig_query.estimator.hll.estimate_cardinality()
                f_of_total = leaf_kmers / query_kmers * sim
            else:
                f_of_total = 0

            if not found and sim < args.threshold:
                print('best match: {}'.format(best_ss.name()))
                print('similarity is {:.5f} of db signature;'.format(sim))
                print('this is below specified threshold => exiting.')
                break

            # subtract found hashes from search hashes, construct new search
            new_mins = set(query.estimator.mh.get_mins())
            found_mins = best_ss.estimator.mh.get_mins()

            # print interim & save
            print('found: {:.2f} {} {}'.format(f_of_total,
                                               len(new_mins),
                                               best_ss.name()))
            found.append((f_of_total, best_ss, sim))
            sum_found += f_of_total

            new_mins -= set(found_mins)
            e = sourmash_lib.Estimators(ksize=args.ksize, n=len(new_mins))
            for m in new_mins:
                e.mh.add_hash(m)
            new_ss = sig.SourmashSignature('foo', e)
            query = new_ss

        print('found {}, total fraction {:.3f}'.format(len(found), sum_found))
        print('')

        if not found:
            sys.exit(0)

        found.sort()
        found.reverse()

        print('Composition:')
        for (frac, leaf_sketch, sim) in found:
            print('{:.2f} {}'.format(frac, leaf_sketch.name()))

        if args.output:
            print('Composition:', file=args.output)
            for (frac, leaf_sketch, sim) in found:
                print('{:.2f} {}'.format(frac, leaf_sketch.name()),
                      file=args.output)

        if args.csv:
            fieldnames = ['fraction', 'name', 'similarity', 'sketch_kmers']
            w = csv.DictWriter(args.csv, fieldnames=fieldnames)

            w.writeheader()
            for (frac, leaf_sketch, sim) in found:
                cardinality = leaf_sketch.estimator.hll.estimate_cardinality()
                w.writerow(dict(fraction=frac, name=leaf_sketch.name(),
                                similarity=sim,
                                sketch_kmers=cardinality))


def main():
    SourmashCommands()
    return 0
