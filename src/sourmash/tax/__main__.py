"""
Command-line entry point for 'python -m sourmash.tax'
"""
import sys
import csv
import os
from collections import defaultdict, Counter
import re

import sourmash
from ..sourmash_args import FileOutputCSV, FileOutput
from sourmash.logging import set_quiet, error, notify, print_results
from sourmash.lca.lca_utils import display_lineage, zip_lineage

from . import tax_utils
from .tax_utils import ClassificationResult, MultiLineageDB

usage='''
sourmash taxonomy <command> [<args>] - manipulate/work with taxonomy information.
or
sourmash tax <command> [<args>]


** Commands can be:

annotate -g <gather_csv> [<gather_csv> ... ] -t [<taxonomy_csv> ...]      - annotate gather CSVs with taxonomic lineages
genome -g <gather_csv> [<gather_csv> ... ] -t [<taxonomy_csv> ...]        - taxonomic classification of genomes from gather results
metagenome -g <gather_csv> [<gather_csv> ... ] -t [<taxonomy_csv> ...]    - summarize taxonomic information for metagenome gather results

** Use '-h' to get subcommand-specific help, e.g.

sourmash taxonomy metagenome -h
'''

_output_type_to_ext = {
    'csv_summary': '.summarized.csv',
    'classification': '.classifications.csv',
    'krona': '.krona.tsv',
    'lineage_summary': '.lineage_summary.tsv',
    'annotate': '.with-lineages.csv',
    'human': '.human.txt',
    'lineage_csv': '.lineage.csv',
    'kreport': ".kreport.txt",
    }

# some utils
def make_outfile(base, output_type, *, output_dir = ""):
    limit_float_decimals=False
    if base == "-":
        limit_float_decimals=True
        return base, limit_float_decimals

    ext = _output_type_to_ext[output_type]

    fname = base+ext
    if output_dir:
        fname = os.path.join(output_dir, fname)
    notify(f"saving '{output_type}' output to '{fname}'.")
    return fname, limit_float_decimals


##### taxonomy command line functions
def metagenome(args):
    """
    summarize taxonomic information for metagenome gather results
    """
    set_quiet(args.quiet)

    # first, load taxonomic_assignments
    try:
        tax_assign = MultiLineageDB.load(args.taxonomy_csv,
                       keep_full_identifiers=args.keep_full_identifiers,
                       keep_identifier_versions=args.keep_identifier_versions,
                       force=args.force)
        available_ranks = tax_assign.available_ranks
    except ValueError as exc:
        error(f"ERROR: {str(exc)}")
        sys.exit(-1)

    if not tax_assign:
        error(f'ERROR: No taxonomic assignments loaded from {",".join(args.taxonomy_csv)}. Exiting.')
        sys.exit(-1)

    if args.rank and args.rank not in available_ranks:
        error(f"ERROR: No taxonomic information provided for rank {args.rank}: cannot summarize at this rank")
        sys.exit(-1)

    # next, collect and load gather results
    gather_csvs = tax_utils.collect_gather_csvs(args.gather_csv, from_file= args.from_file)
    try:
        gather_results, idents_missed, total_missed, _ = tax_utils.check_and_load_gather_csvs(gather_csvs, tax_assign, force=args.force,
                                                                                       fail_on_missing_taxonomy=args.fail_on_missing_taxonomy)
    except ValueError as exc:
        error(f"ERROR: {str(exc)}")
        sys.exit(-1)

    if not gather_results:
        notify('No gather results loaded. Exiting.')
        sys.exit(-1)

    # actually summarize at rank
    summarized_gather = {}
    seen_perfect = set()
    for rank in sourmash.lca.taxlist(include_strain=False):
        try:
            summarized_gather[rank], seen_perfect, _ = tax_utils.summarize_gather_at(rank, tax_assign, gather_results, skip_idents=idents_missed,
                                                                keep_full_identifiers=args.keep_full_identifiers,
                                                                keep_identifier_versions = args.keep_identifier_versions,
                                                                seen_perfect = seen_perfect)

        except ValueError as exc:
            error(f"ERROR: {str(exc)}")
            sys.exit(-1)

    # write summarized output csv
    if "csv_summary" in args.output_format:
        summary_outfile, limit_float = make_outfile(args.output_base, "csv_summary", output_dir=args.output_dir)
        with FileOutputCSV(summary_outfile) as out_fp:
            tax_utils.write_summary(summarized_gather, out_fp, limit_float_decimals=limit_float)

    # write summarized output in human-readable format
    if "human" in args.output_format:
        summary_outfile, limit_float = make_outfile(args.output_base, "human", output_dir=args.output_dir)

        with FileOutput(summary_outfile) as out_fp:
            tax_utils.write_human_summary(summarized_gather, out_fp, args.rank or "species")

    # if lineage summary table
    if "lineage_summary" in args.output_format:
        lineage_outfile, limit_float = make_outfile(args.output_base, "lineage_summary", output_dir=args.output_dir)

        ## aggregate by lineage, by query
        lineageD, query_names, num_queries = tax_utils.aggregate_by_lineage_at_rank(summarized_gather[args.rank], by_query=True)

        with FileOutputCSV(lineage_outfile) as out_fp:
            tax_utils.write_lineage_sample_frac(query_names, lineageD, out_fp, format_lineage=True, sep='\t')

    # write summarized --> krona output tsv
    if "krona" in args.output_format:
        krona_resultslist = tax_utils.format_for_krona(args.rank, summarized_gather)

        krona_outfile, limit_float = make_outfile(args.output_base, "krona", output_dir=args.output_dir)
        with FileOutputCSV(krona_outfile) as out_fp:
            tax_utils.write_krona(args.rank, krona_resultslist, out_fp)

    # write summarized --> kreport output tsv
    if "kreport" in args.output_format:
        kreport_outfile, limit_float = make_outfile(args.output_base, "kreport", output_dir=args.output_dir)

        with FileOutputCSV(kreport_outfile) as out_fp:
            tax_utils.write_kreport(summarized_gather, out_fp)


def genome(args):
    """
    taxonomic classification of genomes from gather results
    """
    set_quiet(args.quiet)

    # first, load taxonomic_assignments
    try:
        tax_assign = MultiLineageDB.load(args.taxonomy_csv,
                       keep_full_identifiers=args.keep_full_identifiers,
                       keep_identifier_versions=args.keep_identifier_versions,
                       force=args.force)
        available_ranks = tax_assign.available_ranks
    except ValueError as exc:
        error(f"ERROR: {str(exc)}")
        sys.exit(-1)

    if not tax_assign:
        error(f'ERROR: No taxonomic assignments loaded from {",".join(args.taxonomy_csv)}. Exiting.')
        sys.exit(-1)

    if args.rank and args.rank not in available_ranks:
        error(f"ERROR: No taxonomic information provided for rank {args.rank}: cannot classify at this rank")
        sys.exit(-1)

    # get gather_csvs from args
    gather_csvs = tax_utils.collect_gather_csvs(args.gather_csv, from_file=args.from_file)

    classifications = defaultdict(list)
    matched_queries=set()
    krona_results = []
    status = "nomatch"
    seen_perfect = set()

    # read in all gather CSVs (queries in more than one gather file will raise error; with --force they will only be loaded once)
    # note: doing one CSV at a time would work and probably be more memory efficient, but we would need to change how we check
    # for duplicated queries
    try:
        gather_results, idents_missed, total_missed, _ = tax_utils.check_and_load_gather_csvs(gather_csvs, tax_assign, force=args.force,
                                                                            fail_on_missing_taxonomy=args.fail_on_missing_taxonomy)

    except ValueError as exc:
        error(f"ERROR: {str(exc)}")
        sys.exit(-1)

    # if --rank is specified, classify to that rank
    estimate_query_ani = True
    if args.rank:
        try:
            best_at_rank, seen_perfect, estimate_query_ani = tax_utils.summarize_gather_at(args.rank, tax_assign, gather_results, skip_idents=idents_missed,
                                                     keep_full_identifiers=args.keep_full_identifiers,
                                                     keep_identifier_versions = args.keep_identifier_versions,
                                                     best_only=True, seen_perfect=seen_perfect, estimate_query_ani=True)

        except ValueError as exc:
            error(f"ERROR: {str(exc)}")
            sys.exit(-1)

       # best at rank is a list of SummarizedGather tuples
        for sg in best_at_rank:
            status = 'nomatch'
            if sg.query_name in matched_queries:
                continue
            if args.ani_threshold and sg.query_ani_at_rank < args.ani_threshold:
                status="below_threshold"
                notify(f"WARNING: classifying query {sg.query_name} at desired rank {args.rank} does not meet query ANI/AAI threshold {args.ani_threshold}")
            elif sg.fraction <= args.containment_threshold: # should this just be less than?
                status="below_threshold"
                notify(f"WARNING: classifying query {sg.query_name} at desired rank {args.rank} does not meet containment threshold {args.containment_threshold}")
            else:
                status="match"
            classif = ClassificationResult(sg.query_name, status, sg.rank, sg.fraction, sg.lineage, sg.query_md5, sg.query_filename, sg.f_weighted_at_rank, sg.bp_match_at_rank, sg.query_ani_at_rank)
            classifications[args.rank].append(classif)
            matched_queries.add(sg.query_name)
            if "krona" in args.output_format:
                lin_list = display_lineage(sg.lineage).split(';')
                krona_results.append((sg.fraction, *lin_list))
    else:
        # classify to the rank/match that passes the containment threshold.
        # To do - do we want to store anything for this match if nothing >= containment threshold?
        for rank in tax_utils.ascending_taxlist(include_strain=False):
            # gets best_at_rank for all queries in this gather_csv
            try:
                best_at_rank, seen_perfect, estimate_query_ani = tax_utils.summarize_gather_at(rank, tax_assign, gather_results, skip_idents=idents_missed,
                                                                                                keep_full_identifiers=args.keep_full_identifiers,
                                                                                                keep_identifier_versions = args.keep_identifier_versions,
                                                                                                best_only=True, seen_perfect=seen_perfect, estimate_query_ani=estimate_query_ani)
            except ValueError as exc:
                error(f"ERROR: {str(exc)}")
                sys.exit(-1)

            for sg in best_at_rank:
                status = 'nomatch'
                if sg.query_name in matched_queries:
                    continue
                if sg.query_ani_at_rank is not None and args.ani_threshold and sg.query_ani_at_rank >= args.ani_threshold:
                    status="match"
                elif sg.fraction >= args.containment_threshold:
                    status = "match"
                if status == "match":
                    classif = ClassificationResult(query_name=sg.query_name, status=status, rank=sg.rank,
                                                    fraction=sg.fraction, lineage=sg.lineage,
                                                    query_md5=sg.query_md5, query_filename=sg.query_filename,
                                                    f_weighted_at_rank=sg.f_weighted_at_rank, bp_match_at_rank=sg.bp_match_at_rank,
                                                    query_ani_at_rank= sg.query_ani_at_rank)
                    classifications[sg.rank].append(classif)
                    matched_queries.add(sg.query_name)
                    continue
                elif rank == "superkingdom" and status == "nomatch":
                    status="below_threshold"
                    classif = ClassificationResult(query_name=sg.query_name, status=status,
                                                   rank="", fraction=0, lineage="",
                                                   query_md5=sg.query_md5, query_filename=sg.query_filename,
                                                   f_weighted_at_rank=sg.f_weighted_at_rank, bp_match_at_rank=sg.bp_match_at_rank, 
                                                   query_ani_at_rank=sg.query_ani_at_rank)
                    classifications[sg.rank].append(classif)

    if not any([classifications, krona_results]):
        notify('No results for classification. Exiting.')
        sys.exit(-1)

    # write outputs
    if "csv_summary" in args.output_format:
        summary_outfile, limit_float = make_outfile(args.output_base, "classification", output_dir=args.output_dir)
        with FileOutputCSV(summary_outfile) as out_fp:
            tax_utils.write_classifications(classifications, out_fp, limit_float_decimals=limit_float)

    # write summarized output in human-readable format
    if "human" in args.output_format:
        summary_outfile, limit_float = make_outfile(args.output_base, "human", output_dir=args.output_dir)

        with FileOutput(summary_outfile) as out_fp:
            tax_utils.write_human_summary(classifications, out_fp, args.rank or "species")

    if "krona" in args.output_format:
        # classifications only at a single rank
        assert len(classifications) == 1
        if args.rank:
            assert args.rank in classifications

        krona_outfile, limit_float = make_outfile(args.output_base, "krona", output_dir=args.output_dir)
        with FileOutputCSV(krona_outfile) as out_fp:
            tax_utils.write_krona(args.rank, krona_results, out_fp)

    if "lineage_csv" in args.output_format:
        # should only classify at a single rank
        assert len(classifications) == 1
        if args.rank:
            assert args.rank in classifications

        lineage_outfile, _ = make_outfile(args.output_base, "lineage_csv",
                                          output_dir=args.output_dir)
        with FileOutputCSV(lineage_outfile) as out_fp:
            tax_utils.write_lineage_csv(classifications, out_fp)


def annotate(args):
    """
    Annotate gather results with taxonomic lineage for each match.

    Produces gather csv with lineage information as the final column.
    """

    set_quiet(args.quiet)

    # first, load taxonomic_assignments
    try:
        tax_assign = MultiLineageDB.load(args.taxonomy_csv,
                       keep_full_identifiers=args.keep_full_identifiers,
                       keep_identifier_versions=args.keep_identifier_versions,
                       force=args.force)
    except ValueError as exc:
        error(f"ERROR: {str(exc)}")
        sys.exit(-1)

    if not tax_assign:
        error(f'ERROR: No taxonomic assignments loaded from {",".join(args.taxonomy_csv)}. Exiting.')
        sys.exit(-1)

    # get gather_csvs from args
    gather_csvs = tax_utils.collect_gather_csvs(args.gather_csv, from_file=args.from_file)

    # handle each gather csv separately
    for n, g_csv in enumerate(gather_csvs):
        gather_results, idents_missed, total_missed, header = tax_utils.check_and_load_gather_csvs(g_csv, tax_assign, force=args.force,
                                                                                 fail_on_missing_taxonomy=args.fail_on_missing_taxonomy)

        if not gather_results:
            continue

        out_base = os.path.basename(g_csv.rsplit('.csv')[0])
        this_outfile, limit_float = make_outfile(out_base, "annotate", output_dir=args.output_dir)

        with FileOutputCSV(this_outfile) as out_fp:
            header.append("lineage")
            w = csv.DictWriter(out_fp, header, delimiter=',')
            w.writeheader()

            # add taxonomy info and then print directly
            for row in gather_results:
                match_ident = row['name']
                lineage = tax_utils.find_match_lineage(match_ident, tax_assign, skip_idents=idents_missed,
                                             keep_full_identifiers=args.keep_full_identifiers,
                                             keep_identifier_versions=args.keep_identifier_versions)
                row['lineage'] = display_lineage(lineage)
                w.writerow(row)


def prepare(args):
    "Combine multiple taxonomy databases into one and/or translate formats."
    notify("loading taxonomies...")
    try:
        tax_assign = MultiLineageDB.load(args.taxonomy_csv,
                                         force=args.force,
                       keep_full_identifiers=args.keep_full_identifiers,
                       keep_identifier_versions=args.keep_identifier_versions)
    except ValueError as exc:
        error("ERROR while loading taxonomies!")
        error(str(exc))
        sys.exit(-1)

    notify(f"...loaded {len(tax_assign)} entries.")

    notify(f"saving to '{args.output}', format {args.database_format}...")
    try:
        tax_assign.save(args.output, args.database_format)
    except ValueError as exc:
        error("ERROR while saving!")
        error(str(exc))
        sys.exit(-1)

    notify("done!")


def grep(args):
    term = args.pattern
    tax_assign = MultiLineageDB.load(args.taxonomy_csv,
                                     force=args.force)

    silent = args.silent or args.count

    notify(f"searching {len(args.taxonomy_csv)} taxonomy files for '{term}'")
    if args.invert_match:
        notify("-v/--invert-match specified; returning only lineages that do not match.")
    if args.rank:
        notify(f"limiting matches to {args.rank} level")

    # build the search pattern
    pattern = args.pattern
    if args.ignore_case:
        pattern = re.compile(pattern, re.IGNORECASE)
    else:
        pattern = re.compile(pattern)

    # determine if lineage matches.
    def find_pattern(lineage, select_rank):
        for (rank, name) in lineage:
            if select_rank is None or rank == select_rank:
                if pattern.search(name):
                    return True
        return False

    if args.invert_match:
        def search_pattern(l, r):
            return not find_pattern(l, r)
    else:
        search_pattern = find_pattern

    match_ident = []
    for ident, lineage in tax_assign.items():
        if search_pattern(lineage, args.rank):
            match_ident.append((ident, lineage))

    if silent:
        notify(f"found {len(match_ident)} matches.")
        notify("(no matches will be saved because of --silent/--count")
    else:
        with FileOutputCSV(args.output) as fp:
            w = csv.writer(fp)

            w.writerow(['ident'] + list(sourmash.lca.taxlist(include_strain=False)))
            for ident, lineage in sorted(match_ident):
                w.writerow([ident] + [ x.name for x in lineage ])

        notify(f"found {len(match_ident)} matches; saved identifiers to picklist file '{args.output}'")


def summarize(args):
    "Summarize multiple taxonomy databases."
    notify("loading taxonomies...")
    try:
        tax_assign = MultiLineageDB.load(args.taxonomy_files,
                                         force=args.force,
                       keep_full_identifiers=args.keep_full_identifiers,
                       keep_identifier_versions=args.keep_identifier_versions)
    except ValueError as exc:
        error("ERROR while loading taxonomies!")
        error(str(exc))
        sys.exit(-1)

    notify(f"...loaded {len(tax_assign)} entries.")

    print_results(f"number of distinct taxonomic lineages: {len(tax_assign)}")

    # count the number of distinct lineage names seen
    rank_counts = defaultdict(int)
    name_seen = set()
    for v in tax_assign.values():
        sofar = []
        for rank, name in v:
            if name not in name_seen:
                rank_counts[rank] += 1
                name_seen.add(name)

    rank_count_items = list(rank_counts.items())
    rank_count_items.sort(key=lambda x: x[1])
    for rank, count in rank_count_items:
        rank_name_str = f"{rank}:"
        print_results(f"rank {rank_name_str:<20s} {count} distinct taxonomic lineages")

    if args.output_lineage_information:
        notify("now calculating detailed lineage counts...")
        lineage_counts = Counter()
        for v in tax_assign.values():
            tup = v
            while tup:
                lineage_counts[tup] += 1
                tup = tup[:-1]
        notify("...done!")

        with FileOutputCSV(args.output_lineage_information) as fp:
            w = csv.writer(fp)
            w.writerow(['rank', 'lineage_count', 'lineage'])

            # output in order of most common
            for lineage, count in lineage_counts.most_common():
                rank = lineage[-1].rank
                lin = ";".join(zip_lineage(lineage, truncate_empty=True))
                w.writerow([rank, str(count), lin])

        n = len(lineage_counts)
        notify(f"saved {n} lineage counts to '{args.output_lineage_information}'")


def main(arglist=None):
    args = sourmash.cli.get_parser().parse_args(arglist)
    submod = getattr(sourmash.cli.sig, args.subcmd)
    mainmethod = getattr(submod, 'main')
    return mainmethod(args)


if __name__ == '__main__':
    main(sys.argv)
