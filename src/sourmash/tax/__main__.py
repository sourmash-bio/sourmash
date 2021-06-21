"""
Command-line entry point for 'python -m sourmash.tax'
"""
import sys
import csv
import os
from collections import defaultdict

import sourmash
from sourmash.lca.lca_utils import display_lineage

from ..sourmash_args import FileOutputCSV

from sourmash.logging import set_quiet, error, notify

from . import tax_utils
from .tax_utils import ClassificationResult

usage='''
sourmash taxonomy <command> [<args>] - manipulate/work with taxonomy information.
or
sourmash tax <command> [<args>]


** Commands can be:

summarize <gather_results> [<gather_results> ... ]        - summarize taxonomic information for metagenome gather results
classify <gather_results> [<gather_results> ... ]         - taxonomic classification of genomes from gather results
label <gather_results> [<gather_results> ... ]            - add taxonomic information to gather results csv(s)

** Use '-h' to get subcommand-specific help, e.g.

sourmash taxonomy summarize -h
'''

# some utils
def make_outfile(base, ext):
    if base == "-":
        return base
    return base + ext


##### taxonomy command line functions
def summarize(args):
    """
    summarize taxonomic information for metagenome gather results
    """
    set_quiet(args.quiet)

    # first, load taxonomic_assignments
    tax_assign = {}
    available_ranks = set()
    for tax_csv in args.taxonomy_csv:

        try:
            this_tax_assign, _, avail_ranks = tax_utils.load_taxonomy_csv(tax_csv, split_identifiers=not args.keep_full_identifiers,
                                              keep_identifier_versions = args.keep_identifier_versions,
                                              force=args.force)
            # maybe check for overlapping tax assignments? currently, later ones will override earlier ones
            tax_assign.update(this_tax_assign)
            available_ranks.update(set(avail_ranks))

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
        summarized_gather[rank], seen_perfect = tax_utils.summarize_gather_at(rank, tax_assign, gather_results, skip_idents=idents_missed,
                                                                split_identifiers=not args.keep_full_identifiers,
                                                                keep_identifier_versions = args.keep_identifier_versions,
                                                                seen_perfect = seen_perfect)

    # write summarized output csv
    if "summary" in args.output_format:
        summary_outfile = make_outfile(args.output_base, ".summarized.csv")
        with FileOutputCSV(summary_outfile) as out_fp:
            tax_utils.write_summary(summarized_gather, out_fp)

    # if lineage summary table
    if "lineage_summary" in args.output_format:
        lineage_outfile = make_outfile(args.output_base, ".lineage_summary.tsv")

        ## aggregate by lineage, by query
        lineageD, query_names, num_queries = tax_utils.aggregate_by_lineage_at_rank(summarized_gather[args.rank], by_query=True)

        with FileOutputCSV(lineage_outfile) as out_fp:
            tax_utils.write_lineage_sample_frac(query_names, lineageD, out_fp, format_lineage=True, sep='\t')

    # write summarized --> krona output tsv
    if "krona" in args.output_format:
        krona_resultslist = tax_utils.format_for_krona(args.rank, summarized_gather)

        krona_outfile = make_outfile(args.output_base, ".krona.tsv")
        with FileOutputCSV(krona_outfile) as out_fp:
            tax_utils.write_krona(args.rank, krona_resultslist, out_fp)


def classify(args):
    """
    taxonomic classification of genomes from gather results
    """
    set_quiet(args.quiet)

    # first, load taxonomic_assignments
    tax_assign = {}
    available_ranks = set()
    for tax_csv in args.taxonomy_csv:

        try:
            this_tax_assign, _, avail_ranks = tax_utils.load_taxonomy_csv(tax_csv, split_identifiers=not args.keep_full_identifiers,
                                              keep_identifier_versions = args.keep_identifier_versions,
                                              force=args.force)
            # maybe check for overlapping tax assignments? currently later ones will override earlier ones
            tax_assign.update(this_tax_assign)
            available_ranks.update(set(avail_ranks))
        except ValueError as exc:
            error(f"ERROR: {str(exc)}")

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
    if args.rank:
        best_at_rank, seen_perfect = tax_utils.summarize_gather_at(args.rank, tax_assign, gather_results, skip_idents=idents_missed,
                                                     split_identifiers=not args.keep_full_identifiers,
                                                     keep_identifier_versions = args.keep_identifier_versions,
                                                     best_only=True, seen_perfect=seen_perfect)

       # best at rank is a list of SummarizedGather tuples
        for sg in best_at_rank:
            status = 'nomatch'
            if sg.query_name in matched_queries:
                continue
            if sg.fraction <= args.containment_threshold:
                status="below_threshold"
                notify(f"WARNING: classifying query {sg.query_name} at desired rank {args.rank} does not meet containment threshold {args.containment_threshold}")
            else:
                status="match"
            classif = ClassificationResult(sg.query_name, status, sg.rank, sg.fraction, sg.lineage)
            classifications[args.rank].append(classif)
            matched_queries.add(sg.query_name)
            if "krona" in args.output_format:
                lin_list = display_lineage(sg.lineage).split(';')
                krona_results.append((sg.query_name, sg.fraction, *lin_list))
    else:
        # classify to the match that passes the containment threshold.
        # To do - do we want to store anything for this match if nothing >= containment threshold?
        for rank in tax_utils.ascending_taxlist(include_strain=False):
            # gets best_at_rank for all queries in this gather_csv
            best_at_rank, seen_perfect = tax_utils.summarize_gather_at(rank, tax_assign, gather_results, skip_idents=idents_missed,
                                                         split_identifiers=not args.keep_full_identifiers,
                                                         keep_identifier_versions = args.keep_identifier_versions,
                                                         best_only=True, seen_perfect=seen_perfect)

            for sg in best_at_rank:
                status = 'nomatch'
                if sg.query_name in matched_queries:
                    continue
                if sg.fraction >= args.containment_threshold:
                    status = "match"
                    classif = ClassificationResult(sg.query_name, status, sg.rank, sg.fraction, sg.lineage)
                    classifications[args.rank].append(classif)
                    matched_queries.add(sg.query_name)
                    if "krona" in args.output_format:
                        lin_list = display_lineage(sg.lineage).split(';')
                        krona_results.append((sg.query_name, sg.fraction, *lin_list))
                    break
                if rank == "superkingdom" and status == "nomatch":
                    status="below_threshold"
                    classif = ClassificationResult(sg.query_name, status, "", 0, "")
                    classifications[args.rank].append(classif)

    if not any([classifications, krona_results]):
        notify('No results for classification. Exiting.')
        sys.exit(-1)

    # write outputs
    if "summary" in args.output_format:
        summary_outfile = make_outfile(args.output_base, ".classifications.csv")
        with FileOutputCSV(summary_outfile) as out_fp:
            tax_utils.write_classifications(classifications, out_fp)

    if "krona" in args.output_format:
        krona_outfile = make_outfile(args.output_base, ".krona.tsv")
        with FileOutputCSV(krona_outfile) as out_fp:
            tax_utils.write_krona(args.rank, krona_results, out_fp)


def annotate(args):
    """
    Annotate gather results with taxonomic lineage for each match.

    Produces gather csv with lineage information as the final column.
    """

    set_quiet(args.quiet)

    # first, load taxonomic_assignments
    tax_assign = {}
    this_tax_assign = None
    for tax_csv in args.taxonomy_csv:

        try:
            this_tax_assign, _, avail_ranks = tax_utils.load_taxonomy_csv(tax_csv, split_identifiers=not args.keep_full_identifiers,
                                              keep_identifier_versions = args.keep_identifier_versions,
                                              force=args.force)
        except ValueError as exc:
            error(f"ERROR: {str(exc)}")
            sys.exit(-1)

        # maybe check for overlapping tax assignments? currently later ones will override earlier ones
        if this_tax_assign:
            tax_assign.update(this_tax_assign)

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
        out_path = os.path.join(args.output_dir, out_base)
        this_outfile = make_outfile(out_path, ".with-lineages.csv")

        with FileOutputCSV(this_outfile) as out_fp:
            header.append("lineage")
            w = csv.DictWriter(out_fp, header, delimiter=',')
            w.writeheader()

            # add taxonomy info and then print directly
            for row in gather_results:
                match_ident = row['name']
                lineage = tax_utils.find_match_lineage(match_ident, tax_assign, skip_idents=idents_missed,
                                             split_identifiers=not args.keep_full_identifiers,
                                             keep_identifier_versions=args.keep_identifier_versions)
                row['lineage'] = display_lineage(lineage)
                w.writerow(row)


def main(arglist=None):
    args = sourmash.cli.get_parser().parse_args(arglist)
    submod = getattr(sourmash.cli.sig, args.subcmd)
    mainmethod = getattr(submod, 'main')
    return mainmethod(args)


if __name__ == '__main__':
    main(sys.argv)
