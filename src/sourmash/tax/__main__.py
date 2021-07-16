"""
Command-line entry point for 'python -m sourmash.tax'
"""
import sys
import csv
import os
from collections import defaultdict

import sourmash
from ..sourmash_args import FileOutputCSV
from sourmash.logging import set_quiet, error, notify
from sourmash.lca.lca_utils import display_lineage

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

# some utils
def make_outfile(base, output_type, *, output_dir = ""):
    limit_float_decimals=False
    if base == "-":
        limit_float_decimals=True
        return base, limit_float_decimals
    ext=""
    if output_type == 'csv_summary':
        ext = '.summarized.csv'
    elif output_type == 'classification':
        ext = '.classifications.csv'
    elif output_type == 'krona':
        ext = '.krona.tsv'
    elif output_type == 'lineage_summary':
        ext = '.lineage_summary.tsv'
    elif output_type == 'annotate':
        ext = '.with-lineages.csv'
    fname = base+ext
    if output_dir:
        fname = os.path.join(output_dir, fname)
    notify(f"saving `{output_type}` output to {fname}.")
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
            summarized_gather[rank], seen_perfect = tax_utils.summarize_gather_at(rank, tax_assign, gather_results, skip_idents=idents_missed,
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
    if args.rank:
        try:
            best_at_rank, seen_perfect = tax_utils.summarize_gather_at(args.rank, tax_assign, gather_results, skip_idents=idents_missed,
                                                     keep_full_identifiers=args.keep_full_identifiers,
                                                     keep_identifier_versions = args.keep_identifier_versions,
                                                     best_only=True, seen_perfect=seen_perfect)
        except ValueError as exc:
            error(f"ERROR: {str(exc)}")
            sys.exit(-1)

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
            classif = ClassificationResult(sg.query_name, status, sg.rank, sg.fraction, sg.lineage, sg.query_md5, sg.query_filename, sg.f_weighted_at_rank, sg.bp_match_at_rank)
            classifications[args.rank].append(classif)
            matched_queries.add(sg.query_name)
            if "krona" in args.output_format:
                lin_list = display_lineage(sg.lineage).split(';')
                krona_results.append((sg.fraction, *lin_list))
    else:
        # classify to the match that passes the containment threshold.
        # To do - do we want to store anything for this match if nothing >= containment threshold?
        for rank in tax_utils.ascending_taxlist(include_strain=False):
            # gets best_at_rank for all queries in this gather_csv
            try:
                best_at_rank, seen_perfect = tax_utils.summarize_gather_at(rank, tax_assign, gather_results, skip_idents=idents_missed,
                                                         keep_full_identifiers=args.keep_full_identifiers,
                                                         keep_identifier_versions = args.keep_identifier_versions,
                                                         best_only=True, seen_perfect=seen_perfect)
            except ValueError as exc:
                error(f"ERROR: {str(exc)}")
                sys.exit(-1)

            for sg in best_at_rank:
                status = 'nomatch'
                if sg.query_name in matched_queries:
                    continue
                if sg.fraction >= args.containment_threshold:
                    status = "match"
                    classif = ClassificationResult(sg.query_name, status, sg.rank, sg.fraction, sg.lineage, sg.query_md5, sg.query_filename, sg.f_weighted_at_rank, sg.bp_match_at_rank)
                    classifications[sg.rank].append(classif)
                    matched_queries.add(sg.query_name)
                    continue
                if rank == "superkingdom" and status == "nomatch":
                    status="below_threshold"
                    classif = ClassificationResult(query_name=sg.query_name, status=status,
                                                   rank="", fraction=0, lineage="",
                                                   query_md5=sg.query_md5, query_filename=sg.query_filename,
                                                   f_weighted_at_rank=sg.f_weighted_at_rank, bp_match_at_rank=sg.bp_match_at_rank)
                    classifications[sg.rank].append(classif)

    if not any([classifications, krona_results]):
        notify('No results for classification. Exiting.')
        sys.exit(-1)

    # write outputs
    if "csv_summary" in args.output_format:
        summary_outfile, limit_float = make_outfile(args.output_base, "classification", output_dir=args.output_dir)
        with FileOutputCSV(summary_outfile) as out_fp:
            tax_utils.write_classifications(classifications, out_fp, limit_float_decimals=limit_float)

    if "krona" in args.output_format:
        krona_outfile, limit_float = make_outfile(args.output_base, "krona", output_dir=args.output_dir)
        with FileOutputCSV(krona_outfile) as out_fp:
            tax_utils.write_krona(args.rank, krona_results, out_fp)


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


def main(arglist=None):
    args = sourmash.cli.get_parser().parse_args(arglist)
    submod = getattr(sourmash.cli.sig, args.subcmd)
    mainmethod = getattr(submod, 'main')
    return mainmethod(args)


if __name__ == '__main__':
    main(sys.argv)
