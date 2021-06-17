"""
Command-line entry point for 'python -m sourmash.tax'
"""
import sys
import csv
import json
import os
from collections import defaultdict

import sourmash
import copy
from sourmash.sourmash_args import FileOutput
from sourmash.lca.lca_utils import pop_to_rank, display_lineage

from ..sourmash_args import FileOutputCSV

from sourmash.logging import set_quiet, error, notify, set_quiet, print_results, debug
from sourmash import sourmash_args

from . import tax_utils

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
    for tax_csv in args.taxonomy_csv:

        this_tax_assign, _ = tax_utils.load_taxonomy_csv(tax_csv, split_identifiers=not args.keep_full_identifiers,
                                              keep_identifier_versions = args.keep_identifier_versions,
                                              force=args.force)
        # to do -- maybe check for overlapping tax assignments? rn later ones will override earlier ones
        tax_assign.update(this_tax_assign)

    if not tax_assign:
        notify(f'No taxonomic assignments loaded from {args.taxonomy_csv}. Exiting.')
        sys.exit(-1)

    # next, collect and load gather results
    gather_csvs = tax_utils.collect_gather_csvs(args.gather_results, from_file= args.from_file)
    gather_results, idents_missed, total_missed, _ = tax_utils.check_and_load_gather_csvs(gather_csvs, tax_assign, force=args.force,
                                                                                       fail_on_missing_taxonomy=args.fail_on_missing_taxonomy)

    if not gather_results:
        notify(f'No gather results loaded. Exiting.')
        sys.exit(-1)

    # actually summarize at rank
    summarized_gather = {}
    for rank in sourmash.lca.taxlist(include_strain=False):
        summarized_gather[rank] = tax_utils.summarize_gather_at(rank, tax_assign, gather_results, skip_idents=idents_missed,
                                                                split_identifiers=not args.keep_full_identifiers,
                                                                keep_identifier_versions = args.keep_identifier_versions)

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
    for tax_csv in args.taxonomy_csv:

        this_tax_assign, _ = tax_utils.load_taxonomy_csv(tax_csv, split_identifiers=not args.keep_full_identifiers,
                                              keep_identifier_versions = args.keep_identifier_versions,
                                              force=args.force)
        # to do -- maybe check for overlapping tax assignments? rn later ones will override earlier ones
        tax_assign.update(this_tax_assign)

    if not tax_assign:
        notify(f'No taxonomic assignments loaded from {args.taxonomy_csv}. Exiting.')
        sys.exit(-1)

    # get gather_csvs from args
    gather_csvs = tax_utils.collect_gather_csvs(args.gather_results, from_file=args.from_file)

    classifications = defaultdict(list)
    matched_queries=set()
    krona_results = []
    num_empty=0
    status = "nomatch"

    # handle each gather result separately
    for n, g_csv in enumerate(gather_csvs):
        gather_results, idents_missed, total_missed, _ = tax_utils.check_and_load_gather_csvs(g_csv, tax_assign, force=args.force,
                                                                                 fail_on_missing_taxonomy=args.fail_on_missing_taxonomy)

        if not gather_results:
            continue

        # if --rank is specified, classify to that rank
        if args.rank:
            best_at_rank = tax_utils.summarize_gather_at(args.rank, tax_assign, gather_results, skip_idents=idents_missed,
                                                         split_identifiers=not args.keep_full_identifiers,
                                                         keep_identifier_versions = args.keep_identifier_versions,
                                                         best_only=True)

           # this now returns list of SummarizedGather tuples
            for (query_name, rank, fraction, lineage) in best_at_rank:
                status = 'nomatch'
                if query_name in matched_queries:
                    notify(f"already matched query {query_name}. Skipping...")
                    continue
                if fraction <= args.containment_threshold:
                    status="below_threshold"
                    notify(f"WARNING: classifying at desired rank {args.rank} does not meet containment threshold {args.containment_threshold}")
                else:
                    status="match"
                classifications[args.rank].append((query_name, status, rank, fraction, lineage))
                matched_queries.add(query_name)
                if "krona" in args.output_format:
                    lin_list = display_lineage(lineage).split(';')
                    krona_results.append((containment, *lin_list))
        else:
            # classify to the match that passes the containment threshold.
            # To do - do we want to store anything for this match if nothing >= containment threshold?
            for rank in tax_utils.ascending_taxlist(include_strain=False):
                # gets best_at_rank for all queries in this gather_csv
                best_at_rank = tax_utils.summarize_gather_at(rank, tax_assign, gather_results, skip_idents=idents_missed,
                                                             split_identifiers=not args.keep_full_identifiers,
                                                             keep_identifier_versions = args.keep_identifier_versions,
                                                             best_only=True)

                for (query_name, rank, fraction, lineage) in best_at_rank:
                    status = 'nomatch'
                    if query_name in matched_queries:
                        notify(f"already matched query {query_name}. Skipping...")
                        continue
                    if fraction >= args.containment_threshold:
                        status = "match"
                        classifications[args.rank].append((query_name, status, rank, fraction, lineage))
                        matched_queries.add(query_name)
                        if "krona" in args.output_format:
                            lin_list = display_lineage(lineage).split(';')
                            krona_results.append((query_name, containment, *lin_list))
                        break
                    if rank == "superkingdom" and status == "nomatch":
                        status="below_threshold"
                        classifications[args.rank].append((query_name, status, "", 0, ""))

    notify(f'loaded {n} gather files for classification.')

    if not any([classifications, krona_results]):
        notify(f'No results for classification. Exiting.')
        sys.exit(-1)

    # write outputs
    if "summary" in args.output_format:
        summary_outfile = make_outfile(args.output_base, ".classifications.csv")
        with FileOutputCSV(summary_outfile) as out_fp:
            #tax_utils.write_summary(classifications, out_fp)
            tax_utils.write_classifications(classifications, out_fp)

    if "krona" in args.output_format:
        krona_outfile = make_outfile(args.output_base, ".krona.tsv")
        with FileOutputCSV(krona_outfile) as out_fp:
            tax_utils.write_krona(args.rank, krona_results, out_fp)


def label(args):
    """
    Integrate lineage information into gather results.

    Produces gather csv with lineage information as the final column.
    """

    set_quiet(args.quiet)

    # first, load taxonomic_assignments
    tax_assign = {}
    for tax_csv in args.taxonomy_csv:

        this_tax_assign, _ = tax_utils.load_taxonomy_csv(tax_csv, split_identifiers=not args.keep_full_identifiers,
                                              keep_identifier_versions = args.keep_identifier_versions,
                                              force=args.force)

        # to do -- maybe check for overlapping tax assignments? rn later ones will override earlier ones
        tax_assign.update(this_tax_assign)

    if not tax_assign:
        notify(f'No taxonomic assignments loaded from {args.taxonomy_csv}. Exiting.')
        sys.exit(-1)

    # get gather_csvs from args
    gather_csvs = tax_utils.collect_gather_csvs(args.gather_results, from_file=args.from_file)

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
