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
from sourmash.lca.command_index import load_taxonomy_assignments

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
combine <summarized_gather_results> [<summarized_gather_results> ... ] - combine outputs of `summarize` for multiple samples
classify <gather_results> [<gather_results> ... ]   - taxonomic classification of genomes from gather results

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
    tax_assign, _ = load_taxonomy_assignments(args.taxonomy_csv, use_headers=True,
                                              split_identifiers=not args.keep_full_identifiers,
                                              keep_identifier_versions = args.keep_identifier_versions,
                                              force=args.force)

    if not tax_assign:
        notify(f'No taxonomic assignments loaded from {args.taxonomy_csv}. Exiting.')
        sys.exit(-1)

    # next, collect and load gather results
    gather_csvs = tax_utils.collect_gather_csvs(args.gather_results, from_file= args.from_file)
    gather_results, idents_missed, total_missed = tax_utils.check_and_load_gather_csvs(gather_csvs, tax_assign, force=args.force,
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
        with FileOutputCSV(summary_outfile) as csv_fp:
            tax_utils.write_summary(summarized_gather, csv_fp)

    # if lineage summary table
    if "lineage_summary" in args.output_format:
        lineage_outfile = make_outfile(args.output_base, ".lineage_summary.tsv")

        ## aggregate by lineage, by query
        lineageD, query_names, num_queries = tax_utils.aggregate_by_lineage_at_rank(summarized_gather[args.rank], by_query=True)

        with FileOutputCSV(lineage_outfile) as csv_fp:
            tax_utils.write_lineage_sample_frac(query_names, lineageD, csv_fp, flatten_lineage=True, sep='\t')

    # write summarized --> krona output csv
    if "krona" in args.output_format:
        krona_resultslist = tax_utils.format_for_krona(args.rank, summarized_gather)

        krona_outfile = make_outfile(args.output_base, ".krona.tsv")
        with FileOutputCSV(krona_outfile) as out_fp:
            tax_utils.write_krona(args.rank, krona_resultslist, out_fp)


def classify(args):
    """
    taxonomic classification of genomes from gather results
    """
    # classify:: summarize at rank, choose best match
    ## currently reports a single rank. do we want to optionally report at all ranks? (no, bc summarize does that?)
    set_quiet(args.quiet)

    # load taxonomy assignments
    tax_assign, _ = load_taxonomy_assignments(args.taxonomy_csv, use_headers=True,
                                              split_identifiers=not args.keep_full_identifiers,
                                              keep_identifier_versions = args.keep_identifier_versions,
                                              force=args.force)

    if not tax_assign:
        notify(f'No taxonomic assignments loaded from {args.taxonomy_csv}. Exiting.')
        sys.exit(-1)

    # get gather_csvs from args
    gather_csvs = tax_utils.collect_gather_csvs(args.gather_results, from_file=args.from_file)

    classifications = defaultdict(list)
    seen_queries=set()
    krona_results = []
    num_empty=0

    # handle each gather result separately
    for n, g_csv in enumerate(gather_csvs):
        gather_results, idents_missed, total_missed = tax_utils.check_and_load_gather_csvs(g_csv, tax_assign, force=args.force,
                                                                                 fail_on_missing_taxonomy=args.fail_on_missing_taxonomy)

        if not gather_results:
            continue

        # if --rank is specified, classify to that rank
        # to do, what to do if don't have gather results at desired rank (e.g. strain)?
        if args.rank:
            # todo: check we have gather results at this rank
            # better idea: return available taxonomic ranks from tax_assign! then check that rank is in these.
            #if not tax_utils.check_taxonomy_exists(tax_assign, args.rank):
            #    notify(f"No taxonomic information at rank {args.rank}: cannot classify at this rank")
            best_at_rank = tax_utils.summarize_gather_at(args.rank, tax_assign, gather_results, skip_idents=idents_missed,
                                                         split_identifiers=not args.keep_full_identifiers,
                                                         keep_identifier_versions = args.keep_identifier_versions,
                                                         best_only=True)
           # this now returns list of SummarizedGather tuples
            for (query_name, rank, fraction, lineage) in best_at_rank:
                if query_name in seen_queries:
                    notify(f"WARNING: duplicate query {query_name}. Skipping...")
                    continue
                if fraction <= args.containment_threshold:
                    notify(f"WARNING: classifying at desired rank {args.rank} does not meet containment threshold {args.containment_threshold}")
                classifications[args.rank].append((query_name, rank, fraction, lineage))
                seen_queries.add(query_name)
                if "krona" in args.output_format:
                    lin_list = display_lineage(lineage).split(';')
                    krona_results.append((containment, *lin_list))
        else:
            # classify to the match that passes the containment threshold.
            # To do - do we want to report anything if nothing >= containment threshold?
            for rank in tax_utils.ascending_taxlist(include_strain=False):
                # gets for all queries at once
                best_at_rank = tax_utils.summarize_gather_at(rank, tax_assign, gather_results, skip_idents=idents_missed,
                                                             split_identifiers=not args.keep_full_identifiers,
                                                             keep_identifier_versions = args.keep_identifier_versions,
                                                             best_only=True)

                for (query_name, rank, fraction, lineage) in best_at_rank:
                    if query_name in seen_queries:
                        notify(f"WARNING: duplicate query {query_name}. Skipping...")
                        continue
                    if fraction >= args.containment_threshold:
                        classifications[args.rank].append((query_name, rank, fraction, lineage))
                        seen_queries.add(query_name)
                        if "krona" in args.output_format:
                            lin_list = display_lineage(lineage).split(';')
                            krona_results.append((query_name, containment, *lin_list))
                        break

    notify(f'loaded {n} gather files for classification.')

    if not any([classifications, krona_results]):
        notify(f'No results for classification. Exiting.')
        sys.exit(-1)

    # write output csv
    if "summary" in args.output_format:
        summary_outfile = make_outfile(args.output_base, ".classifications.csv")
        with FileOutputCSV(summary_outfile) as csv_fp:
            tax_utils.write_summary(classifications, csv_fp)

    if "krona" in args.output_format:
        krona_outfile = make_outfile(args.output_base, ".krona.tsv")
        with FileOutputCSV(krona_outfile) as csv_fp:
            tax_utils.write_krona(args.rank, krona_results, csv_fp)


def combine(args):
    """
    Combine summarize gather results by lineage and sample.

    Takes in one or more output csvs from `sourmash taxonomy summarize`
    and produces a tab-separated file with fractions for each sample.

    Uses the file basename (minus .csv extension) as sample identifier.

    example output:

    lineage    sample1  sample2 sample3
    lin_a     0.4    0.17     0.6
    lin_b     0.0    0.0      0.1
    lin_c     0.3    0.4      0.2

    """

    set_quiet(args.quiet)

    # load summarized gather csvs into lineage dictionary
    sumgather_csvs = tax_utils.collect_gather_csvs(args.summarized_gather_results, from_file=args.from_file)

    linD, all_samples = tax_utils.combine_sumgather_csvs_by_lineage(sumgather_csvs, rank=args.rank, force=args.force)
    if not linD:
        notify(f'No summarized gather results loaded from {args.summarized_gather_results}. Exiting.')
        sys.exit(-1)

    # write output csv
    if "csv" in args.output_format:
        outfile = make_outfile(args.output_base, ".combined.csv")
        with FileOutputCSV(outfile) as csv_fp:
            tax_utils.write_lineage_sample_frac(all_samples, linD, csv_fp, sep=",")
    if "tsv" in args.output_format:
        outfile = make_outfile(args.output_base, ".combined.tsv")
        with FileOutputCSV(outfile) as csv_fp:
            tax_utils.write_lineage_sample_frac(all_samples, linD, csv_fp, sep="\t")


def main(arglist=None):
    args = sourmash.cli.get_parser().parse_args(arglist)
    submod = getattr(sourmash.cli.sig, args.subcmd)
    mainmethod = getattr(submod, 'main')
    return mainmethod(args)


if __name__ == '__main__':
    main(sys.argv)
