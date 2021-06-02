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

    # load gather results and taxonomy assignments
    gather_results = tax_utils.load_gather_results(args.gather_results)
    if not gather_results:
        notify(f'No gather results loaded from {args.gather_results}. Exiting.')
        sys.exit(-1)

    tax_assign, _ = load_taxonomy_assignments(args.taxonomy_csv, use_headers=True,
                                              split_identifiers=not args.keep_full_identifiers,
                                              keep_identifier_versions = args.keep_identifier_versions,
                                              force=args.force)
    if not tax_assign:
        notify(f'No taxonomic assignments loaded from {args.taxonomy_csv}. Exiting.')
        sys.exit(-1)

    # check for match identites not found in lineage spreadsheets
    n_missed, ident_missed = tax_utils.find_missing_identities(gather_results, tax_assign)
    if n_missed:
        notify(f'The following are missing from the taxonomy information: {",".join(ident_missed)}')
        if args.fail_on_missing_taxonomy:
            notify(f'Failing on missing taxonomy, as requested via --fail-on-missing-taxonomy.')
            sys.exit(-1)

    # actually summarize at rank
    summarized_gather = {}
    for rank in sourmash.lca.taxlist(include_strain=False):
        summarized_gather[rank] = tax_utils.summarize_gather_at(rank, tax_assign, gather_results, skip_idents=ident_missed,
                                                                split_identifiers=not args.keep_full_identifiers,
                                                                keep_identifier_versions = args.keep_identifier_versions)

    # write summarized output csv
    if "summary" in args.output_format:
        summary_outfile = make_outfile(args.output_base, ".summarized.csv")
        with FileOutputCSV(summary_outfile) as csv_fp:
            tax_utils.write_summary(summarized_gather, csv_fp)

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

    # load gather results for each genome and summarize with --best-only to classify
    gather_info, cli_gather_res, csv_gather_res = [],[],[]
    query_name = None
    if args.gather_results:
        query_name = args.query_name
        cli_gather_res = [(query_name, args.gather_results)]
    if args.from_csv:
        csv_gather_res, seen_idents = tax_utils.load_gather_files_from_csv(args.from_csv)
        if query_name and query_name in seen_idents:
            notify("query name is also found in --from-csv filelist!")
            if args.force:
                fixed_csv_res = []
                #remove query_name result line from csv_gather_res -- is this a good desired behavior?
                notify(f"--force is set. Removing {query_name} entry from the --from-csv gather results in favor of cli input.")
                for (ident, gather_res) in csv_gather_res:
                    if ident != query_name:
                        fixed_csv_res.append((ident, gather_res))
                csv_gather_res = fixed_csv_res
            else:
                notify('Exiting.')
                sys.exit(-1)

    # full list of (ident,gather_results)
    gather_info = cli_gather_res + csv_gather_res

    classifications = defaultdict(list)
    krona_results = []
    num_empty=0
    for n, (name, g_results) in enumerate(gather_info):

        gather_results = tax_utils.load_gather_results(g_results)
        if not gather_results:
            notify(f'No gather results loaded from {args.gather_results}.')
            num_empty+=1
            if args.force:
                notify('--force is set. Attempting to continue to next set of gather results.')
                continue
            else:
                notify('Exiting.')
                sys.exit(-1)

        # check for match identites not found in lineage spreadsheets
        n_missed, ident_missed = tax_utils.find_missing_identities(gather_results, tax_assign)
        if n_missed:
            notify(f'The following are missing from the taxonomy information: {",".join(ident_missed)}')
            if args.fail_on_missing_taxonomy:
                notify(f'Failing on missing taxonomy, as requested via --fail-on-missing-taxonomy.')
                sys.exit(-1)

        # if --rank is specified, classify to that rank
        # to do, what to do if don't have gather results at desired rank (e.g. strain)?
        if args.rank:
            # todo: check we have gather results at this rank
            #if not tax_utils.check_taxonomy_exists(tax_assign, args.rank):
            #    notify(f"No taxonomic information at rank {args.rank}: cannot classify at this rank")
            best_at_rank = tax_utils.summarize_gather_at(args.rank, tax_assign, gather_results, skip_idents=ident_missed,
                                                         split_identifiers=not args.keep_full_identifiers,
                                                         keep_identifier_versions = args.keep_identifier_versions,
                                                         best_only=True)[0]
            (lineage,containment) = best_at_rank
            if containment <= args.containment_threshold:
                notify(f"WARNING: classifying at desired rank {args.rank} does not meet containment threshold {args.containment_threshold}")
            classifications[args.rank].append((name, best_at_rank))
            if "krona" in args.output_format:
                lin_list = display_lineage(lineage).split(';')
                krona_results.append((containment, *lin_list))
        else:
            # classify to the match that passes the containment threshold. To do - do we want to report anything if nothing >= containment threshold?
            for rank in tax_utils.ascending_taxlist(include_strain=False):
                best_at_rank = tax_utils.summarize_gather_at(rank, tax_assign, gather_results, skip_idents=ident_missed,
                                                             split_identifiers=not args.keep_full_identifiers,
                                                             keep_identifier_versions = args.keep_identifier_versions,
                                                             best_only=True)[0]
                (lineage,containment) = best_at_rank
                if containment >= args.containment_threshold:
                    classifications[rank].append((name, best_at_rank))
                    if "krona" in args.output_format:
                        lin_list = display_lineage(lineage).split(';')
                        krona_results.append((containment, *lin_list))
                    break

    notify(f'loaded {n+1-num_empty} gather files for classification.')

    if not any([classifications,krona_results]):
        notify(f'No results for classification. Exiting.')
        sys.exit(-1)

    # write output csv
    if "summary" in args.output_format:
        summary_outfile = make_outfile(args.output_base, ".classifications.csv")
        with FileOutputCSV(summary_outfile) as csv_fp:
            tax_utils.write_classifications(classifications, csv_fp)

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
    linD, all_samples = tax_utils.combine_sumgather_csvs_by_lineage(args.summarized_gather_results, rank=args.rank)
    #if not linD:
    #    notify(f'No summarized gather results loaded from {args.summarized_gather_results}. Exiting.')
    #    sys.exit(-1)

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
