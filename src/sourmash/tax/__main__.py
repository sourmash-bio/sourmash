"""
Command-line entry point for 'python -m sourmash.tax'
"""
import sys
import csv
import os
from collections import defaultdict, Counter
from dataclasses import asdict, fields
import re

import sourmash
from ..sourmash_args import FileOutputCSV, FileOutput
from sourmash.logging import set_quiet, error, notify, print_results
from sourmash.lca.lca_utils import zip_lineage

from . import tax_utils
from .tax_utils import MultiLineageDB, GatherRow

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

# outfile utils
_output_type_to_ext = {
    'csv_summary': '.summarized.csv',
    'classification': '.classifications.csv',
    'krona': '.krona.tsv',
    'lineage_summary': '.lineage_summary.tsv',
    'annotate': '.with-lineages.csv',
    'human': '.human.txt',
    'lineage_csv': '.lineage.csv',
    'kreport': ".kreport.txt",
    'lingroup_report': ".lingroup_report.tsv"
    }

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
                       force=args.force, LIN_taxonomy=args.LIN_taxonomy)
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
         query_gather_results = tax_utils.check_and_load_gather_csvs(gather_csvs, tax_assign, force=args.force,
                                                                     fail_on_missing_taxonomy=args.fail_on_missing_taxonomy,
                                                                     keep_full_identifiers=args.keep_full_identifiers,
                                                                     keep_identifier_versions = args.keep_identifier_versions,
                                                                     LIN_taxonomy=args.LIN_taxonomy,
                                                                     )
    except ValueError as exc:
        error(f"ERROR: {str(exc)}")
        sys.exit(-1)

    if not query_gather_results:
        notify('No gather results loaded. Exiting.')
        sys.exit(-1)

    single_query_output_formats =  ['csv_summary', 'kreport', "lineage_summary"]
    desired_single_outputs = []
    if len(query_gather_results) > 1: # working with multiple queries
        desired_single_outputs = [x for x in args.output_format if x in single_query_output_formats]
        if desired_single_outputs:
            notify(f"WARNING: found results for multiple gather queries. Can only output multi-query result formats: skipping {', '.join(desired_single_outputs)}")
        # remove single query outputs from output format
        args.output_format = [x for x in args.output_format if x not in single_query_output_formats]
        if not args.output_format: # or do we want to insert `human` here so we always report something?
            error(f"ERROR: No output formats remaining.")
            sys.exit(-1)

    # for each queryResult, actually summarize at rank, reporting any errors that occur.
    for queryResult in query_gather_results:
        try:
            queryResult.build_summarized_result()
        except ValueError as exc:
            error(f"ERROR: {str(exc)}")
            sys.exit(-1)

    # write summarized output in human-readable format
    if "lineage_summary" in args.output_format:
        lineage_outfile, limit_float = make_outfile(args.output_base, "lineage_summary", output_dir=args.output_dir)

        ## aggregate by lineage by query
        lineageD, query_names= tax_utils.aggregate_by_lineage_at_rank(query_gather_results=query_gather_results,
                                                                      rank=args.rank, by_query=True)

        with FileOutputCSV(lineage_outfile) as out_fp:
            tax_utils.write_lineage_sample_frac(query_names, lineageD, out_fp, sep='\t')

    # write summarized --> krona output tsv
    if "krona" in args.output_format:
        krona_results, header =  tax_utils.format_for_krona(query_gather_results, rank=args.rank)

        krona_outfile, limit_float = make_outfile(args.output_base, "krona", output_dir=args.output_dir)
        with FileOutputCSV(krona_outfile) as out_fp:
            tax_utils.write_krona(header, krona_results, out_fp)

    if "human" in args.output_format:
        summary_outfile, limit_float = make_outfile(args.output_base, "human", output_dir=args.output_dir)

        with FileOutput(summary_outfile) as out_fp:
            tax_utils.write_human_summary(query_gather_results, out_fp, args.rank or "species")

    # write summarized output csv
    single_query_results = query_gather_results[0]
    if "csv_summary" in args.output_format:
        summary_outfile, limit_float = make_outfile(args.output_base, "csv_summary", output_dir=args.output_dir)
        with FileOutputCSV(summary_outfile) as out_fp:
            tax_utils.write_summary(query_gather_results, out_fp, limit_float_decimals=limit_float)

    # write summarized --> kreport output tsv
    if "kreport" in args.output_format:
        kreport_outfile, limit_float = make_outfile(args.output_base, "kreport", output_dir=args.output_dir)

        with FileOutputCSV(kreport_outfile) as out_fp:
            header, kreport_results = single_query_results.make_kreport_results()
            tax_utils.write_output(header, kreport_results, out_fp, sep="\t", write_header=False)

    # write summarized --> LINgroup output tsv
    if "LINgroup_report" in args.output_format:
        try:
            lingroups = tax_utils.read_lingroups(args.LINgroups)
        except ValueError as exc:
            error(f"ERROR: {str(exc)}")
            sys.exit(-1)

        lingroup_reportfile, limit_float = make_outfile(args.output_base, "lingroup_report", output_dir=args.output_dir)

        with FileOutputCSV(lingroup_reportfile) as out_fp:
            header, lgreport_results = single_query_results.make_lingroup_results(LINgroupsD = lingroups)
            tax_utils.write_output(header, lgreport_results, out_fp, sep="\t", write_header=True)


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
                       force=args.force, LIN_taxonomy=args.LIN_taxonomy)
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

    try:
         query_gather_results = tax_utils.check_and_load_gather_csvs(gather_csvs, tax_assign, force=args.force,
                                                                                       fail_on_missing_taxonomy=args.fail_on_missing_taxonomy,
                                                                                       keep_full_identifiers=args.keep_full_identifiers,
                                                                                       keep_identifier_versions = args.keep_identifier_versions,
                                                                                       )

    except ValueError as exc:
        error(f"ERROR: {str(exc)}")
        sys.exit(-1)

    if not query_gather_results:
        notify('No results for classification. Exiting.')
        sys.exit(-1)

    # for each queryResult, summarize at rank and classify according to thresholds, reporting any errors that occur.
    for queryResult in query_gather_results:
        try:
            queryResult.build_classification_result(rank=args.rank,
                                                    ani_threshold=args.ani_threshold,
                                                    containment_threshold=args.containment_threshold)

        except ValueError as exc:
            error(f"ERROR: {str(exc)}")
            sys.exit(-1)

    # write outputs
    if "csv_summary" in args.output_format:
        summary_outfile, limit_float = make_outfile(args.output_base, "classification", output_dir=args.output_dir)
        with FileOutputCSV(summary_outfile) as out_fp:
            tax_utils.write_summary(query_gather_results, out_fp, limit_float_decimals=limit_float, classification=True)

    # write summarized output in human-readable format
    if "human" in args.output_format:
        summary_outfile, limit_float = make_outfile(args.output_base, "human", output_dir=args.output_dir)

        with FileOutput(summary_outfile) as out_fp:
            tax_utils.write_human_summary(query_gather_results, out_fp, args.rank or "species", classification=True)

    # The following require a single rank:
    # note: interactive krona can handle mult ranks, do we want to enable?
    if "krona" in args.output_format:
        krona_results, header =  tax_utils.format_for_krona(query_gather_results=query_gather_results, rank=args.rank, classification=True)
        krona_outfile, limit_float = make_outfile(args.output_base, "krona", output_dir=args.output_dir)
        with FileOutputCSV(krona_outfile) as out_fp:
            tax_utils.write_krona(header, krona_results, out_fp)

    if "lineage_csv" in args.output_format:
        lineage_outfile, _ = make_outfile(args.output_base, "lineage_csv",
                                          output_dir=args.output_dir)
        lineage_results = []
        header = None
        for q_res in query_gather_results:
            if not header:
                ranks = list(q_res.ranks)
                if 'strain' in ranks: # maintains prior functionality.. but we could keep strain now, i think?
                    ranks.remove('strain')
                header = ["ident", *ranks]
            lineageD = q_res.classification_result.as_lineage_dict(q_res.query_info, ranks)
            lineage_results.append(lineageD)
        with FileOutputCSV(lineage_outfile) as out_fp:
            tax_utils.write_output(header, lineage_results, out_fp)
    



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
        query_gather_results = tax_utils.check_and_load_gather_csvs(gather_csvs, tax_assign, force=args.force,
                                                                                       fail_on_missing_taxonomy=args.fail_on_missing_taxonomy,
                                                                                       keep_full_identifiers=args.keep_full_identifiers,
                                                                                       keep_identifier_versions = args.keep_identifier_versions,
                                                                                       )

        if not query_gather_results:
            continue

        out_base = os.path.basename(g_csv.rsplit('.csv')[0])
        this_outfile, limit_float = make_outfile(out_base, "annotate", output_dir=args.output_dir)

        header = [field.name for field in fields(GatherRow)]
        with FileOutputCSV(this_outfile) as out_fp:
            header.append("lineage")
            w = csv.DictWriter(out_fp, header, delimiter=',')
            w.writeheader()

            for gather_res in query_gather_results:
                for taxres in gather_res.raw_taxresults:
                    gr = asdict(taxres.raw)
                    write_gr = {key: gr[key] for key in gr if key in header}
                    write_gr['lineage'] = taxres.lineageInfo.display_lineage(truncate_empty=True)
                    w.writerow(write_gr)



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
        for lp in lineage:
            if select_rank is None or lp.rank == select_rank:
                if pattern.search(lp.name):
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
        for vv in v:
            name = vv.name
            rank = vv.rank
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
