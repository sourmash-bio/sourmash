"""
Utility functions for taxonomy analysis tools.
"""
import sys
import csv
from os.path import exists, basename, dirname, abspath
from collections import namedtuple, defaultdict, Counter

__all__ = ['get_ident', 'load_gather_results',
           'summarize_gather_at', 'find_missing_identities']

from sourmash.logging import notify, error, debug
from sourmash.sourmash_args import load_pathlist_from_file

SummarizedGatherResult = namedtuple("SummarizedGatherResult", "query_name, rank, fraction, lineage")

# import lca utils as needed for now
from sourmash.lca import lca_utils
from sourmash.lca.lca_utils import (LineagePair, build_tree, find_lca,
                                    taxlist, count_lca_for_assignments,
                                    zip_lineage, display_lineage,
                                    make_lineage, is_lineage_match,
                                    pop_to_rank)


def get_ident(ident, *, split_identifiers=True, keep_identifier_versions=False):
    # split identifiers = split on whitespace
    # keep identifiers = don't split .[12] from assembly accessions
    "Hack and slash identifiers."
    if split_identifiers:
        ident = ident.split(' ')[0]
        if not keep_identifier_versions:
            ident = ident.split('.')[0]
    return ident


def ascending_taxlist(include_strain=True):
    """
    Provide an ordered list of taxonomic ranks: strain --> superkingdom
    """
    ascending_taxlist = ['species', 'genus', 'family', 'order',
                         'class', 'phylum', 'superkingdom']
    if include_strain:
        ascending_taxlist = ['strain'] + ascending_taxlist
    for k in ascending_taxlist:
        yield k


def collect_gather_csvs(cmdline_gather_input, *, from_file=None):
    """
    collect gather files from cmdline; --from-file input
    """
    gather_csvs = cmdline_gather_input
    if from_file:
        more_files = load_pathlist_from_file(from_file)
        for gf in more_files:
            if gf not in gather_csvs:
                gather_csvs.append(gf)
    return gather_csvs


def load_gather_results(gather_csv):
    "Load a single gather csv"
    header = []
    gather_results = []
    with open(gather_csv, 'rt') as fp:
        r = csv.DictReader(fp)
        #do we want to check for critical column names?
        for n, row in enumerate(r):
            if not header:
                header= list(row.keys())
            gather_results.append(row)
    if not gather_results:
        raise ValueError(f'No gather results loaded from {gather_csv}.')
    else:
        notify(f'loaded {len(gather_results)} gather results.')
    return gather_results, header


def check_and_load_gather_csvs(gather_csvs, tax_assign, *, fail_on_missing_taxonomy=False, force=False):
    '''
    Load gather csvs, checking for empties and ids missing from taxonomic assignments.
    '''
    if not isinstance(gather_csvs, list):
        gather_csvs = [gather_csvs]
    gather_results = []
    total_missed = 0
    all_ident_missed = set()
    header = []
    for gather_csv in gather_csvs:
        these_results = []
        try:
            these_results, header = load_gather_results(gather_csv)
        except ValueError:
            if force:
                notify(f'--force is set. Attempting to continue to next set of gather results.')
                continue
            else:
                notify(f'Exiting.')
                raise

        # check for match identites in these gather_results not found in lineage spreadsheets
        n_missed, ident_missed = find_missing_identities(these_results, tax_assign)
        if n_missed:
            notify(f'The following are missing from the taxonomy information: {",".join(ident_missed)}')
            if fail_on_missing_taxonomy:
                raise ValueError(f'Failing on missing taxonomy, as requested via --fail-on-missing-taxonomy.')

            total_missed += n_missed
            all_ident_missed.update(ident_missed)
        # add these results to gather_results
        gather_results += these_results

    return gather_results, all_ident_missed, total_missed, header


def find_match_lineage(match_ident, tax_assign, *, skip_idents = [], split_identifiers=True, keep_identifier_versions=False):
    lineage=""
    match_ident = get_ident(match_ident, split_identifiers=split_identifiers, keep_identifier_versions=keep_identifier_versions)
    # if identity not in lineage database, and not --fail-on-missing-taxonomy, skip summarizing this match
    if match_ident in skip_idents:
        return lineage
    try:
        lineage = tax_assign[match_ident]
    except KeyError:
        raise ValueError(f"ident {match_ident} is not in the taxonomy database.")
    return lineage


def summarize_gather_at(rank, tax_assign, gather_results, *, skip_idents = [], split_identifiers=True, keep_identifier_versions=False, best_only=False):
    """
    Summarize gather results at specified taxonomic rank
    """
    sum_uniq_weighted = defaultdict(lambda: defaultdict(float))
    for row in gather_results:
        query_name = row['query_name']
        match_ident = row['name']
        # get lineage for match
        lineage = find_match_lineage(match_ident, tax_assign, skip_idents = skip_idents, split_identifiers=split_identifiers, keep_identifier_versions=keep_identifier_versions)
        # ident was in skip_idents
        if not lineage:
            continue

        # summarize at rank!
        lineage = pop_to_rank(lineage, rank)
        assert lineage[-1].rank == rank, lineage[-1]

        f_uniq_weighted = row['f_unique_weighted']
        f_uniq_weighted = float(f_uniq_weighted)
        sum_uniq_weighted[query_name][lineage] += f_uniq_weighted

    # sort and store each as SummarizedGatherResult
    sum_uniq_weighted_sorted = []
    for query_name, lineage_weights in sum_uniq_weighted.items():
        query_results = []
        sumgather_items = list(lineage_weights.items())
        sumgather_items.sort(key = lambda x: -x[1])
        if best_only:
            lineage, fraction = sumgather_items[0]
            sum_uniq_weighted_sorted.append(SummarizedGatherResult(query_name, rank, fraction, lineage))
        else:
            for lineage, fraction in sumgather_items:
                sum_uniq_weighted_sorted.append(SummarizedGatherResult(query_name, rank, fraction, lineage))

    return sum_uniq_weighted_sorted


def find_missing_identities(gather_results, tax_assign):
    """
    Identify match ids/accessions from gather results
    that are not present in taxonomic assignments.
    """
    n_missed = 0
    ident_missed= set()
    for row in gather_results:
        match_ident = row['name']
        match_ident = get_ident(match_ident)
        if match_ident not in tax_assign:
            n_missed += 1
            ident_missed.add(match_ident)

    notify(f'of {len(gather_results)}, missed {n_missed} lineage assignments.')
    return n_missed, ident_missed


# pass ranks; have ranks=[default_ranks]
def make_krona_header(min_rank, *, include_strain=False):
    "make header for krona output"
    header = ["fraction"]
    tl = list(taxlist(include_strain=include_strain))
    try:
        rank_index = tl.index(min_rank)
    except ValueError:
        raise ValueError(f"Rank {min_rank} not present in available ranks!")
    return tuple(header + tl[:rank_index+1])


def aggregate_by_lineage_at_rank(rank_results, *, by_query=False):
    '''
    Aggregate list of rank SummarizedGatherResults,
    keeping query info or aggregating across queries.
    '''
    lineage_summary = defaultdict(float)
    if by_query:
        lineage_summary = defaultdict(dict)
    all_queries = []
    for (query_name, rank, fraction, lineage) in rank_results:
        if query_name not in all_queries:
            all_queries.append(query_name)
        if by_query:
            lineage_summary[lineage][query_name] = fraction
        else:
            lineage_summary[lineage] += fraction
    return lineage_summary, all_queries, len(all_queries)


def format_for_krona(rank, summarized_gather):
    '''
    Aggregate list of SummarizedGatherResults and format for krona output
    '''
    num_queries=0
    for res_rank, rank_results in summarized_gather.items():
        if res_rank == rank:
            lineage_summary, all_queries, num_queries = aggregate_by_lineage_at_rank(rank_results, by_query=False)
    # if multiple_samples, divide fraction by the total number of query files
    for lin, fraction in lineage_summary.items():
        # divide total fraction by total number of queries
        lineage_summary[lin] = fraction/num_queries

    # sort by fraction
    lin_items = list(lineage_summary.items())
    lin_items.sort(key = lambda x: -x[1])

    # reformat lineage for krona_results printing
    krona_results = []
    for lin, fraction in lin_items:
        lin_list = display_lineage(lin).split(';')
        krona_results.append((fraction, *lin_list))

    return krona_results


def write_krona(rank, krona_results, out_fp, *, sep='\t'):
    'write krona output'
    header = make_krona_header(rank)
    tsv_output = csv.writer(out_fp, delimiter='\t')
    tsv_output.writerow(header)
    for res in krona_results:
        tsv_output.writerow(res)


def write_summary(summarized_gather, csv_fp, *, sep='\t'):
    '''
    Write taxonomy-summarized gather results for each rank.
    '''
    header= ["query_name", "rank", "fraction", "lineage"]
    w = csv.writer(csv_fp)
    w.writerow(header)
    for rank, rank_results in summarized_gather.items():
        for (query_name, rank, fraction, lineage) in rank_results:
            w.writerow([query_name, rank, f'{fraction:.3f}', display_lineage(lineage)])


def write_classifications(classifications, csv_fp, *, sep='\t'):
    '''
    Write taxonomy-classifed gather results.
    '''
    header= ["query_name", "status", "rank", "fraction", "lineage"]
    w = csv.writer(csv_fp)
    w.writerow(header)
    for rank, rank_results in classifications.items():
        for (query_name, status, rank, fraction, lineage) in rank_results:
            w.writerow([query_name, status, rank, f'{fraction:.3f}', display_lineage(lineage)])


def combine_sumgather_csvs_by_lineage(gather_csvs, *, rank="species", accept_ranks = list(lca_utils.taxlist(include_strain=False)), force=False):
    '''
    Takes in one or more output csvs from `sourmash taxonomy summarize`
    and combines the results into a nested dictionary with lineages
    as the keys {lineage: {sample1: frac1, sample2: frac2}}.
    Uses the file basename (minus .csv extension) as sample identifier.

    usage:

        linD, all_samples = combine_sumgather_by_lineage(["sample1.csv", "sample2.csv"], rank="genus")

    output:

        linD = {lin_a: {'sample1': 0.4, 'sample2': 0.17, 'sample3': 0.6}
                lin_b: {'sample1': 0.0, 'sample2': 0.0,  'sample3': 0.1}
                lin_c: {'sample1': 0.3, 'sample2': 0.4,  'sample3': 0.2} }

        all_samples = ['sample1','sample2','sample3']

    '''
    if rank not in accept_ranks:
        raise ValueError(f"Rank {rank} not available.")

    sgD = defaultdict(dict)
    all_samples = []
    for g_csv in gather_csvs:
        # collect lineage info for this sample
        lineageD = defaultdict(list)
        with open(g_csv, 'r') as fp:
            r = csv.DictReader(fp)
            for row in r:
                if row["rank"] == rank:
                    query_name = row["query_name"]
                    lin = row["lineage"]
                    frac = row["fraction"]
                    if query_name not in all_samples:
                        all_samples.append(query_name)
                    sgD[lin][query_name] = frac
            fp.close()
    return sgD, all_samples


def write_lineage_sample_frac(sample_names, lineage_dict, out_fp, *, format_lineage=False, sep='\t'):
    '''
    takes in a lineage dictionary with sample counts (output of combine_sumgather_by_lineage)
    and produces a tab-separated file with fractions for each sample.

    input: {lin_a: {sample1: 0.4, sample2: 0.17, sample3: 0.6}
            lin_b: {sample1: 0.0, sample2: 0.0, sample3: 0.1}
            lin_c: {sample1: 0.3, sample2: 0.4, sample3: 0.2}}

    output:

    lineage    sample1	sample2	sample3
    lin_a	  0.4    0.17     0.6
    lin_b	  0.0    0.0      0.1
    lin_c	  0.3    0.4      0.2
    '''

    header = ["lineage"] + sample_names
    w = csv.DictWriter(out_fp, header, delimiter=sep)
    w.writeheader()
    blank_row = {query_name: 0 for query_name in sample_names}
    for lin, sampleinfo in sorted(lineage_dict.items()):
        if format_lineage:
            lin = display_lineage(lin)
        #add lineage and 0 placeholders
        row = {'lineage': lin}
        row.update(blank_row)
        # add info for query_names that exist for this lineage
        row.update(sampleinfo)
        # write row
        w.writerow(row)


def load_taxonomy_csv(filename, *, delimiter=',', force=False,
                              split_identifiers=False,
                              keep_identifier_versions=False):
    """
    Load a taxonomy assignment spreadsheet into a dictionary.

    The 'assignments' dictionary that's returned maps identifiers to
    lineage tuples.
    """
    include_strain=False
    fp = open(filename, newline='')
    r = csv.DictReader(fp, delimiter=delimiter)
    header = r.fieldnames
    if not header:
        raise ValueError(f'Cannot read taxonomy assignments from {filename}. Is file empty?')

    identifier = "ident"
    # check for ident/identifier, handle some common alternatives
    if "ident" not in header:
        # check for ident/identifier, handle some common alternatives
        if 'identifiers' in header:
            identifier = 'identifiers'
            header = ["ident" if "identifiers" == x else x for x in header]
        elif 'accession' in header:
            identifier = 'accession'
            header = ["ident" if "accession" == x else x for x in header]
        else:
            raise ValueError(f'No taxonomic identifiers found.')
    if "strain" in header:
        include_strain=True

   # check that all ranks are in header
    ranks = list(lca_utils.taxlist(include_strain=include_strain))
    if not set(ranks).issubset(header):
        # for now, just raise err if not all ranks are present.
        # in future, we can define `ranks` differently if desired
        # return them from this function so we can check the `available` ranks
        raise ValueError(f'Not all taxonomy ranks present')

    assignments = {}
    num_rows = 0
    n_species = 0
    n_strains = 0

    # now parse and load lineages
    for n, row in enumerate(r):
        if row: #and row[0].strip():        # want non-empty row
            num_rows += 1
            lineage = []
            # read row into a lineage pair
            for rank in lca_utils.taxlist(include_strain=include_strain):
                lineage.append(LineagePair(rank, row[rank]))
            ident = row[identifier]

            # fold, spindle, and mutilate ident?
            if split_identifiers:
                ident = ident.split(' ')[0]

                if not keep_identifier_versions:
                    ident = ident.split('.')[0]

            # clean lineage of null names, replace with 'unassigned'
            lineage = [ (a, lca_utils.filter_null(b)) for (a,b) in lineage ]
            lineage = [ LineagePair(a, b) for (a, b) in lineage ]

            # remove end nulls
            while lineage and lineage[-1].name == 'unassigned':
                lineage = lineage[:-1]

            # store lineage tuple
            if lineage:
                # check duplicates
                if ident in assignments:
                    if assignments[ident] != tuple(lineage):
                        if not force:
                            raise ValueError(f"multiple lineages for identifier {ident}")
                else:
                    assignments[ident] = tuple(lineage)

                    if lineage[-1].rank == 'species':
                        n_species += 1
                    elif lineage[-1].rank == 'strain':
                        n_species += 1
                        n_strains += 1

    fp.close()

    return assignments, num_rows, ranks
