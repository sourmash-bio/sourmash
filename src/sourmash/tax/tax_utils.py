"""
Utility functions for taxonomy analysis tools.
"""
import csv
from os.path import exists, basename
from collections import namedtuple, defaultdict, Counter

__all__ = ['get_ident', 'load_gather_results',
           'summarize_gather_at', 'find_missing_identities']

from sourmash.logging import notify, error, debug

# import lca utils as needed for now
from sourmash.lca import lca_utils
from sourmash.lca.lca_utils import (LineagePair, build_tree, find_lca,
                                    taxlist, count_lca_for_assignments,
                                    zip_lineage, display_lineage,
                                    make_lineage, is_lineage_match,
                                    pop_to_rank)


def get_ident(ident, split_identifiers=True, keep_identifier_versions=False):
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

def load_gather_files_from_csv(from_csv):
    gather_files = []
    seen = set()
    with open(from_csv, 'rt') as fp:
        r = csv.DictReader(fp, fieldnames=['name', 'filepath'])
        for n, row in enumerate(r):
            name = row["name"]
            if name in seen:
                notify(f"found duplicate name: {name}. Ignoring...")
            else:
                seen.add(name)
                gather_files.append((name, row["filepath"]))
    notify(f'loaded {len(gather_files)} gather files from csv input.')
    return gather_files, seen

# load and aggregate all gather results
def load_gather_results(gather_csv):
    gather_results = []
    with open(gather_csv, 'rt') as fp:
        r = csv.DictReader(fp)
        #todo: add a check for all gather column names
        for n, row in enumerate(r):
            gather_results.append(row)
    notify(f'loaded {len(gather_results)} gather results.')
    return gather_results


# this summarizes at a specific rank.
def summarize_gather_at(rank, tax_assign, gather_results, skip_idents = [], split_identifiers=True, keep_identifier_versions=False, best_only=False):
    # collect!
    #sum_uniq_weighted = defaultdict(float)
    # how to modify this to enable multiple gather csvs, AND multigather output? need to use query_name, summarize by query_name
    sum_uniq_weighted = defaultdict(lambda: defaultdict(float))
    for row in gather_results:
        query_name = row['query_name']
        # move these checks to loading function!?
        match_ident = row['name']
        match_ident = get_ident(match_ident, split_identifiers, keep_identifier_versions)
        # if identity not in lineage database, and not --fail-on-missing-taxonomy, skip summarizing this match
        if match_ident in skip_idents:
            continue
        try:
            lineage = tax_assign[match_ident]
        except KeyError:
            raise KeyError(f"ident {match_ident} is not in the taxonomy database.")
        # actual summarization code
        lineage = pop_to_rank(lineage, rank)
        assert lineage[-1].rank == rank, lineage[-1]

        f_uniq_weighted = row['f_unique_weighted']
        f_uniq_weighted = float(f_uniq_weighted)
        sum_uniq_weighted[query_name][lineage] += f_uniq_weighted

    sum_uniq_weighted_sorted = defaultdict(list)
    for query_name, lineage_weights in sum_uniq_weighted.items():
        items = list(lineage_weights.items())
        items.sort(key = lambda x: -x[1])
        if best_only:
            sum_uniq_weighted_sorted[query_name] = [items[0]] # list to keep formatting the same as non best-only
        else:
            sum_uniq_weighted_sorted[query_name] = items

    return sum_uniq_weighted_sorted


def find_missing_identities(gather_results, tax_assign):
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
def make_krona_header(min_rank, include_strain=False):
    header = ["fraction"]
    tl = list(taxlist(include_strain=include_strain))
    try:
        rank_index = tl.index(min_rank)
    except ValueError:
        raise ValueError(f"Rank {min_rank} not present in available ranks!")
    return tuple(header + tl[:rank_index+1])


# this is for summarized results of a single query
def aggregate_by_lineage_at_rank(rank_results):
    #query_lineage_summary = defaultdict(float)
    query_lineage_summary = Counter()
    for lin, fraction in rank_results:
        query_lineage_summary[lin] += fraction
    return query_lineage_summary


def format_and_summarize_for_krona(rank, summarized_gather):
    num_queries=0
    #krona_summary = defaultdict(float)
    krona_summary = Counter()
    for res_rank, query_summarized_gather in summarized_gather.items():
        if res_rank == rank:
            for query, sumgather in query_summarized_gather.items():
                num_queries += 1
                query_lineage_summary = aggregate_by_lineage_at_rank(sumgather)
            #add results from each query
            krona_summary.update(query_lineage_summary)

    # if multiple_samples, divide fraction by the total number of query files
    for lin, fraction in krona_summary.items():
        # add query-specific fraction (fraction/total num queries)
        krona_summary[lin] = fraction/num_queries

    # sort by fraction
    krona_items = list(krona_summary.items())
    krona_items.sort(key = lambda x: -x[1])

    # reformat lineage for krona_results printing
    krona_results = []
    for lin, fraction in krona_items:
        lin_list = display_lineage(lin).split(';')
        krona_results.append((fraction, *lin_list))

    return krona_results


def write_krona(rank, krona_results, out_fp, sep='\t'):
    header = make_krona_header(rank)
    tsv_output = csv.writer(out_fp, delimiter='\t')
    tsv_output.writerow(header)
    for res in krona_results:
        tsv_output.writerow(res)


def write_summary(summarized_gather, csv_fp, sep='\t'):
    header= ["rank", "fraction", "lineage"]
    w = csv.writer(csv_fp)
    w.writerow(header)
    for rank, rank_results in summarized_gather.items():
        for sorted_result in rank_results:
            lin,val = sorted_result
            w.writerow([rank, f'{val:.3f}', display_lineage(lin)])


def write_classifications(classifications, csv_fp, sep='\t'):
    header= ["query_name", "classification_rank", "fraction_matched_at_rank", "lineage"]
    w = csv.writer(csv_fp)
    w.writerow(header)
    for rank, rank_results in classifications.items():
        # do we want to sort the results somehow?
        #items = list(sum_uniq_weighted.items())
        #items.sort(key = lambda x: -x[1])
        for result in rank_results:
            name, (lin,val) = result
            w.writerow([rank, name, f'{val:.3f}', display_lineage(lin)])


def combine_sumgather_csvs_by_lineage(gather_csvs, rank="species", accept_ranks = list(lca_utils.taxlist(include_strain=False)), force=False):
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

    all_samples = [basename(g_csv).rsplit(".csv", 1)[0].rsplit('.summarized')[0] for g_csv in gather_csvs]

    # default dict to store lineage: {sample_id: fraction} info. better way to do this?
    sgD = defaultdict(lambda: {sample_id : 0.0 for sample_id in all_samples})
    for g_csv in gather_csvs:
        sample_id = basename(g_csv).rsplit(".csv", 1)[0].rsplit('.summarized')[0]

        # collect lineage info for this sample
        with open(g_csv, 'r') as fp:
            r = csv.DictReader(fp)
            for n, row in enumerate(r):
                if row["rank"] == rank:
                    lin = row["lineage"]
                    frac = row["fraction"]
                    sgD[lin][sample_id] = frac
            fp.close()
    return sgD, all_samples


def write_lineage_sample_frac(sample_names, lineage_dict, out_fp, sep='\t'):
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
    for lin, sampleinfo in sorted(lineage_dict.items()):
        row = {'lineage': lin}
        row.update(sampleinfo)
        w.writerow(row)


# see https://github.com/luizirber/2020-cami/blob/master/scripts/gather_to_opal.py
def write_cami_profiling_bioboxes_format(sample_id, ranks, taxons, out_fp, *, taxonomy_id=None, program=None, format_version="0.9.1", sep="\t"):
    # init version, not working yet
    header_title = "# Taxonomic Profiling Output"
    sample_info = f"@SampleID:{sample_id}"
    version_info = f" @Version:{format_version}"
    rank_info = f"@Ranks:{ranks}"
    output_lines = [header_title, sample_info, version_info, rank_info]
    if taxonomy_id is not None:
        output_lines.append(f"@TaxonomyID:{taxonomy_id}")
#    if program is not None:
#        output_lines.append(f"@__program__: {program}")
    output_lines.append(f"@@TAXID\tRANK\tTAXPATH\tPERCENTAGE") # actual tsv header

    for tax in taxons.itertuples(index=False, name=None):
        tax_line = "\t".join(str(t) for t in tax)
        output_lines.append(tax_line)

    #write instead of return!
    #return "\n".join(output_lines)
