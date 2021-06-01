"""
Utility functions for taxonomy analysis tools.
"""
import csv
from os.path import exists
from collections import namedtuple, defaultdict, Counter
import itertools

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
    sum_uniq_weighted = defaultdict(float)
    for row in gather_results:
        # move these checks to loading function!
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
        sum_uniq_weighted[lineage] += f_uniq_weighted

    items = list(sum_uniq_weighted.items())
    items.sort(key = lambda x: -x[1])
    if best_only:
        return [items[0]]# return list to keep formatting the same as non best-only
    return items

def find_missing_identities(gather_results, tax_assign):
    n_missed = 0
    ident_missed= []
    for row in gather_results:
        match_ident = row['name']
        match_ident = get_ident(match_ident)
        if match_ident not in tax_assign:
            n_missed += 1
            ident_missed.append(match_ident)

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

def format_for_krona(rank, summarized_gather):
    krona_results = []
    for gather_rank, rank_results in summarized_gather.items():
        if gather_rank == rank:
            for sorted_result in rank_results:
                lin,fraction = sorted_result
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

def format_tax_to_frac(taxonomy_csvs, rank, output_csv):
    '''
    takes the output for sourmash taxonomy summarize and produces a 
    tab-separated file with fractions for each sample. Sample names
    are based on csv file names, with ".csv" removed
    lineage	sample1	sample2	sample3	
    lin_a	.4	.17	.6
    lin_b	0	0	.1
    lin_c	0	.3	0
    lin_d	.2	.1	0
    lin_e	0	0	.01
    lin_f	0	.07	0
    lin_g	0	0	0
    lin_h	.3	.4	.2
    '''
    samples = [csv.split(".")[0] for csv in csvs]

    possible_ranks = ['superkingdom', "phylum", "class", "order", "family", "genus", "species"]
    if rank not in possible_ranks:
        raise ValueError(f"Rank {rank} not available")

    
    lineage_dict = {}
    sample_name_dict = {}
    seen_lineages = set()

    # create dictionary that holds all of the sample names
    for file in csvs:
        sample_name = file.split('.')[0]
        sample_name_dict[sample_name] = 0 

    for file in csvs:
        with open(file, 'r') as fp:
            r = csv.DictReader(fp)
            for n, row in enumerate(r):
                if row["rank"] == rank:
                    seen_lineages.add(row["lineage"])
            fp.close()

    for lineage in seen_lineages:
        lineage_dict[lineage] = sample_name_dict.copy()

    for sample in sample_name_dict:
        with open(sample + ".csv", "r") as fp:
            r = csv.DictReader(fp)
            for n, row in enumerate(r):
                if row["rank"] == rank:
                    lineage = (row["lineage"])
                    fraction = (row["fraction"])
                    lineage_dict[lineage][sample] = fraction
            fp.close()


    samples.insert(0, "lineage")
    with open(output_csv, 'w') as f_output:
        w = csv.DictWriter(f_output, samples)
        w.writeheader()
        for key,val in sorted(lineage_dict.items()):
            row = {'lineage': key}
            row.update(val)
            w.writerow(row)
