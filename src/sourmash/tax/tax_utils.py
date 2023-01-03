"""
Utility functions for taxonomy analysis tools.
"""
import os
import csv
from collections import namedtuple, defaultdict
from collections import abc
import gzip

from sourmash import sqlite_utils, sourmash_args
from sourmash.exceptions import IndexNotSupported
from sourmash.distance_utils import containment_to_distance

import sqlite3


__all__ = ['get_ident', 'ascending_taxlist', 'collect_gather_csvs',
           'load_gather_results', 'check_and_load_gather_csvs',
           'find_match_lineage', 'summarize_gather_at',
           'find_missing_identities', 'make_krona_header',
           'aggregate_by_lineage_at_rank', 'format_for_krona',
           'write_krona', 'write_summary', 'write_classifications',
           'combine_sumgather_csvs_by_lineage', 'write_lineage_sample_frac',
           'MultiLineageDB']

from sourmash.logging import notify
from sourmash.sourmash_args import load_pathlist_from_file

# CTB: these could probably usefully be converted into dataclasses.
QueryInfo = namedtuple("QueryInfo", "query_md5, query_filename, query_bp, query_hashes, total_weighted_hashes")
SummarizedGatherResult = namedtuple("SummarizedGatherResult", "query_name, rank, fraction, lineage, query_md5, query_filename, f_weighted_at_rank, bp_match_at_rank, query_ani_at_rank, total_weighted_hashes")
ClassificationResult = namedtuple("ClassificationResult", "query_name, status, rank, fraction, lineage, query_md5, query_filename, f_weighted_at_rank, bp_match_at_rank, query_ani_at_rank")

# Essential Gather column names that must be in gather_csv to allow `tax` summarization
EssentialGatherColnames = ('query_name', 'name', 'f_unique_weighted', 'f_unique_to_query', 'unique_intersect_bp', 'remaining_bp', 'query_md5', 'query_filename')

# import lca utils as needed for now
from sourmash.lca import lca_utils
from sourmash.lca.lca_utils import (LineagePair, taxlist, display_lineage, pop_to_rank)


def get_ident(ident, *,
              keep_full_identifiers=False, keep_identifier_versions=False):
    # split identifiers = split on whitespace
    # keep identifiers = don't split .[12] from assembly accessions
    "Hack and slash identifiers."
    if not keep_full_identifiers:
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
    gather_csvs = []
    # ignore command line duplicates
    for gf in cmdline_gather_input:
        if gf not in gather_csvs:
            gather_csvs.append(gf)
        else:
            notify(f'ignoring duplicated reference to file: {gf}')
    # ignore pathlist duplicates
    if from_file:
        more_files = load_pathlist_from_file(from_file)
        for gf in more_files:
            if gf not in gather_csvs:
                gather_csvs.append(gf)
            else:
               notify(f'ignoring duplicated reference to file: {gf}')
    return gather_csvs


def load_gather_results(gather_csv, *, essential_colnames=EssentialGatherColnames,
                        seen_queries=None, force=False):
    "Load a single gather csv"
    if not seen_queries:
        seen_queries=set()
    header = []
    gather_results = []
    gather_queries = set()
    with sourmash_args.FileInputCSV(gather_csv) as r:
        header = r.fieldnames
        # check for empty file
        if not header:
            raise ValueError(f"Cannot read gather results from '{gather_csv}'. Is file empty?")

        # check for critical column names used by summarize_gather_at
        if not set(essential_colnames).issubset(header):
            raise ValueError(f"Not all required gather columns are present in '{gather_csv}'.")

        for n, row in enumerate(r):
            query_name = row['query_name']
            # check if we've seen this query already in a different gather CSV
            if query_name in seen_queries:
                # do not allow loading of same query from a second CSV.
                raise ValueError(f"Gather query {query_name} was found in more than one CSV. Cannot load from '{gather_csv}'.")
            else:
                gather_results.append(row)
            # add query name to the gather_queries from this CSV
            if query_name not in gather_queries:
                gather_queries.add(query_name)

    if not gather_results:
        raise ValueError(f'No gather results loaded from {gather_csv}.')
    else:
        notify(f"loaded {len(gather_results)} gather results from '{gather_csv}'.")
    return gather_results, header, gather_queries


def check_and_load_gather_csvs(gather_csvs, tax_assign, *, fail_on_missing_taxonomy=False, force=False):
    '''
    Load gather csvs, checking for empties and ids missing from taxonomic assignments.
    '''
    if not isinstance(gather_csvs, list):
        gather_csvs = [gather_csvs]
    gather_results = []
    total_missed = 0
    all_ident_missed = set()
    seen_queries = set()
    header = []
    n_ignored = 0
    for n, gather_csv in enumerate(gather_csvs):
        these_results = []
        try:
            these_results, header, seen_queries = load_gather_results(gather_csv, seen_queries=seen_queries, force=force)
        except ValueError as exc:
            if force:
                if "found in more than one CSV" in str(exc):
                    notify('Cannot force past duplicated gather query. Exiting.')
                    raise
                notify(str(exc))
                notify('--force is set. Attempting to continue to next set of gather results.')
                n_ignored+=1
                continue
            else:
                notify('Exiting.')
                raise

        # check for match identites in these gather_results not found in lineage spreadsheets
        ident_missed = find_missing_identities(these_results, tax_assign)
        if ident_missed:
            notify(f'The following are missing from the taxonomy information: {",".join(ident_missed)}')
            if fail_on_missing_taxonomy:
                raise ValueError('Failing on missing taxonomy, as requested via --fail-on-missing-taxonomy.')

            total_missed += len(ident_missed)
            all_ident_missed.update(ident_missed)
        # add these results to gather_results
        gather_results += these_results

    num_gather_csvs_loaded = n+1 - n_ignored
    notify(f'loaded {len(gather_results)} results total from {str(num_gather_csvs_loaded)} gather CSVs')

    return gather_results, all_ident_missed, total_missed, header


def find_match_lineage(match_ident, tax_assign, *, skip_idents = [],
                       keep_full_identifiers=False,
                       keep_identifier_versions=False):
    lineage=""
    match_ident = get_ident(match_ident, keep_full_identifiers=keep_full_identifiers, keep_identifier_versions=keep_identifier_versions)
    # if identity not in lineage database, and not --fail-on-missing-taxonomy, skip summarizing this match
    if match_ident in skip_idents:
        return lineage
    try:
        lineage = tax_assign[match_ident]
    except KeyError:
        raise ValueError(f"ident {match_ident} is not in the taxonomy database.")
    return lineage


def summarize_gather_at(rank, tax_assign, gather_results, *, skip_idents = [],
                        keep_full_identifiers=False,
                        keep_identifier_versions=False, best_only=False,
                        seen_perfect=set(),
                        estimate_query_ani=False):
    """
    Summarize gather results at specified taxonomic rank
    """
    # init dictionaries
    sum_uniq_weighted = defaultdict(lambda: defaultdict(float))
    # store together w/ ^ instead?
    sum_uniq_to_query = defaultdict(lambda: defaultdict(float))
    sum_uniq_bp = defaultdict(lambda: defaultdict(float))
    query_info = {}

    set_ksize = False
    ksize, scaled, query_nhashes = None, 0, None

    for row in gather_results:
        # get essential gather info
        if not set_ksize and "ksize" in row.keys():
            set_ksize = True
            ksize = int(row['ksize'])
            scaled = int(row['scaled'])
        
        query_name = row['query_name']
        f_unique_to_query = float(row['f_unique_to_query'])
        f_uniq_weighted = float(row['f_unique_weighted'])
        unique_intersect_bp = int(row['unique_intersect_bp'])
        total_weighted_hashes = int(row.get('total_weighted_hashes', 0))
        query_md5 = row['query_md5']
        query_filename = row['query_filename']
        # get query_bp
        if query_name not in query_info.keys(): #REMOVING THIS AFFECTS GATHER RESULTS!!! BUT query bp should always be same for same query? bug?
            if "query_nhashes" in row.keys():
                query_nhashes = int(row["query_nhashes"])
            if "query_bp" in row.keys():
                query_bp = int(row["query_bp"])
            else:
                query_bp = unique_intersect_bp + int(row['remaining_bp'])
        
        # store query info
        query_info[query_name] = QueryInfo(query_md5=query_md5, query_filename=query_filename, query_bp=query_bp, query_hashes=query_nhashes, total_weighted_hashes=total_weighted_hashes)
        
        if estimate_query_ani and (not ksize or not scaled):
            if not set_ksize:
                estimate_query_ani=False
                notify("WARNING: Please run gather with sourmash >= 4.4 to estimate query ANI at rank. Continuing without ANI...")
        
        match_ident = row['name']

        # 100% match? are we looking at something in the database?
        if f_unique_to_query >= 1.0 and query_name not in seen_perfect: # only want to notify once, not for each rank
            ident = get_ident(match_ident,
                              keep_full_identifiers=keep_full_identifiers,
                              keep_identifier_versions=keep_identifier_versions)
            seen_perfect.add(query_name)
            notify(f'WARNING: 100% match! Is query "{query_name}" identical to its database match, {ident}?')

        # get lineage for match
        lineage = find_match_lineage(match_ident, tax_assign,
                                    skip_idents=skip_idents,
                                    keep_full_identifiers=keep_full_identifiers,
                                    keep_identifier_versions=keep_identifier_versions)
        # ident was in skip_idents
        if not lineage:
            continue

        # summarize at rank!
        lineage = pop_to_rank(lineage, rank)
        assert lineage[-1].rank == rank, lineage[-1]
        # record info
        sum_uniq_to_query[query_name][lineage] += f_unique_to_query
        sum_uniq_weighted[query_name][lineage] += f_uniq_weighted
        sum_uniq_bp[query_name][lineage] += unique_intersect_bp

    # sort and store each as SummarizedGatherResult
    sum_uniq_to_query_sorted = []
    for query_name, lineage_weights in sum_uniq_to_query.items():
        qInfo = query_info[query_name]
        sumgather_items = list(lineage_weights.items())
        sumgather_items.sort(key = lambda x: -x[1])
        query_ani = None
        if best_only:
            lineage, fraction = sumgather_items[0]
            if fraction > 1:
                raise ValueError(f"The tax summary of query '{query_name}' is {fraction}, which is > 100% of the query!! This should not be possible. Please check that your input files come directly from a single gather run per query.")
            elif fraction == 0:
                continue
            f_weighted_at_rank = sum_uniq_weighted[query_name][lineage]
            bp_intersect_at_rank = sum_uniq_bp[query_name][lineage]
            if estimate_query_ani:
                query_ani = containment_to_distance(fraction, ksize, scaled,
                                                    n_unique_kmers= qInfo.query_hashes, sequence_len_bp= qInfo.query_bp).ani
            sres = SummarizedGatherResult(query_name, rank, fraction, lineage, qInfo.query_md5,
                                          qInfo.query_filename, f_weighted_at_rank, bp_intersect_at_rank, query_ani, qInfo.total_weighted_hashes * scaled)
            sum_uniq_to_query_sorted.append(sres)
        else:
            total_f_weighted= 0.0
            total_f_classified = 0.0
            total_bp_classified = 0
            for lineage, fraction in sumgather_items:
                query_ani = None
                if fraction > 1:
                    raise ValueError(f"The tax summary of query '{query_name}' is {fraction}, which is > 100% of the query!! This should not be possible. Please check that your input files come directly from a single gather run per query.")
                elif fraction == 0:
                    continue
                total_f_classified += fraction
                f_weighted_at_rank = sum_uniq_weighted[query_name][lineage]
                total_f_weighted += f_weighted_at_rank
                bp_intersect_at_rank = int(sum_uniq_bp[query_name][lineage])
                total_bp_classified += bp_intersect_at_rank
                if estimate_query_ani:
                    query_ani = containment_to_distance(fraction, ksize, scaled,
                                                        n_unique_kmers=qInfo.query_hashes, sequence_len_bp=qInfo.query_bp).ani
                sres = SummarizedGatherResult(query_name, rank, fraction, lineage, query_md5,
                                              query_filename, f_weighted_at_rank, bp_intersect_at_rank, query_ani, qInfo.total_weighted_hashes * scaled)
                sum_uniq_to_query_sorted.append(sres)

            # record unclassified
            lineage = ()
            query_ani = None
            fraction = 1.0 - total_f_classified
            if fraction > 0:
                f_weighted_at_rank = 1.0 - total_f_weighted
                bp_intersect_at_rank = qInfo.query_bp - total_bp_classified
                sres = SummarizedGatherResult(query_name, rank, fraction, lineage, query_md5,
                                              query_filename, f_weighted_at_rank, bp_intersect_at_rank, query_ani, qInfo.total_weighted_hashes*scaled)
                sum_uniq_to_query_sorted.append(sres)

    return sum_uniq_to_query_sorted, seen_perfect, estimate_query_ani


def find_missing_identities(gather_results, tax_assign):
    """
    Identify match ids/accessions from gather results
    that are not present in taxonomic assignments.
    """
    ident_missed= set()
    for row in gather_results:
        match_ident = row['name']
        match_ident = get_ident(match_ident)
        if match_ident not in tax_assign:
            ident_missed.add(match_ident)

    if ident_missed:
        notify(f'of {len(gather_results)} gather results, missed {len(ident_missed)} lineage assignments.')
    return ident_missed


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
    for res in rank_results:
        if res.query_name not in all_queries:
            all_queries.append(res.query_name)
        if by_query:
            lineage_summary[res.lineage][res.query_name] = res.fraction
        else:
            lineage_summary[res.lineage] += res.fraction
    return lineage_summary, all_queries, len(all_queries)


def format_for_krona(rank, summarized_gather):
    '''
    Aggregate list of SummarizedGatherResults and format for krona output
    '''
    num_queries=0
    for res_rank, rank_results in summarized_gather.items():
        if res_rank == rank:
            lineage_summary, all_queries, num_queries = aggregate_by_lineage_at_rank(rank_results, by_query=False)
    # if aggregating across queries divide fraction by the total number of queries
    for lin, fraction in lineage_summary.items():
        # divide total fraction by total number of queries
        lineage_summary[lin] = fraction/num_queries

    # sort by fraction
    lin_items = list(lineage_summary.items())
    lin_items.sort(key = lambda x: -x[1])

    # reformat lineage for krona_results printing
    krona_results = []
    unclassified_fraction = 0
    for lin, fraction in lin_items:
        # save unclassified fraction for the end
        if lin == ():
            unclassified_fraction = fraction
            continue
        lin_list = display_lineage(lin).split(';')
        krona_results.append((fraction, *lin_list))

    # handle unclassified
    if unclassified_fraction:
        len_unclassified_lin = len(krona_results[-1]) -1
        unclassifed_lin = ["unclassified"]*len_unclassified_lin
        krona_results.append((unclassified_fraction, *unclassifed_lin))

    return krona_results


def write_krona(rank, krona_results, out_fp, *, sep='\t'):
    'write krona output'
    # CTB: do we want to optionally allow restriction to a specific rank
    # & above?
    header = make_krona_header(rank)
    tsv_output = csv.writer(out_fp, delimiter='\t')
    tsv_output.writerow(header)
    for res in krona_results:
        tsv_output.writerow(res)


def write_summary(summarized_gather, csv_fp, *, sep=',', limit_float_decimals=False):
    '''
    Write taxonomy-summarized gather results for each rank.
    '''
    header = SummarizedGatherResult._fields
    w = csv.DictWriter(csv_fp, header, delimiter=sep)
    w.writeheader()
    for rank, rank_results in summarized_gather.items():
        for res in rank_results:
            rD = res._asdict()
            if limit_float_decimals:
                rD['fraction'] = f'{res.fraction:.3f}'
                rD['f_weighted_at_rank'] = f'{res.f_weighted_at_rank:.3f}'
            rD['lineage'] = display_lineage(res.lineage)
            if rD['lineage'] == "":
                rD['lineage'] = "unclassified"
            w.writerow(rD)


def write_kreport(summarized_gather, csv_fp, *, sep='\t'):
    '''
    Write taxonomy-summarized gather results as kraken-style kreport.

    While this format typically records the percent of number of reads assigned to taxa,
    we can create comparable output by reporting the percent of k-mers (percent containment)
    and the total number of k-mers matched.

    standard reads-based `kreport` columns:
    - `Percent Reads Contained in Taxon`: The cumulative percentage of reads for this taxon and all descendants.
    - `Number of Reads Contained in Taxon`: The cumulative number of reads for this taxon and all descendants.
    - `Number of Reads Assigned to Taxon`: The number of reads assigned directly to this taxon (not a cumulative count of all descendants).
    - `Rank Code`: (U)nclassified, (R)oot, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. 
    - `NCBI Taxon ID`: Numerical ID from the NCBI taxonomy database.
    - `Scientific Name`: The scientific name of the taxon.

    Example reads-based `kreport` with all columns:
    ```
    88.41	2138742	193618	K	2	Bacteria
    0.16	3852	818	P	201174	  Actinobacteria
    0.13	3034	0	C	1760	    Actinomycetia
    0.13	3034	45	O	85009	      Propionibacteriales
    0.12	2989	1847	F	31957	        Propionibacteriaceae
    0.05	1142	352	G	1912216	          Cutibacterium
    0.03	790	790	S	1747	            Cutibacterium acnes
    ```

    sourmash `kreport` caveats:
    - `Percent k-mers Contained in Taxon`: weighted by k-mer abundance
    - `Estimated bp Contained in Taxon`: NOT WEIGHTED BY ABUNDANCE
    - `Number of Reads Assigned to Taxon` and `NCBI Taxon ID` will not be reported (blank entries).

    In the future, we may wish to report the NCBI taxid when we can (NCBI taxonomy only).
    '''
    columns = ["percent_containment", "num_bp_contained", "num_bp_assigned", "rank_code", "ncbi_taxid", "sci_name"]
    w = csv.DictWriter(csv_fp, columns, delimiter=sep)

    rankCode = { "superkingdom": "D", "kingdom": "K", "phylum": "P", "class": "C",
                 "order": "O", "family":"F", "genus": "G", "species": "S"} # , "": "U"

    # check - are we using v4.5.0 or later gather CSVs?
    for rank, rank_results in summarized_gather.items():
        for res in rank_results:
            if res.total_weighted_hashes == 0:
                raise ValueError("ERROR: cannot produce 'kreport' format from gather results before sourmash v4.5.0")

    unclassified_written=False
    for rank, rank_results in summarized_gather.items():
        rcode = rankCode[rank]
        for res in rank_results:
            # SummarizedGatherResults have an unclassified lineage at every rank, to facilitate reporting at a specific rank.
            # Here, we only need to report it once, since it will be the same fraction for all ranks
            if not res.lineage:
                rank_sciname = "unclassified"
                rcode = "U"
                # if we've already written the unclassified portion, skip and continue to next loop iteration
                if unclassified_written:
                    continue
                else:
                    unclassified_written=True
            else:
                rank_sciname = res.lineage[-1].name
            kresD = {"rank_code": rcode, "ncbi_taxid": "", "sci_name": rank_sciname,  "num_bp_assigned": 0}
            # total percent containment, weighted to include abundance info
            proportion = res.f_weighted_at_rank * 100
            kresD['percent_containment'] = f'{proportion:.2f}'
            # weighted bp
            kresD["num_bp_contained"] = int(res.f_weighted_at_rank * res.total_weighted_hashes)
            if rank == 'species' or rank_sciname == "unclassified":
                kresD["num_bp_assigned"] = kresD["num_bp_contained"]
            w.writerow(kresD)


def write_human_summary(summarized_gather, out_fp, display_rank):
    '''
    Write human-readable taxonomy-summarized gather results for a specific rank.
    '''
    header = SummarizedGatherResult._fields

    found_ANI = False
    results = [] 
    for rank, rank_results in summarized_gather.items():
        # only show results for a specified rank.
        if rank == display_rank:
            rank_results = list(rank_results)
            rank_results.sort(key=lambda res: -res.f_weighted_at_rank)

            for res in rank_results:
                rD = res._asdict()
                rD['fraction'] = f'{res.fraction:.3f}'
                rD['f_weighted_at_rank'] = f"{res.f_weighted_at_rank*100:>4.1f}%"
                if rD['query_ani_at_rank'] is not None:
                    found_ANI = True
                    rD['query_ani_at_rank'] = f"{res.query_ani_at_rank*100:>3.1f}%"
                else:
                    rD['query_ani_at_rank'] = '-    '
                rD['lineage'] = display_lineage(res.lineage)
                if rD['lineage'] == "":
                    rD['lineage'] = "unclassified"

                results.append(rD)


    if found_ANI:
        out_fp.write("sample name    proportion   cANI   lineage\n")
        out_fp.write("-----------    ----------   ----   -------\n")

        for rD in results:
            out_fp.write("{query_name:<15s}   {f_weighted_at_rank}     {query_ani_at_rank}  {lineage}\n".format(**rD))
    else:
        out_fp.write("sample name    proportion   lineage\n")
        out_fp.write("-----------    ----------   -------\n")

        for rD in results:
            out_fp.write("{query_name:<15s}   {f_weighted_at_rank}     {lineage}\n".format(**rD))


def write_lineage_csv(summarized_gather, csv_fp):
    '''
    Write a lineage-CSV format file suitable for use with sourmash tax ... -t.
    '''
    ranks = lca_utils.taxlist(include_strain=False)
    header = ['ident', *ranks]
    w = csv.DictWriter(csv_fp, header)
    w.writeheader()
    for rank, rank_results in summarized_gather.items():
        for res in rank_results:
            d = {}
            d[rank] = ""
            for rank, name in res.lineage:
                d[rank] = name

            d['ident'] = res.query_name
            w.writerow(d)


def write_classifications(classifications, csv_fp, *, sep=',', limit_float_decimals=False):
    '''
    Write taxonomy-classifed gather results.
    '''
    header = ClassificationResult._fields
    w = csv.DictWriter(csv_fp, header, delimiter=sep)
    w.writeheader()
    for rank, rank_results in classifications.items():
        for res in rank_results:
            rD = res._asdict()
            if limit_float_decimals:
                rD['fraction'] = f'{res.fraction:.3f}'
                rD['f_weighted_at_rank'] = f'{res.f_weighted_at_rank:.3f}'
            rD['lineage'] = display_lineage(res.lineage)
            # needed?
            if rD['lineage'] == "":
                rD['lineage'] = "unclassified"
            w.writerow(rD)


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
    unclassified_row = None
    for lin, sampleinfo in sorted(lineage_dict.items()):
        if format_lineage:
            lin = display_lineage(lin)

        #add lineage and 0 placeholders
        row = {'lineage': lin}
        row.update(blank_row)
        # add info for query_names that exist for this lineage
        row.update(sampleinfo)
        # if unclassified, save this row for the end
        if not lin:
            row.update({'lineage': 'unclassified'})
            unclassified_row = row
            continue
        # write row
        w.writerow(row)
    if unclassified_row:
        w.writerow(unclassified_row)


class LineageDB(abc.Mapping):
    "Base LineageDB class built around an assignments dictionary."
    def __init__(self, assign_d, avail_ranks):
        self.assignments = assign_d
        self.available_ranks = set(avail_ranks)

    def __getitem__(self, ident):
        "Retrieve the lineage tuple for identifer (or raise KeyError)"
        return self.assignments[ident]

    def __iter__(self):
        "Return all identifiers for this db."
        return iter(self.assignments)

    def __len__(self):
        "Return number of lineages"
        return len(self.assignments)

    def __bool__(self):
        "Are there any lineages at all in this database?"
        return bool(self.assignments)

    @classmethod
    def load(cls, filename, *, delimiter=',', force=False,
             keep_full_identifiers=False, keep_identifier_versions=True):
        """
        Load a taxonomy assignment CSV file into a LineageDB.

        'keep_full_identifiers=False' will split identifiers from strings
        using whitespace, e.g. 'IDENT other name stuff' => 'IDENT'

        'keep_identifier_versions=False' will remove trailing versions,
        e.g. 'IDENT.1' => 'IDENT'.
        """
        include_strain=False
        if not keep_identifier_versions and keep_full_identifiers:
            raise ValueError("keep_identifer_versions=False doesn't make sense with keep_full_identifiers=True")

        if not os.path.exists(filename):
            raise ValueError(f"'{filename}' does not exist")

        if os.path.isdir(filename):
            raise ValueError(f"'{filename}' is a directory")

        with sourmash_args.FileInputCSV(filename) as r:
            header = r.fieldnames
            if not header:
                raise ValueError(f'cannot read taxonomy assignments from {filename}')

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
                elif 'name' in header and 'lineage' in header:
                    return cls.load_from_gather_with_lineages(filename,
                                                              force=force)
                else:
                    header_str = ",".join([repr(x) for x in header])
                    raise ValueError(f'No taxonomic identifiers found; headers are {header_str}')
            # is "strain" an available rank?
            if "strain" in header:
                include_strain=True

            # check that all ranks are in header
            ranks = list(lca_utils.taxlist(include_strain=include_strain))
            if not set(ranks).issubset(header):
                # for now, just raise err if not all ranks are present.
                # in future, we can define `ranks` differently if desired
                # return them from this function so we can check the `available` ranks
                raise ValueError('Not all taxonomy ranks present')

            assignments = {}
            num_rows = 0
            n_species = 0
            n_strains = 0

            # now parse and load lineages
            for n, row in enumerate(r):
                num_rows += 1
                lineage = []
                # read row into a lineage pair
                for rank in lca_utils.taxlist(include_strain=include_strain):
                    lin = row[rank]
                    lineage.append(LineagePair(rank, lin))
                ident = row[identifier]

                # fold, spindle, and mutilate ident?
                ident = get_ident(ident,
                                  keep_full_identifiers=keep_full_identifiers,
                                  keep_identifier_versions=keep_identifier_versions)

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

        return LineageDB(assignments, ranks)


    @classmethod
    def load_from_gather_with_lineages(cls, filename, *, force=False):
        """
        Load an annotated gather-with-lineages CSV file produced by
        'tax annotate' into a LineageDB.
        """
        include_strain = False

        if not os.path.exists(filename):
            raise ValueError(f"'{filename}' does not exist")

        if os.path.isdir(filename):
            raise ValueError(f"'{filename}' is a directory")

        with sourmash_args.FileInputCSV(filename) as r:
            header = r.fieldnames
            if not header:
                raise ValueError(f'cannot read taxonomy assignments from {filename}')

            if "name" not in header or "lineage" not in header:
                raise ValueError(f"Expected headers 'name' and 'lineage' not found. Is this a with-lineages file?")

            ranks = list(lca_utils.taxlist(include_strain=include_strain))
            assignments = {}
            num_rows = 0
            n_species = 0
            n_strains = 0

            # now parse and load lineages
            for n, row in enumerate(r):
                num_rows += 1

                name = row['name']
                ident = get_ident(name)
                lineage = row['lineage']
                lineage = lca_utils.make_lineage(lineage)

                # check duplicates
                if ident in assignments:
                    if assignments[ident] != tuple(lineage):
                        # this should not happen with valid
                        # sourmash tax annotate output, but check anyway.
                        if not force:
                            raise ValueError(f"multiple lineages for identifier {ident}")
                else:
                    assignments[ident] = tuple(lineage)

                    if lineage[-1].rank == 'species':
                        n_species += 1
                    elif lineage[-1].rank == 'strain':
                        n_species += 1
                        n_strains += 1

        return LineageDB(assignments, ranks)


class LineageDB_Sqlite(abc.Mapping):
    """
    A LineageDB based on a sqlite3 database with a 'sourmash_taxonomy' table.
    """
    # NOTE: 'order' is a reserved name in sql, so we have to use 'order_'.
    columns = ('superkingdom', 'phylum', 'order_', 'class', 'family',
               'genus', 'species', 'strain')
    table_name = 'sourmash_taxonomy'

    def __init__(self, conn, *, table_name=None):
        self.conn = conn

        # provide for legacy support for pre-sourmash_internal days...
        if table_name is not None:
            self.table_name = table_name

        # check that the right table is there.
        c = conn.cursor()
        try:
            c.execute(f'SELECT * FROM {self.table_name} LIMIT 1')
        except (sqlite3.DatabaseError, sqlite3.OperationalError):
            raise ValueError("not a taxonomy database")
            
        # check: can we do a 'select' on the right table?
        self.__len__()
        c = conn.cursor()

        # get available ranks...
        ranks = set()
        for column, rank in zip(self.columns, taxlist(include_strain=True)):
            query = f'SELECT COUNT({column}) FROM {self.table_name} WHERE {column} IS NOT NULL AND {column} != ""'
            c.execute(query)
            cnt, = c.fetchone()
            if cnt:
                ranks.add(rank)

        self.available_ranks = ranks
        self.cursor = c

    @classmethod
    def load(cls, location):
        "load taxonomy information from an existing sqlite3 database"
        conn = sqlite_utils.open_sqlite_db(location)
        if not conn:
            raise ValueError("not a sqlite taxonomy database")

        table_name = None
        c = conn.cursor()
        try:
            info = sqlite_utils.get_sourmash_internal(c)
        except sqlite3.OperationalError:
            info = {}

        if 'SqliteLineage' in info:
            if info['SqliteLineage'] != '1.0':
                raise IndexNotSupported

            table_name = 'sourmash_taxonomy'
        else:
            # legacy support for old taxonomy DB, pre sourmash_internal.
            try:
                c.execute('SELECT * FROM taxonomy LIMIT 1')
                table_name = 'taxonomy'
            except sqlite3.OperationalError:
                pass

        if table_name is None:
            raise ValueError("not a sqlite taxonomy database")

        return cls(conn, table_name=table_name)

    def _make_tup(self, row):
        "build a tuple of LineagePairs for this sqlite row"
        tup = [ LineagePair(n, r) for (n, r) in zip(taxlist(True), row) ]
        return tuple(tup)

    def __getitem__(self, ident):
        "Retrieve lineage for identifer"
        c = self.cursor
        c.execute(f'SELECT superkingdom, phylum, class, order_, family, genus, species, strain FROM {self.table_name} WHERE ident=?', (ident,))

        # retrieve names list...
        names = c.fetchone()
        if names:
            # ...and construct lineage tuple
            tup = self._make_tup(names)
            while tup and not tup[-1].name:
                tup = tup[:-1]

            return tup

        raise KeyError(ident)

    def __bool__(self):
        "Do we have any info?"
        return bool(len(self))

    def __len__(self):
        "Return number of rows"
        c = self.conn.cursor()
        c.execute(f'SELECT COUNT(DISTINCT ident) FROM {self.table_name}')
        nrows, = c.fetchone()
        return nrows

    def __iter__(self):
        "Return all identifiers"
        # create new cursor so as to allow other operations
        c = self.conn.cursor()
        c.execute(f'SELECT DISTINCT ident FROM {self.table_name}')

        for ident, in c:
            yield ident

    def items(self):
        "return all items in the sqlite database"
        c = self.conn.cursor()

        c.execute(f'SELECT DISTINCT ident, superkingdom, phylum, class, order_, family, genus, species, strain FROM {self.table_name}')

        for ident, *names in c:
            yield ident, self._make_tup(names)


class MultiLineageDB(abc.Mapping):
    "A wrapper for (dynamically) combining multiple lineage databases."

    # NTP: currently, later lineage databases will override earlier ones.
    # Do we want to report/summarize shadowed identifiers?

    def __init__(self):
        self.lineage_dbs = []

    @property
    def available_ranks(self):
        "build the union of available ranks across all databases"
        # CTB: do we need to worry about lineages of shadowed identifiers?
        x = set()
        for db in self.lineage_dbs:
            x.update(db.available_ranks)
        return x

    def add(self, db):
        "Add a new lineage database"
        self.lineage_dbs.insert(0, db)

    def __iter__(self):
        "Return all identifiers (once)"
        seen = set()
        for db in self.lineage_dbs:
            for k in db:
                if k not in seen:
                    seen.add(k)
                    yield k

    def items(self):
        "Return all (identifiers, lineage_tup), masking duplicate idents"
        seen = set()
        for db in self.lineage_dbs:
            for k, v in db.items():
                if k not in seen:
                    seen.add(k)
                    yield k, v

    def shadowed_identifiers(self):
        seen = set()
        dups = set()
        for db in self.lineage_dbs:
            for k, v in db.items():
                if k in seen:
                    dups.add(k)
                else:
                    seen.add(k)
        return seen

    def __getitem__(self, ident):
        "Return lineage tuple for first match to identifier."
        for db in self.lineage_dbs:
            if ident in db:
                return db[ident]

        # not found? KeyError!
        raise KeyError(ident)

    def __len__(self):
        "Return number of distinct identifiers. Currently iterates over all."
        # CTB: maybe we can make this unnecessary?
        x = set(self)
        return len(x)

    def __bool__(self):
        "True if any contained database has content."
        return any( bool(db) for db in self.lineage_dbs )

    def save(self, filename_or_fp, file_format):
        assert file_format in ('sql', 'csv')

        is_filename = False
        try:
            filename_or_fp.write
        except AttributeError:
            is_filename = True

        if file_format == 'sql':
            if not is_filename:
                raise ValueError("file format '{file_format}' requires a filename, not a file handle")
            self._save_sqlite(filename_or_fp)
        elif file_format == 'csv':
            # we need a file handle; open file.
            fp = filename_or_fp
            if is_filename:
                if filename_or_fp.endswith('.gz'):
                    fp = gzip.open(filename_or_fp, 'wt', newline="")
                else:
                    fp = open(filename_or_fp, 'w', newline="")

            try:
                self._save_csv(fp)
            finally:
                # close the file we opened!
                if is_filename:
                    fp.close()

    def _save_sqlite(self, filename, *, conn=None):
        from sourmash import sqlite_utils

        if conn is None:
            db = sqlite3.connect(filename)
        else:
            assert not filename
            db = conn

        cursor = db.cursor()
        try:
            sqlite_utils.add_sourmash_internal(cursor, 'SqliteLineage', '1.0')
        except sqlite3.OperationalError:
            raise ValueError("attempt to write a readonly database")

        try:
            # CTB: could add 'IF NOT EXIST' here; would need tests, too.
            cursor.execute("""

        CREATE TABLE sourmash_taxonomy (
            ident TEXT NOT NULL,
            superkingdom TEXT,
            phylum TEXT,
            class TEXT,
            order_ TEXT,
            family TEXT,
            genus TEXT,
            species TEXT,
            strain TEXT
        )
        """)
            did_create = True
        except sqlite3.OperationalError:
            # already exists?
            raise ValueError(f"taxonomy table already exists in '{filename}'")

        # follow up and create index
        cursor.execute("CREATE UNIQUE INDEX sourmash_taxonomy_ident ON sourmash_taxonomy(ident);")
        for ident, tax in self.items():
            x = [ident, *[ t.name for t in tax ]]

            # fill the taxonomy tuple with empty values until it's the
            # right length for the SQL statement -
            while len(x) < 9:
                x.append('')

            cursor.execute('INSERT INTO sourmash_taxonomy (ident, superkingdom, phylum, class, order_, family, genus, species, strain) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)', x)

        db.commit()

    def _save_csv(self, fp):
        headers = ['identifiers'] + list(taxlist(include_strain=True))
        w = csv.DictWriter(fp, fieldnames=headers)
        w.writeheader()

        for n, (ident, tax) in enumerate(self.items()):
            row = {}
            row['identifiers'] = ident

            # convert tax LineagePairs into dictionary
            for t in tax:
                row[t.rank] = t.name

            # add strain if needed
            if 'strain' not in row:
                row['strain'] = ''

            w.writerow(row)

    @classmethod
    def load(cls, locations, **kwargs):
        "Load one or more taxonomies from the given location(s)"
        force = kwargs.get('force', False)

        if isinstance(locations, str):
            raise TypeError("'locations' should be a list, not a string")

        tax_assign = cls()
        for location in locations:
            # try faster formats first
            loaded = False

            # sqlite db?
            try:
                this_tax_assign = LineageDB_Sqlite.load(location)
                loaded = True
            except ValueError:
                pass

            # CSV file?
            if not loaded:
                try:
                    this_tax_assign = LineageDB.load(location, **kwargs)
                    loaded = True
                except (ValueError, csv.Error) as exc:
                    # for the last loader, just pass along ValueError...
                    if not force:
                        raise ValueError(f"cannot read taxonomy assignments from '{location}': {str(exc)}")

            # nothing loaded, goodbye!
            if not loaded and not force:
                raise ValueError(f"cannot read taxonomy assignments from '{location}'")

            if loaded:
                tax_assign.add(this_tax_assign)

        return tax_assign
