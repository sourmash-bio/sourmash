from collections import namedtuple

import sourmash_lib
from .signature import SourmashSignature
from .sbtmh import (search_minhashes,
                                search_minhashes_containment)
from .sbtmh import SearchMinHashesFindBest

# generic SearchResult across individual signatures + SBTs.
SearchResult = namedtuple('SearchResult',
                          'similarity, match_sig, md5, filename, name')

def search_databases(query, databases, threshold, do_containment, best_only):
    # set up the search & score function(s) - similarity vs containment
    search_fn = search_minhashes
    query_match = lambda x: query.similarity(x, downsample=True)
    if do_containment:
        search_fn = search_minhashes_containment
        query_match = lambda x: query.contained_by(x, downsample=True)

    results = []
    found_md5 = set()
    for (sbt_or_siglist, filename, is_sbt) in databases:
        if is_sbt:
            if best_only:            # this needs to be reset for each SBT
                search_fn = SearchMinHashesFindBest().search

            tree = sbt_or_siglist
            for leaf in tree.find(search_fn, query, threshold):
                similarity = query_match(leaf.data)
                if similarity >= threshold and \
                       leaf.data.md5sum() not in found_md5:
                    sr = SearchResult(similarity=similarity,
                                      match_sig=leaf.data,
                                      md5=leaf.data.md5sum(),
                                      filename=filename,
                                      name=leaf.data.name())
                    found_md5.add(sr.md5)
                    results.append(sr)

        else: # list of signatures
            for ss in sbt_or_siglist:
                similarity = query_match(ss)
                if similarity >= threshold and \
                       ss.md5sum() not in found_md5:
                    sr = SearchResult(similarity=similarity,
                                      match_sig=ss,
                                      md5=ss.md5sum(),
                                      filename=filename,
                                      name=ss.name())
                    found_md5.add(sr.md5)
                    results.append(sr)


    # sort results on similarity (reverse)
    results.sort(key=lambda x: -x.similarity)

    return results


def gather_databases(query, databases, threshold_bp):
    from sourmash_lib.sbtmh import SearchMinHashesFindBestIgnoreMaxHash

    orig_query = query
    orig_mins = orig_query.minhash.get_hashes()

    # calculate the band size/resolution R for the genome
    R_metagenome = orig_query.minhash.scaled

    # define a function to do a 'best' search and get only top match.
    def find_best(dblist, query):
        results = []
        for (sbt_or_siglist, filename, is_sbt) in dblist:
            search_fn = SearchMinHashesFindBestIgnoreMaxHash().search

            if is_sbt:
                tree = sbt_or_siglist

                for leaf in tree.find(search_fn, query, 0.0):
                    leaf_e = leaf.data.minhash
                    similarity = query.minhash.similarity_ignore_maxhash(leaf_e)
                    if similarity > 0.0:
                        results.append((similarity, leaf.data))
            else:
                for ss in sbt_or_siglist:
                    similarity = query.minhash.similarity_ignore_maxhash(ss.minhash)
                    if similarity > 0.0:
                        results.append((similarity, ss))

        if not results:
            return None, None, None

        # take the best result
        results.sort(key=lambda x: -x[0])   # reverse sort on similarity
        best_similarity, best_leaf = results[0]
        return best_similarity, best_leaf, filename


    # define a function to build new signature object from set of mins
    def build_new_signature(mins, template_sig):
        e = template_sig.minhash.copy_and_clear()
        e.add_many(mins)
        return SourmashSignature(e)

    # construct a new query that doesn't have the max_hash attribute set.
    new_mins = query.minhash.get_hashes()
    query = build_new_signature(new_mins, orig_query)

    sum_found = 0.
    GatherResult = namedtuple('GatherResult',
                               'intersect_bp, f_orig_query, f_match, f_unique_to_query, filename, name, md5, leaf')
    while 1:
        best_similarity, best_leaf, filename = find_best(databases, query)
        if not best_leaf:          # no matches at all!
            break

        # subtract found hashes from search hashes, construct new search
        query_mins = set(query.minhash.get_hashes())
        found_mins = best_leaf.minhash.get_hashes()

        # figure out what the resolution of the banding on the genome is,
        # based either on an explicit --scaled parameter, or on genome
        # cardinality (deprecated)
        if not best_leaf.minhash.max_hash:
            error('Best hash match in sbt_gather has no max_hash')
            error('Please prepare database of sequences with --scaled')
            sys.exit(-1)

        R_genome = best_leaf.minhash.scaled

        # pick the highest R / lowest resolution
        R_comparison = max(R_metagenome, R_genome)

        # CTB: these could probably be replaced by minhash.downsample_scaled.
        new_max_hash = sourmash_lib.MAX_HASH / float(R_comparison)
        query_mins = set([ i for i in query_mins if i < new_max_hash ])
        found_mins = set([ i for i in found_mins if i < new_max_hash ])
        orig_mins = set([ i for i in orig_mins if i < new_max_hash ])

        # calculate intersection:
        intersect_mins = query_mins.intersection(found_mins)
        intersect_orig_mins = orig_mins.intersection(found_mins)
        intersect_bp = R_comparison * len(intersect_orig_mins)

        if intersect_bp < threshold_bp:   # hard cutoff for now
            notify('found less than {} in common. => exiting',
                   format_bp(intersect_bp))
            break

        # calculate fractions wrt first denominator - genome size
        genome_n_mins = len(found_mins)
        f_match = len(intersect_mins) / float(genome_n_mins)
        f_orig_query = len(intersect_orig_mins) / float(len(orig_mins))

        # calculate fractions wrt second denominator - metagenome size
        query_n_mins = len(orig_query.minhash.get_hashes())
        f_unique_to_query = len(intersect_mins) / float(query_n_mins)

        result = GatherResult(intersect_bp=intersect_bp,
                              f_orig_query=f_orig_query,
                              f_match=f_match,
                              f_unique_to_query=f_unique_to_query,
                              filename=filename,
                              md5=best_leaf.md5sum(),
                              name=best_leaf.name(),
                              leaf=best_leaf)

        # construct a new query, minus the previous one.
        query_mins -= set(found_mins)
        query = build_new_signature(query_mins, orig_query)

        yield result, len(intersect_mins), new_max_hash, query
