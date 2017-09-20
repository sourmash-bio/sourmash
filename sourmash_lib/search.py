from collections import namedtuple

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
