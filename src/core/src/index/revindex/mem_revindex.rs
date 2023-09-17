use std::sync::atomic::{AtomicUsize, Ordering};

use camino::Utf8Path as Path;
use camino::Utf8PathBuf as PathBuf;
use log::{debug, info};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::collection::Collection;
use crate::encodings::{Colors, Idx};
use crate::index::linear::LinearIndex;
use crate::index::revindex::HashToColor;
use crate::index::{GatherResult, Index, SigCounter};
use crate::prelude::*;
use crate::signature::{Signature, SigsTrait};
use crate::sketch::minhash::KmerMinHash;
use crate::sketch::Sketch;
use crate::Result;

pub struct RevIndex {
    linear: LinearIndex,
    hash_to_color: HashToColor,
    colors: Colors,
}

impl LinearIndex {
    fn index(
        self,
        threshold: usize,
        merged_query: Option<KmerMinHash>,
        queries: Option<&[KmerMinHash]>,
    ) -> RevIndex {
        let processed_sigs = AtomicUsize::new(0);

        #[cfg(feature = "parallel")]
        let sig_iter = self.collection().par_iter();

        #[cfg(not(feature = "parallel"))]
        let sig_iter = self.collection().iter();

        let filtered_sigs = sig_iter.enumerate().filter_map(|(dataset_id, _)| {
            let i = processed_sigs.fetch_add(1, Ordering::SeqCst);
            if i % 1000 == 0 {
                info!("Processed {} reference sigs", i);
            }

            let search_sig = self
                .collection()
                .sig_for_dataset(dataset_id as Idx)
                .expect("Error loading sig")
                .into();

            RevIndex::map_hashes_colors(
                dataset_id as Idx,
                &search_sig,
                queries,
                &merged_query,
                threshold,
                self.template(),
            )
        });

        #[cfg(feature = "parallel")]
        let (hash_to_color, colors) = filtered_sigs.reduce(
            || (HashToColor::new(), Colors::default()),
            HashToColor::reduce_hashes_colors,
        );

        #[cfg(not(feature = "parallel"))]
        let (hash_to_color, colors) = filtered_sigs.fold(
            (HashToColor::new(), Colors::default()),
            HashToColor::reduce_hashes_colors,
        );

        RevIndex {
            hash_to_color,
            colors,
            linear: self,
        }
    }
}

impl RevIndex {
    pub fn new(
        search_sigs: &[PathBuf],
        selection: &Selection,
        threshold: usize,
        queries: Option<&[KmerMinHash]>,
        _keep_sigs: bool,
    ) -> Result<RevIndex> {
        // If threshold is zero, let's merge all queries and save time later
        let merged_query = queries.and_then(|qs| Self::merge_queries(qs, threshold));

        let collection = Collection::from_paths(search_sigs)?.select(&selection)?;
        let linear = LinearIndex::from_collection(collection.try_into()?);

        Ok(linear.index(threshold, merged_query, queries))
    }

    pub fn from_zipfile<P: AsRef<Path>>(
        zipfile: P,
        selection: &Selection,
        threshold: usize,
        queries: Option<&[KmerMinHash]>,
        _keep_sigs: bool,
    ) -> Result<RevIndex> {
        // If threshold is zero, let's merge all queries and save time later
        let merged_query = queries.and_then(|qs| Self::merge_queries(qs, threshold));

        let collection = Collection::from_zipfile(zipfile)?.select(&selection)?;
        let linear = LinearIndex::from_collection(collection.try_into()?);

        Ok(linear.index(threshold, merged_query, queries))
    }

    fn merge_queries(qs: &[KmerMinHash], threshold: usize) -> Option<KmerMinHash> {
        if threshold == 0 {
            let mut merged = qs[0].clone();
            for query in &qs[1..] {
                merged.merge(query).unwrap();
            }
            Some(merged)
        } else {
            None
        }
    }

    pub fn new_with_sigs(
        search_sigs: Vec<Signature>,
        selection: &Selection,
        threshold: usize,
        queries: Option<&[KmerMinHash]>,
    ) -> Result<RevIndex> {
        // If threshold is zero, let's merge all queries and save time later
        let merged_query = queries.and_then(|qs| Self::merge_queries(qs, threshold));

        let collection = Collection::from_sigs(search_sigs)?.select(selection)?;
        let linear = LinearIndex::from_collection(collection.try_into()?);

        let idx = linear.index(threshold, merged_query, queries);

        Ok(idx)
    }

    fn map_hashes_colors(
        dataset_id: Idx,
        search_sig: &Signature,
        queries: Option<&[KmerMinHash]>,
        merged_query: &Option<KmerMinHash>,
        threshold: usize,
        template: &Sketch,
    ) -> Option<(HashToColor, Colors)> {
        let mut search_mh = None;
        if let Some(Sketch::MinHash(mh)) = search_sig.select_sketch(template) {
            search_mh = Some(mh);
        }

        let search_mh = search_mh.expect("Couldn't find a compatible MinHash");
        let mut hash_to_color = HashToColor::new();
        let mut colors = Colors::default();

        if let Some(qs) = queries {
            if let Some(ref merged) = merged_query {
                let (matched_hashes, intersection) = merged.intersection(search_mh).unwrap();
                if !matched_hashes.is_empty() || intersection > threshold as u64 {
                    hash_to_color.add_to(&mut colors, dataset_id, matched_hashes);
                }
            } else {
                for query in qs {
                    let (matched_hashes, intersection) = query.intersection(search_mh).unwrap();
                    if !matched_hashes.is_empty() || intersection > threshold as u64 {
                        hash_to_color.add_to(&mut colors, dataset_id, matched_hashes);
                    }
                }
            }
        } else {
            let matched = search_mh.mins();
            let size = matched.len() as u64;
            if !matched.is_empty() || size > threshold as u64 {
                hash_to_color.add_to(&mut colors, dataset_id, matched);
            }
        };

        if hash_to_color.is_empty() {
            None
        } else {
            Some((hash_to_color, colors))
        }
    }

    pub fn search(
        &self,
        counter: SigCounter,
        similarity: bool,
        threshold: usize,
    ) -> Result<Vec<String>> {
        self.linear.search(counter, similarity, threshold)
    }

    pub fn gather(
        &self,
        mut counter: SigCounter,
        threshold: usize,
        query: &KmerMinHash,
    ) -> Result<Vec<GatherResult>> {
        let mut match_size = usize::max_value();
        let mut matches = vec![];

        while match_size > threshold && !counter.is_empty() {
            let (dataset_id, size) = counter.most_common()[0];
            match_size = if size >= threshold { size } else { break };
            let result = self
                .linear
                .gather_round(dataset_id, match_size, query, matches.len())?;
            if let Some(Sketch::MinHash(match_mh)) =
                result.match_.select_sketch(self.linear.template())
            {
                // Prepare counter for finding the next match by decrementing
                // all hashes found in the current match in other datasets
                for hash in match_mh.iter_mins() {
                    if let Some(color) = self.hash_to_color.get(hash) {
                        counter.subtract(self.colors.indices(color).cloned());
                    }
                }
                counter.remove(&dataset_id);
                matches.push(result);
            } else {
                unimplemented!()
            }
        }
        Ok(matches)
    }

    pub fn template(&self) -> Sketch {
        self.linear.template().clone()
    }

    // TODO: mh should be a sketch, or even a sig...
    pub(crate) fn find_signatures(
        &self,
        mh: &KmerMinHash,
        threshold: f64,
        containment: bool,
        _ignore_scaled: bool,
    ) -> Result<Vec<(f64, Signature, String)>> {
        // TODO: proper threshold calculation
        let threshold: usize = (threshold * (mh.size() as f64)) as _;

        let counter = self.counter_for_query(mh);

        debug!(
            "number of matching signatures for hashes: {}",
            counter.len()
        );

        let mut results = vec![];
        for (dataset_id, size) in counter.most_common() {
            let match_size = if size >= threshold { size } else { break };

            let match_sig = self.linear.sig_for_dataset(dataset_id)?;
            let match_path = self
                .linear
                .collection()
                .record_for_dataset(dataset_id)?
                .internal_location();

            let mut match_mh = None;
            if let Some(Sketch::MinHash(mh)) = match_sig.select_sketch(self.linear.template()) {
                match_mh = Some(mh);
            }
            let match_mh = match_mh.unwrap();

            if size >= threshold {
                let score = if containment {
                    size as f64 / mh.size() as f64
                } else {
                    size as f64 / (mh.size() + match_size - size) as f64
                };
                let filename = match_path.to_string();
                let mut sig: Signature = match_sig.clone().into();
                sig.reset_sketches();
                sig.push(Sketch::MinHash(match_mh.clone()));
                results.push((score, sig, filename));
            } else {
                break;
            };
        }
        Ok(results)
    }

    pub fn counter_for_query(&self, query: &KmerMinHash) -> SigCounter {
        query
            .iter_mins()
            .filter_map(|hash| self.hash_to_color.get(hash))
            .flat_map(|color| self.colors.indices(color))
            .cloned()
            .collect()
    }
}

impl<'a> Index<'a> for RevIndex {
    type Item = Signature;

    fn insert(&mut self, _node: Self::Item) -> Result<()> {
        unimplemented!()
    }

    fn save<P: AsRef<std::path::Path>>(&self, _path: P) -> Result<()> {
        unimplemented!()
    }

    fn load<P: AsRef<std::path::Path>>(_path: P) -> Result<()> {
        unimplemented!()
    }

    fn len(&self) -> usize {
        self.linear.len()
    }

    fn signatures(&self) -> Vec<Self::Item> {
        self.linear
            .signatures()
            .into_iter()
            .map(|sig| sig.into())
            .collect()
    }

    fn signature_refs(&self) -> Vec<&Self::Item> {
        unimplemented!()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use crate::index::revindex::prepare_query;
    use crate::sketch::minhash::max_hash_for_scaled;
    use crate::Result;

    #[test]
    fn revindex_new() -> Result<()> {
        let selection = Selection::builder().ksize(31).scaled(10000).build();
        let search_sigs = [
            "../../tests/test-data/gather/GCF_000006945.2_ASM694v2_genomic.fna.gz.sig".into(),
            "../../tests/test-data/gather/GCF_000007545.1_ASM754v1_genomic.fna.gz.sig".into(),
        ];
        let index = RevIndex::new(&search_sigs, &selection, 0, None, false)?;
        assert_eq!(index.colors.len(), 3);

        Ok(())
    }

    #[test]
    fn revindex_many() -> Result<()> {
        let selection = Selection::builder().ksize(31).scaled(10000).build();
        let search_sigs = [
            "../../tests/test-data/gather/GCF_000006945.2_ASM694v2_genomic.fna.gz.sig".into(),
            "../../tests/test-data/gather/GCF_000007545.1_ASM754v1_genomic.fna.gz.sig".into(),
            "../../tests/test-data/gather/GCF_000008105.1_ASM810v1_genomic.fna.gz.sig".into(),
        ];

        let index = RevIndex::new(&search_sigs, &selection, 0, None, false)?;
        //dbg!(&index.linear.collection().manifest);
        /*
        dbg!(&index.colors.colors);
         0: 86
         1: 132
         2: 91
         (0, 1): 53
         (0, 2): 90
         (1, 2): 26
         (0, 1, 2): 261
         union: 739

        */
        //assert_eq!(index.colors.len(), 3);
        assert_eq!(index.colors.len(), 7);

        Ok(())
    }

    #[test]
    fn revindex_from_sigs() -> Result<()> {
        let selection = Selection::builder().ksize(31).scaled(10000).build();
        let search_sigs: Vec<Signature> = [
            "../../tests/test-data/gather/GCF_000006945.2_ASM694v2_genomic.fna.gz.sig",
            "../../tests/test-data/gather/GCF_000007545.1_ASM754v1_genomic.fna.gz.sig",
            "../../tests/test-data/gather/GCF_000008105.1_ASM810v1_genomic.fna.gz.sig",
        ]
        .into_iter()
        .map(|path| Signature::from_path(path).unwrap().swap_remove(0))
        .collect();

        let index = RevIndex::new_with_sigs(search_sigs, &selection, 0, None)?;
        /*
         dbg!(&index.colors.colors);
         0: 86
         1: 132
         2: 91
         (0, 1): 53
         (0, 2): 90
         (1, 2): 26
         (0, 1, 2): 261
         union: 739
        */
        //assert_eq!(index.colors.len(), 3);
        assert_eq!(index.colors.len(), 7);

        Ok(())
    }

    #[test]
    fn revindex_from_zipstorage() -> Result<()> {
        let selection = Selection::builder()
            .ksize(19)
            .scaled(100)
            .moltype(crate::encodings::HashFunctions::murmur64_protein)
            .build();
        let index = RevIndex::from_zipfile(
            "../../tests/test-data/prot/protein.zip",
            &selection,
            0,
            None,
            false,
        )
        .expect("error building from ziptorage");

        assert_eq!(index.colors.len(), 3);

        let query_sig = Signature::from_path(
            "../../tests/test-data/prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig",
        )
        .expect("Error processing query")
        .swap_remove(0)
        .select(&selection)?;

        let mut query_mh = None;
        if let Some(q) = prepare_query(query_sig, &selection) {
            query_mh = Some(q);
        }
        let query_mh = query_mh.expect("Couldn't find a compatible MinHash");

        let counter_rev = index.counter_for_query(&query_mh);
        let counter_lin = index.linear.counter_for_query(&query_mh);

        let results_rev = index.search(counter_rev, false, 0).unwrap();
        let results_linear = index.linear.search(counter_lin, false, 0).unwrap();
        assert_eq!(results_rev, results_linear);

        let counter_rev = index.counter_for_query(&query_mh);
        let counter_lin = index.linear.counter_for_query(&query_mh);

        let results_rev = index.gather(counter_rev, 0, &query_mh).unwrap();
        let results_linear = index.linear.gather(counter_lin, 0, &query_mh).unwrap();
        assert_eq!(results_rev.len(), 1);
        assert_eq!(results_rev, results_linear);

        Ok(())
    }
}
