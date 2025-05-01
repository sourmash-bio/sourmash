use std::sync::atomic::{AtomicUsize, Ordering};

use camino::Utf8Path as Path;
use camino::Utf8PathBuf as PathBuf;
use log::{debug, info};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::collection::Collection;
use crate::collection::CollectionSet;
use crate::encodings::{Colors, Idx};
use crate::index::linear::LinearIndex;
use crate::index::revindex::{
    self as module, CounterGather, DatasetPicklist, Datasets, DbStats, HashToColor, QueryColors,
    RevIndexOps,
};
use crate::index::{GatherResult, Index, SigCounter};
use crate::prelude::*;
use crate::signature::{Signature, SigsTrait};
use crate::sketch::minhash::{KmerMinHash, KmerMinHashBTree};
use crate::sketch::Sketch;
use crate::Result;
use crate::ScaledType;

pub struct MemRevIndex {
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
    ) -> MemRevIndex {
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

            MemRevIndex::map_hashes_colors(
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

        MemRevIndex {
            hash_to_color,
            colors,
            linear: self,
        }
    }
}

impl MemRevIndex {
    #[allow(clippy::new_ret_no_self)]
    pub fn new(
        search_sigs: &[PathBuf],
        selection: &Selection,
        threshold: usize,
        queries: Option<&[KmerMinHash]>,
    ) -> Result<module::RevIndex> {
        // If threshold is zero, let's merge all queries and save time later
        let merged_query = queries.and_then(|qs| Self::merge_queries(qs, threshold));

        let collection = Collection::from_paths(search_sigs)?.select(selection)?;
        let linear = LinearIndex::from_collection(collection.try_into()?);

        let idx = linear.index(threshold, merged_query, queries);

        Ok(module::RevIndex::Mem(idx))
    }

    pub fn from_zipfile<P: AsRef<Path>>(
        zipfile: P,
        selection: &Selection,
        threshold: usize,
        queries: Option<&[KmerMinHash]>,
    ) -> Result<module::RevIndex> {
        // If threshold is zero, let's merge all queries and save time later
        let merged_query = queries.and_then(|qs| Self::merge_queries(qs, threshold));

        let collection = Collection::from_zipfile(zipfile)?.select(selection)?;
        let linear = LinearIndex::from_collection(collection.try_into()?);

        let idx = linear.index(threshold, merged_query, queries);
        Ok(module::RevIndex::Mem(idx))
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
    ) -> Result<module::RevIndex> {
        // If threshold is zero, let's merge all queries and save time later
        let merged_query = queries.and_then(|qs| Self::merge_queries(qs, threshold));

        let collection = Collection::from_sigs(search_sigs)?.select(selection)?;
        let linear = LinearIndex::from_collection(collection.try_into()?);

        let idx = linear.index(threshold, merged_query, queries);

        Ok(module::RevIndex::Mem(idx))
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

    pub fn template(&self) -> Sketch {
        self.linear.template().clone()
    }

    pub fn scaled(&self) -> ScaledType {
        if let Sketch::MinHash(mh) = self.linear.template() {
            mh.scaled()
        } else {
            unimplemented!()
        }
    }
}

impl RevIndexOps for MemRevIndex {
    fn location(&self) -> &str {
        ""
    }

    fn counter_for_query(
        &self,
        query: &KmerMinHash,
        picklist: Option<DatasetPicklist>,
    ) -> SigCounter {
        query
            .iter_mins()
            .filter_map(|hash| self.hash_to_color.get(hash))
            .flat_map(|color| self.colors.indices(color))
            .filter_map(|idx| {
                if let Some(pl) = &picklist {
                    if pl.dataset_ids.contains(idx) {
                        Some(idx)
                    } else {
                        None
                    }
                } else {
                    Some(idx)
                }
            })
            .cloned()
            .collect()
    }

    /// build a CounterGather struct for a particular query
    fn prepare_gather_counters(
        &self,
        query: &KmerMinHash,
        picklist: Option<DatasetPicklist>,
    ) -> CounterGather {
        let counter = self.counter_for_query(query, picklist);
        let hash_to_color = self.hash_to_color.clone();

        // restrict hash_to_color to hashes contained in query
        let hash_to_color: HashToColor = query
            .iter_mins()
            .filter_map(|&hash| {
                let color = hash_to_color.get(&hash);
                color.map(|c| (hash, *c))
            })
            .collect();

        // build a list of colors for the query
        let query_colors: QueryColors = query
            .iter_mins()
            .filter_map(|hash| hash_to_color.get(hash))
            .map(|color| (*color, self.colors.indices(color)))
            .map(|(color, indices)| (color, indices.cloned().collect::<Vec<u32>>()))
            // CTB: could we add a 'from' to Datasets for this?
            .map(|(color, indices)| (color, Datasets::new(&indices)))
            .collect();

        //eprintln!("query_colors: {:?}", query_colors);

        CounterGather {
            counter,
            query_colors,
            hash_to_color,
        }
    }

    fn gather(
        &self,
        mut cg: CounterGather,
        threshold: usize,
        orig_query: &KmerMinHash,
        _selection: Option<Selection>,
    ) -> Result<Vec<GatherResult>> {
        let match_size = usize::MAX;
        let mut matches = vec![];

        let mut running_query = KmerMinHashBTree::from(orig_query.clone());
        while match_size > threshold && !cg.is_empty() {
            let next_match = cg.peek(threshold);
            if next_match.is_none() {
                break;
            }
            let (dataset_id, match_size) = next_match.unwrap();

            // eprintln!("dataset_id: {} {}", dataset_id, match_size);

            let query_mh = KmerMinHash::from(running_query.clone());
            let result = self.linear.gather_round(
                dataset_id,
                match_size,
                &query_mh,
                matches.len(),
                orig_query,
            )?;
            if let Some(Sketch::MinHash(match_mh)) =
                result.match_.select_sketch(self.linear.template())
            {
                let (matched_hashes, _) = match_mh.intersection(&query_mh)?;
                let mut isect_mh = match_mh.clone();
                isect_mh.clear();
                isect_mh.add_many(&matched_hashes)?;

                cg.consume(&isect_mh);
                matches.push(result);
                running_query.remove_many(isect_mh.iter_mins().copied())?;
            } else {
                unimplemented!()
            }
        }
        Ok(matches)
    }

    fn update(self, _collection: CollectionSet) -> Result<module::RevIndex> {
        Ok(module::RevIndex::Mem(self))
    }

    fn check(&self, _quick: bool) -> DbStats {
        unimplemented!()
    }

    fn compact(&self) {}

    fn flush(&self) -> Result<()> {
        Ok(())
    }

    fn collection(&self) -> &CollectionSet {
        self.linear.collection()
    }

    fn internalize_storage(&mut self) -> Result<()> {
        Ok(())
    }

    fn convert(&self, _output_db: module::RevIndex) -> Result<()> {
        todo!()
    }

    fn find_signatures(
        &self,
        mh: &KmerMinHash,
        threshold: f64,
        picklist: Option<DatasetPicklist>,
    ) -> Result<Vec<(f64, Signature, String)>> {
        let index_scaled = self.scaled();
        let query_scaled = mh.scaled();

        let query_mh = if query_scaled < index_scaled {
            mh.clone()
                .downsample_scaled(index_scaled)
                .expect("cannot downsample query")
        } else {
            mh.clone()
        };

        let threshold: usize = (threshold * (query_mh.size() as f64)) as _;

        let counter = self.counter_for_query(&query_mh, picklist);

        debug!(
            "number of matching signatures for hashes: {}",
            counter.len()
        );

        let mut results = vec![];
        for (dataset_id, size) in counter.most_common() {
            if size < threshold {
                break;
            };

            let match_sig = self.linear.sig_for_dataset(dataset_id)?;
            let match_path = self.location();

            let match_mh = match match_sig.select_sketch(self.linear.template()) {
                Some(Sketch::MinHash(mh)) => mh,
                _ => unimplemented!(),
            };

            if size >= threshold {
                let score = query_mh
                    .jaccard(match_mh)
                    .expect("cannot calculate Jaccard");

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
}

impl Index<'_> for MemRevIndex {
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
    use crate::Result;

    #[test]
    fn mem_revindex_new() -> Result<()> {
        let selection = Selection::builder().ksize(31).scaled(10000).build();
        let search_sigs = [
            "../../tests/test-data/gather/GCF_000006945.2_ASM694v2_genomic.fna.gz.sig".into(),
            "../../tests/test-data/gather/GCF_000007545.1_ASM754v1_genomic.fna.gz.sig".into(),
        ];
        let index = MemRevIndex::new(&search_sigs, &selection, 0, None)?;

        let index = match index {
            module::RevIndex::Mem(idx) => idx,
            _ => unimplemented!(),
        };

        assert_eq!(index.colors.len(), 3);

        Ok(())
    }

    #[test]
    fn mem_revindex_many() -> Result<()> {
        let selection = Selection::builder().ksize(31).scaled(10000).build();
        let search_sigs = [
            "../../tests/test-data/gather/GCF_000006945.2_ASM694v2_genomic.fna.gz.sig".into(),
            "../../tests/test-data/gather/GCF_000007545.1_ASM754v1_genomic.fna.gz.sig".into(),
            "../../tests/test-data/gather/GCF_000008105.1_ASM810v1_genomic.fna.gz.sig".into(),
        ];

        let index = MemRevIndex::new(&search_sigs, &selection, 0, None)?;
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
        let index = match index {
            module::RevIndex::Mem(idx) => idx,
            _ => unimplemented!(),
        };

        assert_eq!(index.colors.len(), 7);

        Ok(())
    }

    #[test]
    fn mem_revindex_from_sigs() -> Result<()> {
        let selection = Selection::builder().ksize(31).scaled(10000).build();
        let search_sigs: Vec<Signature> = [
            "../../tests/test-data/gather/GCF_000006945.2_ASM694v2_genomic.fna.gz.sig",
            "../../tests/test-data/gather/GCF_000007545.1_ASM754v1_genomic.fna.gz.sig",
            "../../tests/test-data/gather/GCF_000008105.1_ASM810v1_genomic.fna.gz.sig",
        ]
        .into_iter()
        .map(|path| Signature::from_path(path).unwrap().swap_remove(0))
        .collect();

        let index = MemRevIndex::new_with_sigs(search_sigs, &selection, 0, None)?;
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
        let index = match index {
            module::RevIndex::Mem(idx) => idx,
            _ => unimplemented!(),
        };

        assert_eq!(index.colors.len(), 7);

        Ok(())
    }

    #[test]
    fn mem_revindex_from_zipstorage() -> Result<()> {
        let selection = Selection::builder()
            .ksize(19)
            .scaled(100)
            .moltype(crate::encodings::HashFunctions::Murmur64Protein)
            .build();
        let index = MemRevIndex::from_zipfile(
            "../../tests/test-data/prot/protein.zip",
            &selection,
            0,
            None,
        )
        .expect("error building from ziptorage");

        let index = match index {
            module::RevIndex::Mem(idx) => idx,
            _ => unimplemented!(),
        };

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

        let counter_rev = index.counter_for_query(&query_mh, None);
        let counter_lin = index.linear.counter_for_query(&query_mh);

        let results_rev = index.search(counter_rev, false, 0).unwrap();
        let results_linear = index.linear.search(counter_lin, false, 0).unwrap();
        assert_eq!(results_rev, results_linear);

        let counter_rev = index.prepare_gather_counters(&query_mh, None);
        let counter_lin = index.linear.counter_for_query(&query_mh);

        let results_rev = index.gather(counter_rev, 0, &query_mh, None).unwrap();
        let results_linear = index.linear.gather(counter_lin, 0, &query_mh).unwrap();
        assert_eq!(results_rev.len(), 1);
        assert_eq!(results_rev, results_linear);

        Ok(())
    }

    #[test]
    fn mem_revindex_test_gather_2() -> Result<()> {
        let selection = Selection::builder().ksize(31).scaled(100000).build();
        let search_sigs: Vec<Signature> = [
            "../../tests/test-data/2.fa.sig",
            "../../tests/test-data/47.fa.sig",
        ]
        .into_iter()
        .map(|path| Signature::from_path(path).unwrap().swap_remove(0))
        .collect();

        let query_sig = Signature::from_path("../../tests/test-data/63.fa.sig")
            .expect("error processing query")
            .swap_remove(0)
            .select(&selection)
            .expect("error getting compatible sig");

        let query_mh = prepare_query(query_sig, &selection).expect("can't get compatible MinHash");

        let index = MemRevIndex::new_with_sigs(search_sigs, &selection, 0, None)?;

        let index = match index {
            module::RevIndex::Mem(idx) => idx,
            _ => unimplemented!(),
        };

        let gather_cg = index.prepare_gather_counters(&query_mh, None);
        // eprintln!("gather_cg: {:?}", gather_cg);
        let results = index.gather(gather_cg, 0, &query_mh, None).unwrap();

        assert_eq!(results.len(), 1);

        Ok(())
    }

    #[test]
    fn mem_revindex_test_gather_3() -> Result<()> {
        let selection = Selection::builder().ksize(31).scaled(100000).build();
        let search_sigs: Vec<Signature> = [
            "../../tests/test-data/2.fa.sig",
            "../../tests/test-data/47.fa.sig",
            "../../tests/test-data/63.fa.sig",
        ]
        .into_iter()
        .map(|path| Signature::from_path(path).unwrap().swap_remove(0))
        .collect();

        let query_sig = Signature::from_path("../../tests/test-data/SRR606249.sig.gz")
            .expect("error processing query")
            .swap_remove(0)
            .select(&selection)
            .expect("error getting compatible sig");

        let query_mh = prepare_query(query_sig, &selection).expect("can't get compatible MinHash");

        let index = MemRevIndex::new_with_sigs(search_sigs, &selection, 0, None)?;

        let index = match index {
            module::RevIndex::Mem(idx) => idx,
            _ => unimplemented!(),
        };

        // run the CounterGather-style gather:
        let gather_cg = index.prepare_gather_counters(&query_mh, None);
        // eprintln!("gather_cg: {:?}", gather_cg);
        let results = index.gather(gather_cg, 0, &query_mh, None).unwrap();
        assert_eq!(results.len(), 3);

        // compare to linear gather.
        let counter_lin = index.linear.counter_for_query(&query_mh);
        let results_linear = index.linear.gather(counter_lin, 0, &query_mh).unwrap();
        assert_eq!(results_linear.len(), 3);
        assert_eq!(results, results_linear);

        Ok(())
    }

    #[test]
    fn mem_revindex_load_and_gather_2() -> Result<()> {
        let selection = Selection::builder().ksize(21).scaled(10000).build();

        let mut basedir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        basedir.push("../../tests/test-data/gather/");

        let against = vec![
            "GCF_000006945.2_ASM694v2_genomic.fna.gz.sig",
            "GCF_000007545.1_ASM754v1_genomic.fna.gz.sig",
            "GCF_000008105.1_ASM810v1_genomic.fna.gz.sig",
            "GCF_000008545.1_ASM854v1_genomic.fna.gz.sig",
            "GCF_000009085.1_ASM908v1_genomic.fna.gz.sig",
            "GCF_000009505.1_ASM950v1_genomic.fna.gz.sig",
            "GCF_000009525.1_ASM952v1_genomic.fna.gz.sig",
            "GCF_000011885.1_ASM1188v1_genomic.fna.gz.sig",
            "GCF_000016045.1_ASM1604v1_genomic.fna.gz.sig",
            "GCF_000016785.1_ASM1678v1_genomic.fna.gz.sig",
            "GCF_000018945.1_ASM1894v1_genomic.fna.gz.sig",
            "GCF_000195995.1_ASM19599v1_genomic.fna.gz.sig",
        ];
        let against: Vec<_> = against
            .into_iter()
            .map(|sig| {
                let mut path = basedir.clone();
                path.push(sig);
                Signature::from_path(path).unwrap().swap_remove(0)
            })
            .collect();

        // build 'against' sketches into a revindex
        let index = MemRevIndex::new_with_sigs(against, &selection, 0, None)?;

        let mut query = None;
        let mut query_filename = basedir.clone();
        query_filename.push("combined.sig");
        let query_sig = Signature::from_path(query_filename)?
            .swap_remove(0)
            .select(&selection)?;

        if let Some(q) = prepare_query(query_sig, &selection) {
            query = Some(q);
        }
        let query = query.unwrap();

        let cg = index.prepare_gather_counters(&query, None);

        let matches = index.gather(
            cg,
            5, // 50kb threshold
            &query,
            Some(selection),
        )?;

        // should be 11, based on test_gather_metagenome_num_results
        assert_eq!(matches.len(), 11);

        fn round5(a: f64) -> f64 {
            (a * 1e5).round() / 1e5
        }

        let match_ = &matches[0];
        let names: Vec<&str> = match_.name().split(' ').take(1).collect();
        assert_eq!(names[0], "NC_003198.1");
        assert_eq!(match_.f_match(), 1.0);
        assert_eq!(round5(match_.f_unique_to_query()), round5(0.33219645));

        let match_ = &matches[1];
        let names: Vec<&str> = match_.name().split(' ').take(1).collect();
        assert_eq!(names[0], "NC_000853.1");
        assert_eq!(match_.f_match(), 1.0);
        assert_eq!(round5(match_.f_unique_to_query()), round5(0.13096862));

        let match_ = &matches[2];
        let names: Vec<&str> = match_.name().split(' ').take(1).collect();
        assert_eq!(names[0], "NC_011978.1");
        assert_eq!(match_.f_match(), 0.898936170212766);
        assert_eq!(round5(match_.f_unique_to_query()), round5(0.115279));

        let match_ = &matches[3];
        let names: Vec<&str> = match_.name().split(' ').take(1).collect();
        assert_eq!(names[0], "NC_002163.1");
        assert_eq!(match_.f_match(), 1.0);
        assert_eq!(round5(match_.f_unique_to_query()), round5(0.10709413));

        let match_ = &matches[4];
        let names: Vec<&str> = match_.name().split(' ').take(1).collect();
        assert_eq!(names[0], "NC_003197.2");
        assert_eq!(round5(match_.f_match()), round5(0.31340206));
        assert_eq!(round5(match_.f_unique_to_query()), round5(0.103683));

        let match_ = &matches[5];
        dbg!(match_);
        let names: Vec<&str> = match_.name().split(' ').take(1).collect();
        assert_eq!(names[0], "NC_009486.1");
        assert_eq!(round5(match_.f_match()), round5(0.4842105));
        assert_eq!(round5(match_.f_unique_to_query()), round5(0.0627557));

        let match_ = &matches[6];
        dbg!(match_);
        let names: Vec<&str> = match_.name().split(' ').take(1).collect();
        assert_eq!(names[0], "NC_006905.1");
        assert_eq!(round5(match_.f_match()), round5(0.161016949152542));
        assert_eq!(
            round5(match_.f_unique_to_query()),
            round5(0.0518417462482947)
        );

        let match_ = &matches[7];
        dbg!(match_);
        let names: Vec<&str> = match_.name().split(' ').take(1).collect();
        assert_eq!(names[0], "NC_011080.1");
        assert_eq!(round5(match_.f_match()), round5(0.125799573560768));
        assert_eq!(
            round5(match_.f_unique_to_query()),
            round5(0.04024556616643930)
        );

        let match_ = &matches[8];
        dbg!(match_);
        let names: Vec<&str> = match_.name().split(' ').take(1).collect();
        assert_eq!(names[0], "NC_011274.1");
        assert_eq!(round5(match_.f_match()), round5(0.0919037199124727));
        assert_eq!(
            round5(match_.f_unique_to_query()),
            round5(0.0286493860845839)
        );

        let match_ = &matches[9];
        dbg!(match_);
        let names: Vec<&str> = match_.name().split(' ').take(1).collect();
        assert_eq!(names[0], "NC_006511.1");
        assert_eq!(round5(match_.f_match()), round5(0.0725995316159251));
        assert_eq!(
            round5(match_.f_unique_to_query()),
            round5(0.021145975443383400)
        );

        let match_ = &matches[10];
        dbg!(match_);
        let names: Vec<&str> = match_.name().split(' ').take(1).collect();
        assert_eq!(names[0], "NC_011294.1");
        assert_eq!(round5(match_.f_match()), round5(0.0148619957537155));
        assert_eq!(
            round5(match_.f_unique_to_query()),
            round5(0.0047748976807639800)
        );

        Ok(())
    }

    #[test]
    fn revindex_load_and_test_counter_gather() -> Result<()> {
        let mut basedir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        basedir.push("../../tests/test-data/gather/");

        let against = vec![
            "GCF_000006945.2_ASM694v2_genomic.fna.gz.sig",
            "GCF_000007545.1_ASM754v1_genomic.fna.gz.sig",
            "GCF_000008105.1_ASM810v1_genomic.fna.gz.sig",
            "GCF_000008545.1_ASM854v1_genomic.fna.gz.sig",
            "GCF_000009085.1_ASM908v1_genomic.fna.gz.sig",
            "GCF_000009505.1_ASM950v1_genomic.fna.gz.sig",
            "GCF_000009525.1_ASM952v1_genomic.fna.gz.sig",
            "GCF_000011885.1_ASM1188v1_genomic.fna.gz.sig",
            "GCF_000016045.1_ASM1604v1_genomic.fna.gz.sig",
            "GCF_000016785.1_ASM1678v1_genomic.fna.gz.sig",
            "GCF_000018945.1_ASM1894v1_genomic.fna.gz.sig",
            "GCF_000195995.1_ASM19599v1_genomic.fna.gz.sig",
        ];
        let against: Vec<PathBuf> = against
            .iter()
            .map(|sig| {
                let mut filename = basedir.clone();
                filename.push(sig);
                filename.into()
            })
            .collect();

        // build 'against' sketches into a revindex
        let selection = Selection::builder().ksize(21).scaled(10000).build();

        let index = MemRevIndex::new(&against[..], &selection, 0, None)?;

        let mut query = None;
        let mut query_filename = basedir.clone();
        query_filename.push("combined.sig");
        let query_sig = Signature::from_path(query_filename)?
            .swap_remove(0)
            .select(&selection)?;

        if let Some(q) = prepare_query(query_sig, &selection) {
            query = Some(q);
        }
        let query = query.unwrap();

        let cg = index.prepare_gather_counters(&query, None);

        let idxlist = cg.dataset_ids();
        assert_eq!(idxlist.len(), 12);

        let found_mh = cg.found_hashes(&query);
        assert_eq!(found_mh.size(), 1466);

        Ok(())
    }

    #[test]
    fn revindex_load_and_gather_picklist() -> Result<()> {
        let mut basedir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        basedir.push("../../tests/test-data/gather/");

        let against = vec![
            "GCF_000006945.2_ASM694v2_genomic.fna.gz.sig",
            "GCF_000007545.1_ASM754v1_genomic.fna.gz.sig",
            "GCF_000008105.1_ASM810v1_genomic.fna.gz.sig",
            "GCF_000008545.1_ASM854v1_genomic.fna.gz.sig",
            "GCF_000009085.1_ASM908v1_genomic.fna.gz.sig",
            "GCF_000009505.1_ASM950v1_genomic.fna.gz.sig",
            "GCF_000009525.1_ASM952v1_genomic.fna.gz.sig",
            "GCF_000011885.1_ASM1188v1_genomic.fna.gz.sig",
            "GCF_000016045.1_ASM1604v1_genomic.fna.gz.sig",
            "GCF_000016785.1_ASM1678v1_genomic.fna.gz.sig",
            "GCF_000018945.1_ASM1894v1_genomic.fna.gz.sig",
            "GCF_000195995.1_ASM19599v1_genomic.fna.gz.sig",
        ];
        let against: Vec<PathBuf> = against
            .iter()
            .map(|sig| {
                let mut filename = basedir.clone();
                filename.push(sig);
                filename.into()
            })
            .collect();

        // build 'against' sketches into a revindex
        let selection = Selection::builder().ksize(21).scaled(10000).build();

        let index = MemRevIndex::new(&against[..], &selection, 0, None)?;

        let mut query = None;
        let mut query_filename = basedir.clone();
        query_filename.push("combined.sig");
        let query_sig = Signature::from_path(query_filename)?
            .swap_remove(0)
            .select(&selection)?;

        if let Some(q) = prepare_query(query_sig, &selection) {
            query = Some(q);
        }
        let query = query.unwrap();

        // build a picklist with only one match
        let pl = DatasetPicklist {
            dataset_ids: vec![0].into_iter().collect(),
        };

        let cg = index.prepare_gather_counters(&query, Some(pl.clone()));

        let matches = index.gather(
            cg,
            5, // 50kb threshold
            &query,
            Some(selection),
        )?;

        // should be 1, b/c of picklist.
        assert_eq!(matches.len(), 1);

        // also do a basic test of containment with picklists -
        let counter = index.counter_for_query(&query, Some(pl.clone()));
        let matches = index.matches_from_counter(counter, 0);
        assert_eq!(matches, [("NC_003197.2 Salmonella enterica subsp. enterica serovar Typhimurium str. LT2, complete genome".into(), 485)]);

        let counter = index.counter_for_query(&query, Some(pl));
        let records = index.records_from_counter(counter, 0);
        assert_eq!(records.len(), 1);

        Ok(())
    }

    #[test]
    fn mem_revindex_find_signatures() -> Result<()> {
        let selection = Selection::builder().ksize(31).scaled(100000).build();
        let search_sigs: Vec<Signature> = [
            "../../tests/test-data/2.fa.sig",
            "../../tests/test-data/47.fa.sig",
            "../../tests/test-data/63.fa.sig",
        ]
        .into_iter()
        .map(|path| Signature::from_path(path).unwrap().swap_remove(0))
        .collect();

        let index = MemRevIndex::new_with_sigs(search_sigs, &selection, 0, None)?;

        let query_sig = Signature::from_path("../../tests/test-data/63.fa.sig")
            .expect("error processing query")
            .swap_remove(0)
            .select(&selection)
            .expect("error getting compatible sig");

        let query_mh = prepare_query(query_sig, &selection).expect("can't get compatible MinHash");

        let results = index.find_signatures(&query_mh, 0.0, None)?;
        assert_eq!(results.len(), 2);

        let results = index.find_signatures(&query_mh, 1.0, None)?;
        assert_eq!(results.len(), 1);

        // build a picklist with only one Idx (2.fa) => no match
        let pl = DatasetPicklist {
            dataset_ids: vec![0].into_iter().collect(),
        };
        let results = index.find_signatures(&query_mh, 0.0, Some(pl))?;
        assert_eq!(results.len(), 0);

        // build a picklist with only one Idx (47.fa) => one match
        let pl = DatasetPicklist {
            dataset_ids: vec![1].into_iter().collect(),
        };
        let results = index.find_signatures(&query_mh, 0.0, Some(pl))?;
        assert_eq!(results.len(), 1);

        Ok(())
    }
}
