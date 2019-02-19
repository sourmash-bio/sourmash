use std::fs::File;
use std::io::{BufReader, Read};
use std::mem;
use std::path::Path;
use std::path::PathBuf;
use std::rc::Rc;

use failure::Error;
use lazy_init::Lazy;
use serde_derive::{Deserialize, Serialize};
use typed_builder::TypedBuilder;

use crate::index::storage::{FSStorage, ReadData, Storage, StorageInfo, ToWriter};
use crate::index::{Comparable, Dataset, DatasetInfo, Index};

#[derive(TypedBuilder)]
pub struct LinearIndex<L> {
    #[builder(default)]
    storage: Option<Rc<dyn Storage>>,

    #[builder(default)]
    pub(crate) datasets: Vec<L>,
}

#[derive(Serialize, Deserialize)]
struct LinearInfo<L> {
    version: u32,
    storage: StorageInfo,
    leaves: Vec<L>,
}

impl<L> Index for LinearIndex<L>
where
    L: Clone + Comparable<L>,
{
    type Item = L;

    fn find<F>(
        &self,
        search_fn: F,
        sig: &Self::Item,
        threshold: f64,
    ) -> Result<Vec<&Self::Item>, Error>
    where
        F: Fn(&dyn Comparable<Self::Item>, &Self::Item, f64) -> bool,
    {
        Ok(self
            .datasets
            .iter()
            .flat_map(|node| {
                if search_fn(node, sig, threshold) {
                    Some(node)
                } else {
                    None
                }
            })
            .collect())
    }

    fn insert(&mut self, node: &L) -> Result<(), Error> {
        self.datasets.push(node.clone());
        Ok(())
    }

    fn save<P: AsRef<Path>>(&self, _path: P) -> Result<(), Error> {
        /*
        let file = File::create(path)?;
        match serde_json::to_writer(file, &self) {
            Ok(_) => Ok(()),
            Err(_) => Err(SourmashError::SerdeError.into()),
        }
        */
        unimplemented!()
    }

    fn load<P: AsRef<Path>>(_path: P) -> Result<(), Error> {
        unimplemented!()
    }

    fn datasets(&self) -> Vec<Self::Item> {
        self.datasets.to_vec()
    }
}

impl<L> LinearIndex<Dataset<L>>
where
    L: std::marker::Sync + ToWriter,
    Dataset<L>: ReadData<L>,
{
    pub fn save_file<P: AsRef<Path>>(
        &mut self,
        path: P,
        storage: Option<Rc<dyn Storage>>,
    ) -> Result<(), Error> {
        let ref_path = path.as_ref();
        let mut basename = ref_path.file_name().unwrap().to_str().unwrap().to_owned();
        if basename.ends_with(".sbt.json") {
            basename = basename.replace(".sbt.json", "");
        }
        let location = ref_path.parent().unwrap();

        let storage = match storage {
            Some(s) => s,
            None => {
                let subdir = format!(".linear.{}", basename);
                Rc::new(FSStorage::new(location.to_str().unwrap(), &subdir))
            }
        };

        let args = storage.args();
        let storage_info = StorageInfo {
            backend: "FSStorage".into(),
            args,
        };

        let info: LinearInfo<DatasetInfo> = LinearInfo {
            storage: storage_info,
            version: 5,
            leaves: self
                .datasets
                .iter_mut()
                .map(|l| {
                    // Trigger data loading
                    let _: &L = (*l).data().unwrap();

                    // set storage to new one
                    mem::replace(&mut l.storage, Some(Rc::clone(&storage)));

                    let filename = (*l).save(&l.filename).unwrap();
                    let new_node = DatasetInfo {
                        filename: filename,
                        name: l.name.clone(),
                        metadata: l.metadata.clone(),
                    };
                    new_node
                })
                .collect(),
        };

        let file = File::create(path)?;
        serde_json::to_writer(file, &info)?;

        Ok(())
    }

    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<LinearIndex<Dataset<L>>, Error> {
        let file = File::open(&path)?;
        let mut reader = BufReader::new(file);

        // TODO: match with available Storage while we don't
        // add a function to build a Storage from a StorageInfo
        let mut basepath = PathBuf::new();
        basepath.push(path);
        basepath.canonicalize()?;

        let linear =
            LinearIndex::<Dataset<L>>::from_reader(&mut reader, &basepath.parent().unwrap())?;
        Ok(linear)
    }

    pub fn from_reader<R, P>(rdr: &mut R, path: P) -> Result<LinearIndex<Dataset<L>>, Error>
    where
        R: Read,
        P: AsRef<Path>,
    {
        // TODO: check https://serde.rs/enum-representations.html for a
        // solution for loading v4 and v5
        let linear: LinearInfo<DatasetInfo> = serde_json::from_reader(rdr)?;

        // TODO: support other storages
        let mut st: FSStorage = (&linear.storage.args).into();
        st.set_base(path.as_ref().to_str().unwrap());
        let storage: Rc<dyn Storage> = Rc::new(st);

        Ok(LinearIndex {
            storage: Some(Rc::clone(&storage)),
            datasets: linear
                .leaves
                .into_iter()
                .map(|l| {
                    let new_node = Dataset {
                        filename: l.filename,
                        name: l.name,
                        metadata: l.metadata,
                        storage: Some(Rc::clone(&storage)),
                        data: Rc::new(Lazy::new()),
                    };
                    new_node
                })
                .collect(),
        })
    }

    pub fn storage(&self) -> Option<Rc<dyn Storage>> {
        self.storage.clone()
    }
}
