use getset::{CopyGetters, Getters, Setters};
use typed_builder::TypedBuilder;

use crate::encodings::HashFunctions;
use crate::manifest::Record;
use crate::Result;
use crate::sketch::Sketch;

#[derive(Default, Debug, TypedBuilder, Clone)]
pub struct Selection {
    #[builder(default, setter(strip_option))]
    ksize: Option<u32>,

    #[builder(default, setter(strip_option))]
    abund: Option<bool>,

    #[builder(default, setter(strip_option))]
    num: Option<u32>,

    #[builder(default, setter(strip_option))]
    scaled: Option<u32>,

    #[builder(default, setter(strip_option))]
    containment: Option<bool>,

    #[builder(default, setter(strip_option))]
    moltype: Option<HashFunctions>,

    #[builder(default, setter(strip_option))]
    picklist: Option<Picklist>,

    #[builder(default, setter(strip_option))]
    sketchtype: Option<Sketch>,
}

#[derive(Default, TypedBuilder, CopyGetters, Getters, Setters, Clone, Debug)]
pub struct Picklist {
    #[getset(get = "pub", set = "pub")]
    #[builder(default = "".into())]
    coltype: String,

    #[getset(get = "pub", set = "pub")]
    #[builder(default = "".into())]
    pickfile: String,

    #[getset(get = "pub", set = "pub")]
    #[builder(default = "".into())]
    column_name: String,

    #[getset(get = "pub", set = "pub")]
    #[builder]
    pickstyle: PickStyle,
}

#[derive(Clone, Default, Debug)]
#[repr(u32)]
pub enum PickStyle {
    #[default]
    Include = 1,
    Exclude = 2,
}

pub trait Select {
    fn select(self, selection: &Selection) -> Result<Self>
    where
        Self: Sized;
}

impl Selection {
    pub fn ksize(&self) -> Option<u32> {
        self.ksize
    }

    pub fn set_ksize(&mut self, ksize: u32) {
        self.ksize = Some(ksize);
    }

    pub fn abund(&self) -> Option<bool> {
        self.abund
    }

    pub fn set_abund(&mut self, value: bool) {
        self.abund = Some(value);
    }

    pub fn num(&self) -> Option<u32> {
        self.num
    }

    pub fn set_num(&mut self, num: u32) {
        self.num = Some(num);
    }

    pub fn scaled(&self) -> Option<u32> {
        self.scaled
    }

    pub fn set_scaled(&mut self, scaled: u32) {
        self.scaled = Some(scaled);
    }

    pub fn containment(&self) -> Option<bool> {
        self.containment
    }

    pub fn set_containment(&mut self, containment: bool) {
        self.containment = Some(containment);
    }

    pub fn moltype(&self) -> Option<HashFunctions> {
        self.moltype.clone()
    }

    pub fn set_moltype(&mut self, value: HashFunctions) {
        self.moltype = Some(value);
    }

    pub fn picklist(&self) -> Option<Picklist> {
        self.picklist.clone()
    }

    pub fn set_picklist(&mut self, value: Picklist) {
        self.picklist = Some(value);
    }

    pub fn set_sketchtype(&mut self, value: Sketch) {
        self.sketchtype = Some(value);
    }

    pub fn from_record(row: &Record) -> Result<Self> {
        Ok(Self {
            ksize: Some(*row.ksize()),
            abund: Some(*row.with_abundance()),
            moltype: Some(row.moltype()),
            num: None,
            scaled: None,
            containment: None,
            picklist: None,
            sketchtype: None,
        })
    }
}
