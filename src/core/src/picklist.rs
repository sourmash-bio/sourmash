use getset::{CopyGetters, Getters, Setters};
use typed_builder::TypedBuilder;

#[derive(Default, TypedBuilder, CopyGetters, Getters, Setters, Clone)]
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

#[derive(Default, Clone)]
#[repr(u32)]
pub enum PickStyle {
    #[default]
    Include = 1,
    Exclude = 2,
}
