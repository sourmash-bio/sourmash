use sourmash::index::sbt::{Node, SBT};
use sourmash::signature::Signature;

use crate::sketch::ukhs::FlatUKHS;

pub mod ukhs;

pub type UKHSTree = SBT<Node<FlatUKHS>, Signature>;
