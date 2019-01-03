use crate::index::Comparable;

pub fn search_minhashes<L>(node: &dyn Comparable<L>, query: &L, threshold: f64) -> bool {
    node.similarity(query) > threshold
}

pub fn search_minhashes_containment<L>(
    node: &dyn Comparable<L>,
    query: &L,
    threshold: f64,
) -> bool {
    node.containment(query) > threshold
}
