#[test]
fn test_murmur() {
    assert_eq!(sourmash::_hash_murmur(b"ACG", 42), 1731421407650554201)
}
