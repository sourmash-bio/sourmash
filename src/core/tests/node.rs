#![cfg(all(target_arch = "wasm32", target_os = "unknown"))]

use wasm_bindgen_test::*;

#[wasm_bindgen_test]
fn pass() {
    assert_eq!(1, 1);
}
