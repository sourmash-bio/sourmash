#![cfg(all(target_arch = "wasm32", target_os = "unknown"))]

use wasm_bindgen_test::wasm_bindgen_test_configure;

wasm_bindgen_test_configure!(run_in_shared_worker);
