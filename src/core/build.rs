use std::env;

fn main() {
    let crate_dir = env::var("CARGO_MANIFEST_DIR").unwrap();
    copy_c_bindings(&crate_dir);
}

#[cfg(not(feature = "maturin"))]
fn copy_c_bindings(_crate_dir: &str) {}

#[cfg(feature = "maturin")]
fn copy_c_bindings(crate_dir: &str) {
    use std::path::Path;

    let header_path = Path::new(crate_dir)
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("include")
        .join("sourmash.h");
    let header = std::fs::read_to_string(header_path).unwrap_or_else(|_| {
        let header_path = Path::new(crate_dir).join("include").join("sourmash.h");
        std::fs::read_to_string(header_path).expect("error reading header")
    });

    // strip directives, not supported by the cffi C parser
    let new_header: String = header
        .lines()
        .filter_map(|s| if s.starts_with("#") { None } else { Some(s) })
        .collect();

    std::fs::write(
        Path::new(crate_dir).join("target").join("header.h"),
        new_header,
    )
    .expect("error writing header");
}
