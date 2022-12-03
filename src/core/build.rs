use std::env;
use std::path::Path;

fn main() {
    let crate_dir = env::var("CARGO_MANIFEST_DIR").unwrap();
    generate_c_bindings(&crate_dir);
}

#[cfg(feature = "cbindgen")]
fn generate_c_bindings(crate_dir: &str) {
    cbindgen::generate(crate_dir)
        .expect("Unable to generate bindings")
        .write_to_file("");
    bindings.write_to_file(Path::new("target").join("header.h"));
}

#[cfg(not(feature = "cbindgen"))]
fn generate_c_bindings(crate_dir: &str) {
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
