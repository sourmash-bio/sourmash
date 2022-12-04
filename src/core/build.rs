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

    let out_dir = env::var("OUT_DIR").unwrap();

    // Hack to try to find header back in workspace root.
    // Ideally wouldn't need to specify the crate Cargo.toml,
    // but maturin doesn't work well with workspaces yet.
    let header_path = Path::new(crate_dir)
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("include")
        .join("sourmash.h");
    let header = std::fs::read_to_string(header_path).unwrap_or_else(|_| {
        // Fallback to find header if workspace Cargo.toml is used.
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
        &new_header,
    )
    .unwrap_or_else(|_| {
        // Need this hack to support editable installations.
        // Ends up finding the `target` dir based on OUT_DIR.
        // Hack fixable with better cargo workspace support in maturin.
        let dir = Path::new(&out_dir)
            .parent()
            .unwrap()
            .parent()
            .unwrap()
            .parent()
            .unwrap()
            .parent()
            .unwrap();
        std::fs::create_dir_all(&dir).expect("error creating target dir");

        std::fs::write(dir.join("header.h"), new_header).expect("error writing header");
    });
}
