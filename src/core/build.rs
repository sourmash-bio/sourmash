use std::env;

fn main() {
    let crate_dir = env::var("CARGO_MANIFEST_DIR").unwrap();
    copy_c_bindings(&crate_dir);
}

#[cfg(not(feature = "maturin"))]
fn copy_c_bindings(_crate_dir: &str) {}

#[cfg(feature = "maturin")]
fn copy_c_bindings(crate_dir: &str) {
    use std::path::{Path, PathBuf};

    fn find_root_dir(crate_dir: &str) -> &Path {
        let root_dir = Path::new(crate_dir);

        if root_dir.join("pyproject.toml").is_file() {
            return root_dir;
        }

        let root_dir = Path::new(crate_dir).parent().unwrap().parent().unwrap();
        if root_dir.join("pyproject.toml").is_file() {
            return root_dir;
        }

        panic!("Couldn't find pyproject.toml to determine root dir");
    }

    fn find_target_dir(out_dir: &str) -> PathBuf {
        use std::ffi::OsStr;

        let mut components = Path::new(out_dir).iter();

        while let Some(dir) = components.next_back() {
            if dir == OsStr::new("target") {
                break;
            }
        }
        let mut dir: PathBuf = components.collect();

        if dir.as_os_str().is_empty() {
            panic!("Couldn't find target dir based on OUT_DIR");
        } else {
            dir.push("target");
            dir
        }
    }

    let root_dir = find_root_dir(crate_dir);
    let header_path = root_dir.join("include").join("sourmash.h");
    let header = std::fs::read_to_string(header_path).expect("error reading header");

    // strip directives, not supported by the cffi C parser
    let new_header: String = header
        .lines()
        .filter_map(|s| {
            if s.starts_with("#") {
                None
            } else {
                Some({
                    let mut s = s.to_owned();
                    s.push_str("\n");
                    s
                })
            }
        })
        .collect();

    let out_dir = env::var("OUT_DIR").unwrap();
    let target_dir = find_target_dir(&out_dir);
    std::fs::create_dir_all(&target_dir).expect("error creating target dir");
    let out_path = target_dir.join("header.h");
    std::fs::write(out_path, &new_header).expect("error writing header");
}
