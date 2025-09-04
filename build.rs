use std::{env, fs, path::PathBuf};

fn main() {
    // Compile bcg729 C sources and link them statically
    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let c_src_dir = manifest_dir.join("src/bcg729/src");
    let c_include_dir = manifest_dir.join("src/bcg729/include");

    let mut build = cc::Build::new();
    build.include(&c_include_dir);
    build.flag_if_supported("-std=c99");
    build.flag_if_supported("-fvisibility=hidden");
    build.flag_if_supported("-fPIC");
    build.define("BCG729_STATIC", None);

    let entries = fs::read_dir(&c_src_dir).expect("failed to read bcg729/src directory");
    let mut found = false;
    for entry in entries {
        let path = entry.unwrap().path();
        if let Some(ext) = path.extension() {
            if ext == "c" {
                // Skip files that are not meant to be compiled standalone (none in upstream)
                build.file(&path);
                found = true;
            }
        }
    }
    assert!(found, "No C sources found in {}", c_src_dir.display());

    build.compile("bcg729");
}
