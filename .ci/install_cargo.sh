#! /bin/sh
curl https://sh.rustup.rs -sSf | sh -s -- -y --default-toolchain=stable
rustup show
export PATH="$HOME/.cargo/bin:$PATH"
rustc -V
rustup target add aarch64-apple-darwin

# update crates.io index without updating Cargo.lock
export CARGO_NET_GIT_FETCH_WITH_CLI=true
cargo update --dry-run
