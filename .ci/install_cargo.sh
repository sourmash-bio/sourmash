#! /bin/sh
curl https://sh.rustup.rs -sSf | sh -s -- -y --default-toolchain stable
export PATH="$HOME/.cargo/bin:$PATH"
rustc -V
rustup target add aarch64-apple-darwin
