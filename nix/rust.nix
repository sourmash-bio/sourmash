# nix/rust.nix
{ sources ? import ./sources.nix }:
let
  pkgs =
    import sources.nixpkgs { overlays = [ (import sources.rust-overlay) ]; };
  rustVersion = pkgs.rust-bin.stable.latest.rust.override {
    #extensions = [ "rust-src" ];
    #targets = [ "x86_64-unknown-linux-musl" ];
    targets = [ "wasm32-wasi" "wasm32-unknown-unknown" ];
  };
in
pkgs.makeRustPlatform {
  cargo = rustVersion;
  rustc = rustVersion;
}
