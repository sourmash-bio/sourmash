{

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
    utils.url = "github:numtide/flake-utils";

    rust-overlay = {
      url = "github:oxalica/rust-overlay";
      inputs = {
        nixpkgs.follows = "nixpkgs";
      };
    };
  };

  outputs =
    {
      self,
      nixpkgs,
      rust-overlay,
      utils,
    }:
    utils.lib.eachDefaultSystem (
      system:
      let
        overlays = [ (import rust-overlay) ];
        pkgs = import nixpkgs {
          inherit system overlays;
        };
        rustVersion = pkgs.rust-bin.stable.latest.default.override {
          #extensions = [ "rust-src" ];
          extensions = [ "llvm-tools-preview" ];
          #targets = [ "x86_64-unknown-linux-musl" ];
          targets = [
            "wasm32-unknown-unknown"
            "wasm32-unknown-emscripten"
          ];
        };
        rustPlatform = pkgs.makeRustPlatform {
          cargo = rustVersion;
          rustc = rustVersion;
        };

        inherit (pkgs) lib;

        python = pkgs.python314Packages;
        rocksdb = pkgs.rocksdb;

        commonArgs = {
          src = ./.;
          preConfigure = lib.optionalString pkgs.stdenv.isDarwin ''
            export MACOSX_DEPLOYMENT_TARGET=10.14
          '';

          nativeBuildInputs =
            with rustPlatform;
            [
              cargoSetupHook
              maturinBuildHook
              bindgenHook
            ]
            ++ [
              rocksdb
            ];

          env = {
            ROCKSDB_INCLUDE_DIR = "${rocksdb}/include";
            ROCKSDB_LIB_DIR = "${rocksdb}/lib";
          };
        };

      in

      with pkgs;
      {
        packages = {

          lib = rustPlatform.buildRustPackage (
            commonArgs
            // {
              name = "libsourmash";
              copyLibs = true;
              cargoLock.lockFile = ./Cargo.lock;
              nativeBuildInputs = with rustPlatform; [ bindgenHook ];
            }
          );

          sourmash = python.buildPythonPackage (
            commonArgs
            // rec {
              pname = "sourmash";
              version = "4.9.4";
              format = "pyproject";

              cargoDeps = rustPlatform.importCargoLock {
                lockFile = ./Cargo.lock;
              };

              propagatedBuildInputs = with python; [
                cffi
                deprecation
                cachetools
                bitstring
                numpy
                scipy
                matplotlib
                screed
              ];
            }
          );

          docker =
            let
              bin = self.defaultPackage.${system};
            in
            pkgs.dockerTools.buildLayeredImage {
              name = bin.pname;
              tag = bin.version;
              contents = [ bin ];

              config = {
                Cmd = [ "/bin/sourmash" ];
                WorkingDir = "/";
              };
            };
        };

        defaultPackage = self.packages.${system}.sourmash;

        devShells.default = pkgs.mkShell (
          commonArgs
          // {
            nativeBuildInputs = with rustPlatform; [ bindgenHook ];

            buildInputs = [
              rustVersion
              openssl
              pkg-config

              git
              pkgs.stdenv.cc.cc.lib
              (python314.withPackages (
                ps: with ps; [
                  virtualenv
                  tox
                  cffi
                ]
              ))
              (python313.withPackages (ps: with ps; [ virtualenv ]))
              (python312.withPackages (ps: with ps; [ virtualenv ]))
              (python311.withPackages (ps: with ps; [ virtualenv ]))

              #rust-cbindgen
              maturin

              wasmtime
              wasm-pack
              nodejs_20
              #emscripten

              #py-spy
              #heaptrack
              cargo-all-features
              cargo-watch
              cargo-limit
              cargo-outdated
              cargo-udeps
              cargo-deny
              cargo-nextest
              cargo-llvm-cov
              cargo-component
              cargo-codspeed
              #cargo-semver-checks
              nixpkgs-fmt
            ];

            shellHook = ''
              export MACOSX_DEPLOYMENT_TARGET=10.14
            '';

            # Needed for matplotlib
            LD_LIBRARY_PATH = lib.makeLibraryPath [ pkgs.stdenv.cc.cc.lib ];

            # workaround for https://github.com/NixOS/nixpkgs/blob/48dfc9fa97d762bce28cc8372a2dd3805d14c633/doc/languages-frameworks/python.section.md#python-setuppy-bdist_wheel-cannot-create-whl
            SOURCE_DATE_EPOCH = 315532800; # 1980

            # exporting to fix doc building errors in sphinx
            LC_ALL = "C.utf8";
          }
        );
      }
    );
}
