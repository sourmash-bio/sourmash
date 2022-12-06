{

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
    utils.url = "github:numtide/flake-utils";

    rust-overlay = {
      url = "github:oxalica/rust-overlay";
      inputs = {
        nixpkgs.follows = "nixpkgs";
        flake-utils.follows = "utils";
      };
    };

    naersk = {
      url = "github:nix-community/naersk";
      inputs = {
        nixpkgs.follows = "nixpkgs";
        flake-utils.follows = "utils";
      };
    };
  };

  outputs = { self, nixpkgs, naersk, rust-overlay, utils }:
    utils.lib.eachDefaultSystem (system:
      let
        overlays = [ (import rust-overlay) ];
        pkgs = import nixpkgs {
          inherit system overlays;
          config.packageOverrides = pkgs:
            {
              maturin = pkgs.rustPlatform.buildRustPackage rec {
                pname = "maturin";
                version = "0.14.3";
                src = pkgs.fetchFromGitHub {
                  owner = "PyO3";
                  repo = "maturin";
                  rev = "v0.14.3";
                  hash = "sha256-n+z05My49LNl5nUaCMw90nLk6EiiAWB/HUCZANvrVuY=";
                };
                cargoHash = "sha256-xSUdFzYCbsdyg2uowMRmFSJL1s76KSU9ZJcPYwZsN3Q=";
                nativeBuildInputs = [ pkgs.pkg-config ];
#                buildInputs = lib.optionals stdenv.isLinux [ dbus ]
#                  ++ lib.optionals stdenv.isDarwin [ Security libiconv ];
                # Requires network access, fails in sandbox.
                doCheck = false;
#                passthru.tests.pyo3 = callPackage ./pyo3-test { };
                /*
                meta = with lib; {
                  description = "Build and publish Rust crates Python packages";
                  longDescription = ''
                    Build and publish Rust crates with PyO3, rust-cpython, and
                    cffi bindings as well as Rust binaries as Python packages.
                    This project is meant as a zero-configuration replacement for
                    setuptools-rust and Milksnake. It supports building wheels for
                    Python and can upload them to PyPI.
                  '';
                  homepage = "https://github.com/PyO3/maturin";
                  license = licenses.asl20;
                  maintainers = [ ];
                };
                */
              }
              ;
            };
        };
        rustVersion = pkgs.rust-bin.stable.latest.default.override {
          #extensions = [ "rust-src" ];
          #targets = [ "x86_64-unknown-linux-musl" ];
          targets = [ "wasm32-wasi" "wasm32-unknown-unknown" "wasm32-unknown-emscripten" ];
        };
        rustPlatform = pkgs.makeRustPlatform {
          cargo = rustVersion;
          rustc = rustVersion;
        };
        naersk-lib = naersk.lib."${system}".override {
          cargo = rustPlatform.rust.cargo;
          rustc = rustPlatform.rust.rustc;
        };

        python = pkgs.python310Packages;

        screed = python.buildPythonPackage rec {
          pname = "screed";
          version = "1.1";
          #format = "pyproject";

          src = pkgs.fetchFromGitHub {
            owner = "dib-lab";
            repo = "screed";
            rev = "v1.1";
            hash = "sha256-g1FZJx94RGBPoTiLfwttdYqCJ02pxtOKK708WA63kHE=";
          };

          SETUPTOOLS_SCM_PRETEND_VERSION = "1.1";
          propagatedBuildInputs = with python; [setuptools bz2file setuptools_scm ];
          doCheck = false;
        };

      in

      with pkgs;
      {
        packages = {

          lib = naersk-lib.buildPackage {
            pname = "libsourmash";
            root = ./.;
            copyLibs = true;
          };

          sourmash = python.buildPythonPackage rec {
            pname = "sourmash";
            version = "4.6.1";
            format = "pyproject";
            #disabled = lib.pythonOlder "3.6";

            src = ./.;
            cargoDeps = rustPlatform.fetchCargoTarball {
              inherit src;
              name = "${pname}-${version}";
              hash = "sha256-Zne9ERU3zQLK1anTEIhUyNx92Gr4YMFCf3cJU6U+BO4=";
            };

            nativeBuildInputs = with rustPlatform; [ cargoSetupHook maturinBuildHook ];

            buildInputs = lib.optionals stdenv.isDarwin [ libiconv ];
            propagatedBuildInputs = with python; [ cffi deprecation cachetools bitstring numpy scipy matplotlib screed ];

            DYLD_LIBRARY_PATH = "${self.packages.${system}.lib}/lib";
            NO_BUILD = "1";
          };

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

        devShell = mkShell {
          nativeBuildInputs = [
            clang_13
          ];

          buildInputs = [
            rustPlatform.rust.cargo
            openssl
            pkgconfig

            git
            stdenv.cc.cc.lib
            (python310.withPackages (ps: with ps; [ virtualenv tox setuptools ]))
            (python311.withPackages (ps: with ps; [ virtualenv setuptools ]))
            (python39.withPackages (ps: with ps; [ virtualenv setuptools ]))
            (python38.withPackages (ps: with ps; [ virtualenv setuptools ]))

            rust-cbindgen
            maturin

            wasmtime
            wasm-pack
            nodejs-16_x

            py-spy
            heaptrack
            cargo-watch
            cargo-limit
            cargo-outdated
            cargo-udeps
            nixpkgs-fmt

            llvmPackages_13.libclang
            llvmPackages_13.libcxxClang
          ];

          BINDGEN_EXTRA_CLANG_ARGS = "-isystem ${llvmPackages_13.libclang.lib}/lib/clang/${lib.getVersion clang}/include";
          LIBCLANG_PATH = "${llvmPackages_13.libclang.lib}/lib";
          LD_LIBRARY_PATH = "${stdenv.cc.cc.lib}/lib64:$LD_LIBRARY_PATH";

          # workaround for https://github.com/NixOS/nixpkgs/blob/48dfc9fa97d762bce28cc8372a2dd3805d14c633/doc/languages-frameworks/python.section.md#python-setuppy-bdist_wheel-cannot-create-whl
          SOURCE_DATE_EPOCH = 315532800; # 1980
        };
      });
}
