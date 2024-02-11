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
  };

  outputs = { self, nixpkgs, rust-overlay, utils }:
    utils.lib.eachDefaultSystem (system:
      let
        overlays = [ (import rust-overlay) ];
        pkgs = import nixpkgs {
          inherit system overlays;
        };
        rustVersion = pkgs.rust-bin.stable.latest.default.override {
          extensions = [ "rust-src" "llvm-tools-preview" ];
          #targets = [ "x86_64-unknown-linux-musl" ];
          targets = [ "wasm32-wasi" "wasm32-unknown-unknown" "wasm32-unknown-emscripten" ];
        };
        rustPlatform = pkgs.makeRustPlatform {
          cargo = rustVersion;
          rustc = rustVersion;
        };

        inherit (pkgs) lib;

        python = pkgs.python311Packages;

        stdenv = if pkgs.stdenv.isDarwin then pkgs.overrideSDK pkgs.stdenv "11.0" else pkgs.stdenv;

        commonArgs = {
          src = ./.;
          stdenv = stdenv;
          preConfigure = lib.optionalString stdenv.isDarwin ''
            export MACOSX_DEPLOYMENT_TARGET=10.14
          '';

          # Work around https://github.com/NixOS/nixpkgs/issues/166205.
          env = lib.optionalAttrs stdenv.cc.isClang {
            NIX_LDFLAGS = "-l${stdenv.cc.libcxx.cxxabi.libName}";
          };

          buildInputs = lib.optionals stdenv.isDarwin [ pkgs.libiconv pkgs.darwin.apple_sdk.frameworks.Security ];

          nativeBuildInputs = with rustPlatform; [ cargoSetupHook maturinBuildHook bindgenHook ];
        };

      in

      with pkgs;
      {
        packages = {

          lib = rustPlatform.buildRustPackage ( commonArgs // {
            name = "libsourmash";
            copyLibs = true;
            cargoLock.lockFile = ./Cargo.lock;
            nativeBuildInputs = with rustPlatform; [ bindgenHook ];
          });

          sourmash = python.buildPythonPackage ( commonArgs // rec {
            pname = "sourmash";
            version = "4.8.6";
            format = "pyproject";

            cargoDeps = rustPlatform.importCargoLock {
              lockFile = ./Cargo.lock;
            };

            propagatedBuildInputs = with python; [ cffi deprecation cachetools bitstring numpy scipy matplotlib screed ];

            DYLD_LIBRARY_PATH = "${self.packages.${system}.lib}/lib";
          });

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

        devShells.default = pkgs.mkShell.override { stdenv = stdenv; } (commonArgs // {
          nativeBuildInputs = with rustPlatform; [ bindgenHook ];

          buildInputs = [
            rustVersion
            openssl
            pkg-config

            git
            stdenv.cc.cc.lib
            (python312.withPackages (ps: with ps; [ virtualenv ]))
            (python311.withPackages (ps: with ps; [ virtualenv tox cffi ]))
            (python310.withPackages (ps: with ps; [ virtualenv ]))

            rust-cbindgen
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
            #cargo-semver-checks
            nixpkgs-fmt
            cargo-llvm-cov
          ];

          shellHook = ''
              export MACOSX_DEPLOYMENT_TARGET=10.14
            '';

          # Needed for matplotlib
          LD_LIBRARY_PATH = lib.makeLibraryPath [ pkgs.stdenv.cc.cc.lib ];

          # workaround for https://github.com/NixOS/nixpkgs/blob/48dfc9fa97d762bce28cc8372a2dd3805d14c633/doc/languages-frameworks/python.section.md#python-setuppy-bdist_wheel-cannot-create-whl
          SOURCE_DATE_EPOCH = 315532800; # 1980
        });
      });
}
