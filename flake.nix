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

    mach-nix = {
      url = "github:DavHau/mach-nix/3.5.0";
      inputs.nixpkgs.follows = "nixpkgs";
      inputs.flake-utils.follows = "utils";
      inputs.pypi-deps-db.follows = "pypi-deps-db";
    };

    pypi-deps-db = {
      url = "github:DavHau/mach-nix/3.5.0";
    };
  };

  outputs = { self, nixpkgs, naersk, rust-overlay, mach-nix, pypi-deps-db, utils }:
    utils.lib.eachDefaultSystem (system:
      let
        overlays = [ (import rust-overlay) ];
        pkgs = import nixpkgs {
          inherit system overlays;
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

        python = "python39";
        mach-nix-wrapper = import mach-nix { inherit pkgs python; };
      in

      with pkgs;
      {
        packages = {
          lib = naersk-lib.buildPackage {
            pname = "libsourmash";
            root = ./.;
            copyLibs = true;
          };
          sourmash = mach-nix-wrapper.buildPythonPackage {
            src = ./.;
            pname = "sourmash";
            version = "4.4.0";
            requirements = ''
              screed>=1.0.5
              cffi>=1.14.0
              numpy
              matplotlib
              scipy
              deprecation>=2.0.6
              cachetools<6,>=4
              bitstring<4,>=3.1.9
            '';
            requirementsExtra = ''
              setuptools
              milksnake
              setuptools_scm[toml] >= 4, <6
            '';
            SETUPTOOLS_SCM_PRETEND_VERSION = "4.4.0";
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
            (python39.withPackages (ps: with ps; [ virtualenv setuptools ]))
            (python38.withPackages (ps: with ps; [ virtualenv setuptools ]))

            rust-cbindgen

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
