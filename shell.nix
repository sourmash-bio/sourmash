# shell.nix
let
  sources = import ./nix/sources.nix;
  pkgs = import sources.nixpkgs { overlays = [ (import sources.rust-overlay) ]; };
  mach-nix = import sources.mach-nix {
    python = "python38";
  };

  customPython = mach-nix.mkPython {
    requirements = ''
       screed>=0.9
       cffi>=1.14.0
       numpy
       matplotlib
       scipy
       deprecation>=2.0.6
       cachetools >=4,<5
       setuptools>=38.6.0
       milksnake
       setuptools_scm>=3.2.0
       setuptools_scm_git_archive
       pytest
       pytest-cov
       hypothesis
       tox
    '';
  };

in
  with pkgs;

  pkgs.mkShell {
    buildInputs = [
      rust-bin.stable.latest.rust
      #customPython
      git
      stdenv.cc.cc.lib
      (python38.withPackages(ps: with ps; [ virtualenv tox ]))
      python39
      python37
      python36
    ];

    shellHook = ''
       # workaround for https://github.com/NixOS/nixpkgs/blob/48dfc9fa97d762bce28cc8372a2dd3805d14c633/doc/languages-frameworks/python.section.md#python-setuppy-bdist_wheel-cannot-create-whl
       export SOURCE_DATE_EPOCH=315532800 # 1980
       export LD_LIBRARY_PATH="${stdenv.cc.cc.lib}/lib64:$LD_LIBRARY_PATH";
    '';
  }
