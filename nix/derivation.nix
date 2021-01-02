# TODO: check https://github.com/lightspeed/palisade for suggestions on setting
# up a nix dev env
# also https://rycwo.xyz/2019/02/16/nixos-series-dev-env

let
  moz_overlay = import (builtins.fetchTarball https://github.com/mozilla/nixpkgs-mozilla/archive/master.tar.gz);
  nixpkgs = import <nixpkgs> { overlays = [ moz_overlay ]; };
  ruststable = (nixpkgs.latest.rustChannels.stable.rust);

  mach-nix = import (
    builtins.fetchGit {
      url = "https://github.com/DavHau/mach-nix/";
      ref = "2.0.0";
    }
  );

  customPython = mach-nix.mkPython {
    python = nixpkgs.python38;
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
  with nixpkgs;

  nixpkgs.mkShell {
    buildInputs = [
      #customPython
      (python38.withPackages(ps: with ps; [ virtualenv tox ]))
      python39
      python37
      python36
      git
#      stdenv
      ruststable
      stdenv.cc.cc.lib
#      openssl
#      pkg-config
#      nasm
#      rustup
#      cmake
#      zlib
   ];

  shellHook = ''
     # workaround for https://github.com/NixOS/nixpkgs/blob/48dfc9fa97d762bce28cc8372a2dd3805d14c633/doc/languages-frameworks/python.section.md#python-setuppy-bdist_wheel-cannot-create-whl
     export SOURCE_DATE_EPOCH=315532800 # 1980
     export LD_LIBRARY_PATH="${stdenv.cc.cc.lib}/lib64:$LD_LIBRARY_PATH";
  '';
}
