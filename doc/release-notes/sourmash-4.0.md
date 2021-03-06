# sourmash v4.0 release notes

```{contents}
   :depth: 2
```

We are pleased to announce release 4.0 of sourmash! This release
contains many feature improvements and new functionality, as well as
many breaking changes with sourmash 2.x and 3.x.

Please see
[our migration guide](../support.md#migrating-from-sourmash-v3-x-to-sourmash-v4-x)
for guidance on updating to sourmash v4, and post questions about
migrating to sourmash 4.0 in the
[sourmash issue tracker](https://github.com/dib-lab/sourmash/issues/new).

## Major changes for 4.0

### Numerical output and search results are unchanged

There are no changes to numerical output or search results in this
release; you should get the same results with v4 as you get with v3,
except where command-line parameters need to be adjusted as noted
below (see: protein ksize #1277, lca summarize changes #1175, sourmash
gather on signatures without abundance #1328). Please
[file an issue](https://github.com/dib-lab/sourmash/issues) if your
results change!

### New or changed behavior

* default SBT storage is now .sbt.zip (#1174, #1170)
* add `sourmash sketch` command for creating signatures (#1159)
* protein ksizes in MinHash are now divided by 3, except in `sourmash compute` (#1277)
* refactor MinHash API and implementation: add, iadd, merge, hashes, and max_hash (#1282, #1154, #1139, #1301)
* add HyperLogLog implementation (#1223)
* `SourmashSignature.name` is now a property (not a method): use `str(sig)` instead of `name()` (#1179, #1232)
*  `lca summarize` no longer merges all signatures, and uses hash abundance by default (#1175)
* `index `and `lca index` (#1186, #1222) now support `--from-file` and no longer require signature files on command line
* `--traverse-directory` is now on by default for signature loading behavior (#1178)
* `sourmash sketch` and `sourmash compute` no longer create empty signatures from empty files and stdin (#1347)
* `sourmash sketch` and `sourmash compute` set `sig.filename` to empty string when filename is `-` (#1347)

### Feature removal

* remove Python 2.7 support (& end Python 2 compatibility) (#1145, #1144)
* remove `lca gather` (#1307)
* remove 10x support from `sourmash compute` (#1229)
* remove `dump` command (#1157)

### Feature/function deprecations 

* deprecate `sourmash compute` (#1159)
* deprecate `load_signatures`, `sourmash.load_one_signature`, `create_sbt_index`, and `load_sbt_index` (#1279, #1304)
* deprecate `import_csv` in favor of new `sourmash sig import --csv` (#1281)

## Refactoring, improvements, and minor bug fixes:

* accept file list in `sourmash sig cat` (#1236)
* add unique_intersect_bp and gather_result_rank to gather CSV output (#1219)
* remove deprecated minhash functions (#1149)
* fix Rust panic error in signature creation (#1172)
* cache nodes in SBT during search (#1161)
* fix two bugs in gather `--output-unassigned` (#1156)

## Documentation updates

* major update and cleanup of docs given new functionality; add sourmash sketch documentation (#1283)
* add information about versioning, migrations, etc to the docs (#1153, #1283)

## Infrastructure and CI changes:

* update finch requirement from 0.3.0 to 0.4.1 (#1290)
* update rand for test, and activate "js" feature for getrandom (#1275)
* dev updates (configs and doc) (#1298)
* move wheel building from Travis to GitHub Actions (#1295)
* fix new clippy warnings from Rust 1.49 (#1267)
* use tox for running tests locally (#696)
* CI: small build fixes (#1252)
* CI: Fix releases in GitHub Actions (#1250)
* update build_wheel action paths
* CI: moving python tests from travis to GH actions (#1249)
* CI: move wheel building to GitHub actions (#1244)
* remove last .rst file from docs (#1185)
* update CI for latest branch name change (#1150)
