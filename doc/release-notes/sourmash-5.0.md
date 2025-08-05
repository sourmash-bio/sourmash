# sourmash v5.0 _draft_ release notes

```{contents}
   :depth: 2
```

Draft text: We are pleased to announce release 5.0 of sourmash! This
release contains many feature improvements and new functionality, as
well as many breaking changes with sourmash 4.0

Please see
[our migration guide](../support.md#migrating-from-sourmash-v4x-to-sourmash-v5x)
for guidance on updating to sourmash v5, and post questions about
migrating to sourmash 5.0 in the
[sourmash issue tracker](https://github.com/sourmash-bio/sourmash/issues/new)!

## Breaking changes for 5.0

### Numerical output and search results are mostly unchanged

* several `tax metagenome` output formats will produce different results: in particular, the `lineage_summary` and `krona` formats now use abundances, and `tax metagenome` requires abundance-weighted gather results. You can use `--ignore-abund` to recover the v4 behavior. See [#3711](https://github.com/sourmash-bio/sourmash/pull/3711) for details.

### New or changed behavior

* the default output format for `tax metagenome` is now `human`, instead of `csv_summary`. You can use `-F human` to recover the v4 behavior. ([#3711](https://github.com/sourmash-bio/sourmash/pull/3711))
* `sig manifest` now defaults to `--no-rebuild-manifest` - see [#3074](https://github.com/sourmash-bio/sourmash/pull/3074).
* `sig check` and `sig collect` now default to `--relpath` - see [#3071](https://github.com/sourmash-bio/sourmash/issues/3071).

### Feature removal

### Feature/function deprecations 

## Refactoring, improvements, and minor bug fixes

## Documentation updates

## Infrastructure and CI changes
