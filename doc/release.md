# Releasing a new version of sourmash

These are adapted from the khmer release docs, originally written by
Michael Crusoe.

## Checklist

Here's a checklist to copy/paste into an issue:

```

Release candidate testing:
- [ ] Command line tests pass for a release candidate
- [ ] All eight release candidate wheels are built

Releasing to PyPI:

- [ ] RC tag(s)s deleted on github
- [ ] Release tag cut
- [ ] Release notes written
- [ ] All eight release wheels built
- [ ] Release wheels uploaded to pypi
- [ ] tar.gz distribution uploaded to pypi

After release to PyPI and conda-forge/bioconda packages built:

- [ ] [PyPI page](https://pypi.org/project/sourmash/) updated
- [ ] Zenodo DOI successfully minted upon new github release - [see search results](https://zenodo.org/search?page=1&size=20&q=sourmash)
- [ ] `pip install sourmash` installs the correct version
- [ ] `mamba create -n smash-release -y sourmash` installs the correct version
```

## Creating the build environment with conda

You can most easily set up your build environment with conda.

Your conda version will need to be at least `v4.9.0`. You can check your
conda version with `conda --version` and update with `conda update conda`.

Create the basic build environment:

```
mamba create -y -n sourmash-rc python=3.10 pip \
    cxx-compiler make twine tox tox-conda rust
```

Then activate it with `conda activate sourmash-rc`.

## Writing release notes

Draft release notes can be created with `git log --oneline
v4.6.1..latest`, but should then be edited manually. We suggest
putting PRs in the following categories:

```
Major new features:

Minor new features:

Bug fixes:

Cleanup and documentation updates:

Developer updates:

Dependabot updates:
```

A convenient way to edit release notes is to put them in a [hackmd.io](https://hackmd.io) document and edit/display them there; then, create a "draft release notes for v..." issue and paste the markdown into the issue.

## Testing a release

0\. First things first: check if Read the Docs is building properly for `latest`.
The build for the `latest` branch on [Read the Docs] should be passing,
and also the [rendered docs] should be up to date.

[Read the Docs]: https://readthedocs.org/projects/sourmash/builds/
[rendered docs]: https://sourmash.readthedocs.io/en/latest/

1\. The below should be done in a clean checkout:
```
cd $(mktemp -d)
git clone https://github.com/sourmash-bio/sourmash
cd sourmash
```

2\. Set your new version number and release candidate.
You might want to check [the releases page] for next version number,
or you can run `make last-tag` and check the output.
```
new_version=4.X.X
rc=rc1
```

Next create a new branch to work on release candidates and the version bump:
```
git checkout -b release/v${new_version}
```
and update the version number in `pyproject.toml` and `flake.nix`:
```
sed -i -e "s|version = .*$|version = \"${new_version}\"|g" pyproject.toml flake.nix
```

Commit the changes and push the branch:
```
git add pyproject.toml
git commit -m "${new_version} release"
git push -u origin release/v${new_version}
```
and then open a PR for the new branch by following the link printed by
```
echo "https://github.com/sourmash-bio/sourmash/pull/new/release/v${new_version}"
```

[the releases page]: https://github.com/sourmash-bio/sourmash/releases

Once the checks for the PR work, let's trigger the automatic wheel building
by creating a tag:

```
git tag -a v${new_version}${rc} -m "${new_version} release candidate ${rc}"
git push origin refs/tags/v${new_version}${rc}
```

3\. Test the release candidate. Bonus: repeat on macOS:
```
python -m pip install -U pip

cd ..
python -m venv testenv1
python -m venv testenv2
python -m venv testenv3
python -m venv testenv4

# First we test the tag

cd testenv1
source bin/activate
git clone --depth 1 --branch v${new_version}${rc} https://github.com/sourmash-bio/sourmash.git
cd sourmash
python -m pip install -r requirements.txt
pytest && cargo test

# Secondly we test via pip

cd ../../testenv2
deactivate
source bin/activate
python -m pip install build
python -m pip install -e git+https://github.com/sourmash-bio/sourmash.git@v${new_version}${rc}#egg=sourmash[test]
cd src/sourmash
pytest && cargo test
make dist
cp dist/sourmash*tar.gz ../../../testenv3/

# Is the distribution in testenv2 complete enough to build another
# functional distribution?

cd ../../../testenv3/
deactivate
source bin/activate
python -m pip install sourmash*tar.gz
tar xzf sourmash-${new_version}${rc}.tar.gz
cd sourmash-${new_version}${rc}
python -m pip install -r requirements.txt
cp -a ../../sourmash/tests/test-data tests/  ## We don't ship the test data, so let's copy it here
pytest && cargo test
```

4\. Do any final testing:

 * check that the binder demo notebook is up to date

5\. Wait for GitHub Actions to finish running on the release candidate tag.

Wait for the
[various cibuildwheel actions](https://github.com/sourmash-bio/sourmash/actions)
to finish and upload; the
[latest release](https://github.com/sourmash-bio/sourmash/releases)
should have eight wheel files attached to it.

6\. Remove release candidate tags

NOTE: If you delete the rc tag before the rc wheels are done building, they
may get added to the wrong release.

```
cd ../sourmash
git tag -d v${new_version}${rc}
git push --delete origin v${new_version}${rc}
```

## How to make a final release

When you've got a thoroughly tested release candidate,
cut a release like so:

1\. Merge the pull request bumping the version. Once the PR is merged,
change back to the `latest` branch and pull the new commit:

```
git checkout latest
git pull --rebase
```

2\. Create the final tag and push to GitHub:

```
git tag -a v${new_version} -m "${new_version} release"
git push --tags origin
```

(make sure to be in the `latest` branch when creating the final tag!)

3\. Upload wheels from GitHub Releases to PyPI

[GitHub Actions will automatically build wheels and upload them to GitHub Releases](https://github.com/sourmash-bio/sourmash/actions?query=workflow%3Acibuildwheel).
This will take about 45 minutes, or more. After they're built, they must be
copied over to PyPI manually.

You can do this in two ways: you can manually download all the files
from [the releases page], or, if you have the
[`GitHub CLI`](https://cli.github.com/), you can use that to download the
packages.

Download the wheels with the `GitHub CLI`:
```
mkdir -p wheel && cd wheel
gh release download v${new_version}
```
or download them manually.

Once you have them downloaded, upload them to PyPI like so:
```
twine upload *.whl
```
twine will correctly determine the version from the filenames.

4\. Once the wheels are uploaded, publish the new release on PyPI (requires an authorized account).
```
cd ..
make dist
twine upload dist/sourmash-${new_version}.tar.gz
```

(This must be done *after* the wheels are available, because some of
the conda package build steps require the source dist and are automatically
triggered when a new version shows up on PyPI.)

5\. Edit the release on GitHub; there will already be one associated
with the tag you pushed. Copy and paste in the release notes.

Note that there will also be releases associated with the Rust `core`
package, which is versioned differently than `sourmash`.  These will
be of the form `rXX.YY.ZZ`, e.g. `r0.9.0`. Please just ignore them :)

## Conda-forge

The [sourmash-minimal feedstock](https://github.com/conda-forge/sourmash-minimal-feedstock/)
in [conda-forge](https://conda-forge.org/) picks up new versions from
PyPI (need the sdist to be published) and opens a new PR.

Check if there are any dependency changes,
with special attention to the minimum supported Rust version.

After tests pass,
merge it and wait for the `sourmash-minimal` package to show up in conda-forge:
```
conda search sourmash-minimal={new_version}
```

[An example conda-forge PR for `4.6.0`](https://github.com/conda-forge/sourmash-minimal-feedstock/pull/37).

[An example bioconda PR for `4.6.0`](https://github.com/bioconda/bioconda-recipes/pull/38205).

## Bioconda

The BiocondaBot has an `autobump` feature that should pick up new releases from PyPI,
and open a PR in Bioconda. Review any changes
(especially dependency versions, since these don't get picked up).

Note that you need to wait for the `sourmash-minimal` package
prepared in the previous section to be available for installation,
and tests are going to fail in Bioconda before that.

An example PR for [`3.4.0`](https://github.com/bioconda/bioconda-recipes/pull/23171).

## Announce it!

If a bioinformatics software is released and no one tweets, is it really released?

Examples:

- [3.4.1](https://twitter.com/ctitusbrown/status/1286652952828993537)
- [3.4.0](https://twitter.com/luizirber/status/1283157954598858752)
- [3.3.0](https://twitter.com/ctitusbrown/status/1257418140729868291)
- [3.2.0](https://twitter.com/luizirber/status/1221923762523623425)
- [3.1.0](https://twitter.com/luizirber/status/1217639572202409984)
- [3.0.0](https://twitter.com/luizirber/status/1213588144458649600)
- [2.3.0](https://twitter.com/luizirber/status/1198027116396171264)
- [2.2.0](https://twitter.com/luizirber/status/1179126660911661057)
- [2.1.0](https://twitter.com/luizirber/status/1166910335120314369)
- [2.0.1](https://twitter.com/luizirber/status/1136786447518711808)
- [2.0.0](https://twitter.com/luizirber/status/1108846466502520832)
