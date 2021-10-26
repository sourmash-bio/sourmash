# Releasing a new version of sourmash

These are adapted from the khmer release docs, originally written by
Michael Crusoe.

## Required build environment

The basic build environment needed below can be created as follows:

```
conda create -y -n sourmash-rc python=3.8 pip cxx-compiler make twine tox tox-conda setuptools_scm
```

Then activate it with `conda activate sourmash-rc`.

You will also need a Rust installation (see
[Development Environment](developer.md)); be sure to update it to the
latest version with `rustup update`.

Your conda version will need to be at least `v4.9.0`. You can check your
conda version with `conda --version` and update with `conda update conda`.

## Testing a release

0\. First things first: check if Read the Docs is building properly for `latest`.
The build on [Read the Docs] should be passing,
and also check if the [rendered docs] are updated.

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
new_version=3.0.0
rc=rc1
```
and then tag the release candidate with the new version number prefixed by the letter 'v':
```
git tag -a v${new_version}${rc} -m "${new_version} release candidate ${rc}"
git push --tags origin
```

[the releases page]: https://github.com/sourmash-bio/sourmash/releases

3\. Test the release candidate. Bonus: repeat on macOS:
```
python -m pip install -U setuptools pip virtualenv wheel

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
python -m pip install -U setuptools pip wheel setuptools_scm
python -m pip install -r requirements.txt
make test

# Secondly we test via pip

cd ../../testenv2
deactivate
source bin/activate
python -m pip install -U setuptools pip wheel setuptools_scm
python -m pip install -e git+https://github.com/sourmash-bio/sourmash.git@v${new_version}${rc}#egg=sourmash[test]
cd src/sourmash
make test
make dist
cp dist/sourmash*tar.gz ../../../testenv3/

# Is the distribution in testenv2 complete enough to build another
# functional distribution?

cd ../../../testenv3/
deactivate
source bin/activate
python -m pip install -U setuptools pip wheel setuptools_scm
python -m pip install sourmash*tar.gz
python -m pip install pytest
tar xzf sourmash-${new_version}${rc}.tar.gz
cd sourmash-${new_version}${rc}
python -m pip install -r requirements.txt
make dist
cp -a ../../sourmash/tests/test-data tests/  ## We don't ship the test data, so let's copy it here
make test
```

4\. Publish the new release on the testing PyPI server.
You will need to [change your PyPI credentials].
We will be using `twine` to upload the package to TestPyPI and verify
everything works before sending it to PyPI:
```
python -m pip install twine
twine upload --repository-url https://test.pypi.org/legacy/ dist/sourmash-${new_version}${rc}.tar.gz
```
Test the PyPI release in a new virtualenv:
```
cd ../../testenv4
deactivate
source bin/activate
python -m pip install -U setuptools pip wheel
# install as much as possible from non-test server!
python -m pip install screed pytest numpy matplotlib scipy bam2fasta deprecation cffi
python -m pip install -i https://test.pypi.org/simple --pre sourmash
sourmash info  # should print "sourmash version ${new_version}${rc}"
```

[change your PyPI credentials]: https://packaging.python.org/tutorials/packaging-projects/#uploading-the-distribution-archives

5\. Do any final testing:

 * check that the binder demo notebook is up to date

## How to make a final release

When you've got a thoroughly tested release candidate,
cut a release like so:

1\. Create the final tag. Write the changes from previous version in the tag commit message. `git log --oneline` can be useful here, because it can be used to compare the two versions (and hopefully we used descriptive PR names and commit messages). An example comparing `2.2.0` to `2.1.0`:
`git log --oneline v2.1.0..v2.2.0`

```
cd ../sourmash
git tag -a v${new_version}
```

2\. Delete the release candidate tag and push the tag updates to GitHub:
```
git tag -d v${new_version}${rc}
git push --tags origin
git push --delete origin v${new_version}${rc}
```

3\. Upload wheels from GitHub Releases to PyPI

[GitHub Actions will automatically build wheels and upload them to GitHub Releases](https://github.com/sourmash-bio/sourmash/actions?query=workflow%3Acibuildwheel).
This will take about 45 minutes, or more. After they're built, they must be
copied over to PyPI manually.

You can do this in two ways: you can manually download all the files
from [the releases page], or, if you have
[`hub`](https://hub.github.com/), you can use that to download the
packages.

Download the wheels with hub:
```
mkdir -p wheel && cd wheel
hub release download v${new_version}
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

(This should be done *after* the wheels are available, because some of
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

An example PR for [`3.4.0`](https://github.com/conda-forge/sourmash-minimal-feedstock/pull/7).

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
