
# Releasing a new version of sourmash


These are adapted from the khmer release docs, originally written by
Michael Crusoe.

Remember to update release numbers/RC in:

* this document
* setup.py

## Testing a release


 1\. The below should be done in a clean checkout:
```
cd $(mktemp -d)
git clone git@github.com:dib-lab/sourmash.git
cd sourmash
```
2\. Set your new version number and release candidate:
```
        new_version=1.0
        rc=rc1
```
 and then tag the release candidate with the new version number prefixed by
   the letter 'v':
```
        git tag -a v${new_version}-${rc}
        git push --tags git@github.com:dib-lab/sourmash.git
```
3\. Test the release candidate. Bonus: repeat on Mac OS X:
```
        cd ..
        virtualenv testenv1
        virtualenv testenv2
        virtualenv testenv3
        virtualenv testenv4
        # First we test the tag

        cd testenv1
        source bin/activate
        git clone --depth 1 --branch v${new_version}-${rc} https://github.com/dib-lab/sourmash.git
        cd sourmash
        pip install -r requirements.txt
        make test
        pip uninstall -y sourmash; pip uninstall -y sourmash; make install
        mkdir ../not-sourmash # if there is a subdir named 'sourmash' py.test will execute tests
        # there instead of the installed sourmash module's tests
        pushd ../not-sourmash; py.test --pyargs sourmash; popd


        # Secondly we test via pip

        cd ../../testenv2
        source bin/activate
        pip install -U setuptools==3.4.1
        pip install -e git+https://github.com/dib-lab/sourmash.git@v${new_version}-${rc}#egg=sourmash[test]
        cd src/sourmash
        make test
        make dist
        cp dist/sourmash*tar.gz ../../../testenv3/
        pip uninstall -y sourmash; pip uninstall -y sourmash; make install
        cd ../.. # no subdir named sourmash here, safe for testing installed sourmash module
        py.test --pyargs sourmash

        # Is the distribution in testenv2 complete enough to build another
        # functional distribution?

        cd ../testenv3/
        source bin/activate
        pip install -U setuptools==3.4.1
        pip install sourmash*tar.gz
        pip install pytest
        tar xzf sourmash*tar.gz
        cd sourmash*
        pip install -r requirements.txt
        make dist
        make test
        pip uninstall -y sourmash; pip uninstall -y sourmash; make install
        mkdir ../not-sourmash
        pushd ../not-sourmash ; py.test  --pyargs sourmash ; popd
```
4\. Publish the new release on the testing PyPI server.  You will need
   to change your PyPI credentials as documented here:
   https://wiki.python.org/moin/TestPyPI.  You may need to re-register:
```
        python setup.py register --repository test
```
  Now, upload the new release:
```
        python setup.py sdist upload -r test
```
   Test the PyPI release in a new virtualenv:
```
        cd ../../testenv4
        source bin/activate
        pip install -U setuptools==3.4.1
        # install as much as possible from non-test server!
        pip install screed pytest numpy matplotlib scipy
        pip install -i https://testpypi.python.org/pypi --pre --no-clean sourmash
        py.test --pyargs sourmash
```
5\. Do any final testing:

   * check that the binder demo notebook is up to date

## How to make a final release


When you've got a thoroughly tested release candidate, cut a release like
so:

1\.Create the final tag and publish the new release on PyPI (requires an
   authorized account):
```
        cd ../sourmash
        git tag -a v${new_version}
        python setup.py register sdist upload
```
2\. Delete the release candidate tag and push the tag updates to GitHub:
```
        git tag -d v${new_version}-${rc}
        git push git@github.com:dib-lab/sourmash.git
        git push --tags git@github.com:dib-lab/sourmash.git
```
3\. Add the release on GitHub, using the tag you just pushed.  Name
   it 'version X.Y.Z', and copy and paste in the release notes:

4\. Make a binary wheel on OS X:
```
        virtualenv build
        cd build
        source bin/activate
        pip install -U setuptools==3.4.1 wheel
        pip install --no-clean sourmash==${new_version}
        cd ../
        python ./setup.py bdist_wheel upload
```
## To test on a blank Ubuntu system


```

   apt-cache update && apt-get -y install python-dev libfreetype6-dev && \
   pip install sourmash[test] && py.test --pyargs sourmash
```
