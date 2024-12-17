***************************************************
Developer Notes: Releasing a new version of STPSF
***************************************************

Prerequisites
=============
* If you are making a release for `poppy` at the same time as a release in STPSF, do that first. Update the dependency requirement to the new version of poppy, in ``stpsf/pyproject.toml``.
* Is the `develop` build `passing on Github Actions? <https://github.com/spacetelescope/stpsf/actions>`_ with all desired release items included?
* Does the latest release install properly?
    * ``$ conda env remove --name ENV_NAME``
	* ``$ conda create -n ENV_NAME python=3.xx``
	* ``$ conda activate ENV_NAME``
	* ``$ pip install -e .``
	* ``$ pip uninstall tox``
	* ``$ pip install tox tox-conda``

Releasing new data packages
===========================

#. Run ``dev_utils/master_data_release.sh`` to make a gzipped tarred archive of the STPSF data:
    #. If you are on the Institute network
        #. ``$ cd stpsf/dev_utils/``
        #. ``$ ./master_data_release.sh 0.X.Y``
    #. If you're working from a local data root
        #. ``$ cd stpsf/dev_utils/``
        #. ``$ DATAROOT="/Users/you/stpsf-data-sources/" ./make-data-sdist.sh 0.X.Y``
        #. ``$ cp ./stpsf-data-0.X.Y.tar.gz /where/ever/you/want/``
#. If the new data package is **required** (meaning you can't run STPSF without it, or you can run but may get incorrect results), you must bump ``DATA_VERSION_MIN`` in ``__init__.py`` to ``(0, X, Y)``
#. Extract the resulting data archive and check that you can run the STPSF tests with ``STPSF_PATH`` pointing to it
#. Copy the data archive into public web space. This now means on Box. The following steps need to be performed in this sequence in order to preserve the naming conventions.
    #. Find ``stpsf-data-LATEST.tar.gz``, and click on "more options" and "Update Version".  Choose the newest version of ``stpsf-data-#.#.#.tar.gz``
    #. This will change the name of ``stpsf-data-LATEST.tar.gz`` to be what you just uploaded, rename the file back to ``stpsf-data-LATEST.tar.gz``
    #. Upload to Box a separate version of ``stpsf-data-#.#.#.tar.gz`` shared data folder for future storage.
    #. Find ``minimal-stpsf-data-LATEST.tar.gz``, and click on "more options" and "Update Version".  Choose the newest version of ``minimal-stpsf-data-#.#.#.tar.gz``
    #. This will change the name of ``minimal-stpsf-data-LATEST.tar.gz`` to be what you just uploaded, rename the file back to ``minimal-stpsf-data-LATEST.tar.gz``
    #. Upload to Box a separate version of ``minimal-stpsf-data-#.#.#.tar.gz`` shared data folder for future storage.
    #. Verify the shared link of ``stpsf-data-latest.tar.gz`` is the same that exists in ``docs/installation.rst`` ("copy shared link" then "link settings")
#. A shared copy will be automatically configured in STScI Central Store with updated symlink ``/grp/stpsf/stpsf-data``
#. Verify code base is still up to date with box links and version names (they should be)
    #. Verify ``installation.rst`` with link to box data (this shouldn't need to change the box link, but verify it hasn't changed)
    #. update minimal in the ci setup (``stpsf/.github/workflows/download_data.yml``) (this also shouldn't need to change as the box link is for latest)
    #. update ``stpsf/stpsf/__init__.py`` with version number  (DATA_VERSION_MIN)
    #. CITATIONS.cff with new version
#. Generate the release notes
    #. You can do a draft release on github to autogenerate the relnotes
    #. In github set new release to be pre-release and make a release candidate tag -  1.3.0.rc1
    #. Auto-generate the release notes
    #. Copy output to ``docs/relnotes.rst``
#. Verify everything on sphinx locally on your computer (in docs directory $make html)
    #. make sure sphinx is installed properly with its needed dependencies
    #. navigate to docs folder
    #. ``$ pip install -r requirements.txt``
    #. ``$ pip uninstall graphviz``
    #. ``$ conda install graphviz``
    #. ``$ make clean`` (not needed for initial run, just to reset everything)
    #. ``$ make html``
#. Merge in any changes you've made including release notes
    #. NOTE: If the builds fail, this may be because we dont have the data cache fetched yet if so:
    #. Go to https://github.com/spacetelescope/stpsf/actions/workflows/download_data.yml
    #. Open "Run Workflow" box, change it to point to your pre-release-xxx branch (or whatever your pre release branch is named), run it
    #. It should then make a cache for you of the new version data
    #.	If you look at https://github.com/spacetelescope/stpsf/actions/caches, you will see what is available that the CI can pick from
    #.	THIS may still not solve the problem, as the cache is only in your branch, not the actual PR.  so if the branch passes, the PR should
    #. Theoretically be fine (despite its failure).  You can merge, and then re-run the "dowload data" action for develop, and then re-run your failed jobs in develop.
    #. Do not release unless develop is passing all tests
#. Update the test_readthedocs branch.  Force development there.  Test it on readthedocs (it should be hidden on the actual site).
    #. Checkout the branch you want to overwrite (test_readthedocs)
        #. ``$git checkout test_readthedocs``
    #. Reset the target branch to match the source branch (develop)
        #. ``$git reset --hard develop``
    #. Push to the github repo (probably upstream, may be origin, just dont do your personal one)
        #. ``$git push upstream test_readthedocs --force``
#. Once readthe docs looks all good test your release on test pypi.
    #. Create new env and install STPSF
    #. ``$ pip install build twine``
    #. ``$ python -m build``
    #. ``$ twine check dist/*``
    #. ``$ twine upload --repository-url https://test.pypi.org/legacy/ dist/* --verbose``  (NOTE: API token is the password in your ~/.pypirc testpypi token)
    #. test that you can download and install in fresh env (have pypi as backup for libraries that aren't on testpypi):
        #. ``$ pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ stpsf==<VERSION>``
#. Tag a version in develop and push it to git (do it through local terminal, not through website)
    #. ``$ git tag -a <release-tag> -m “webbpsf v1.4.0”`` (<release-tag> is just the version number --> 1.4.0)
    #. ``$ git push upstream <release-tag>``
#. Go to stable branch, and look at where it says how many commits behind it is from develop. Click that to generate a pull request (do not squash when you merge here)
#. When tests pass merge them to stable
#. Release on Github:
    #. On Github, click on ``[N] Releases``
    #. Select ``Draft a new release``.
    #. Specify the version number, title, and brief description of the release.
    #. Press ``Publish Release``
    #. Release to PyPI should now happen automatically on GitHub Actions. This will be triggered by a GitHub Actions build of a tagged commit on the `stable` branch.
#. Verify that files stored in ``/grp/stpsf/stpsf-data`` (symlink directory) have the correct permissions.
    #. ``$ cd /grp/stpsf/``
    #. ``$ find . -type f -exec chmod 755 {} \;`` (current and all subdirectories should be rwxr-xr-x)
#. Email an announcement to ``stpsf-users@maillist.stsci.edu``
