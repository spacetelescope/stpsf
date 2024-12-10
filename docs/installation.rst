.. _installation:

Requirements & Installation
===========================

.. note::

   This is entirely optional, but you may wish to sign up to the mailing list ``stpsf-users@maillist.stsci.edu``. This is a very low-traffic moderated announce-only list, to which we will periodically post announcements of updates to this software.

   To subscribe, visit  the `maillist.stsci.edu server <https://maillist.stsci.edu/scripts/wa.exe?SUBED1=STPSF-users&A=1>`_


.. NOTE: installation with conda is unavailable as of v1.1.0. uncomment and edit the following section once it's back.
  .. _install_with_conda:

  Recommended: Installing with conda
  ----------------------------------

  If you already use ``conda``, but do not want to install the full suite of STScI software, you can simply add the AstroConda *channel* and install STPSF as follows (creating a new environment named ``stpsf-env``)::

    $ conda config --add channels http://ssb.stsci.edu/astroconda
    $ conda create -n stpsf-env stpsf
    $ conda activate stpsf-env

  Upgrading to the latest version is done with ``conda update -n stpsf-env --all``.

  .. warning::

     You *must* install STPSF into a specific environment (e.g. ``stpsf-env``); our conda package will not work if installed into the default "root" environment.

.. _install_pip:

Installing with pip
-------------------

STPSF and its underlying optical library POPPY may be installed from the `Python Package Index <http://pypi.python.org/pypi>`_ in the usual manner for Python packages. ::

    $ pip install --upgrade stpsf
    [... progress report ...]

    Successfully installed stpsf

Note that ``pip install stpsf`` only installs the program code. **If you install via pip, you must manually download and install the data files, as** :ref:`described <data_install>` **below.**
To obtain source spectra for calculations, you should also follow :ref:`installation instructions for synphot <synphot_install>`.

.. note::
  Installation through conda is not available as of STPSF version 1.1.0. Conda
  users should instead follow the insructions in the preceding section to
  install via pip.


.. _synphot_install:

Installing or updating synphot
--------------------------------

Stsynphot is an optional dependency, but is highly recommended. Its installation instructions can be found in `the synphot docs <https://synphot.readthedocs.io/en/latest/#installation-and-setup>`_ or `a discussion in the POPPY docs <http://poppy-optics.readthedocs.io/en/stable/installation.html#installing-or-updating-synphot>`_.

.. _data_install:

Installing the Required Data Files
----------------------------------

*If you install via pip or manually*, you must install the data files yourself.

.. (If you install via Conda, the data files are automatically installed, in
    which case you can skip this section.) [uncomment once conda installation is
    available again]

Files containing such information as the JWST pupil shape, instrument throughputs, and aperture positions are distributed separately from STPSF. To run STPSF, you must download these files and tell STPSF where to find them using the ``STPSF_PATH`` environment variable.

1. Download the following file:  `stpsf-data-LATEST.tar.gz <https://stsci.box.com/shared/static/qxpiaxsjwo15ml6m4pkhtk36c9jgj70k.gz>`_  [approx. 70 MB]
2. Untar ``stpsf-data-LATEST.tar.gz`` into a directory of your choosing.
3. Set the environment variable ``STPSF_PATH`` to point to that directory. e.g. ::

    export STPSF_PATH=$HOME/data/stpsf-data

for bash. (You will probably want to add this to your ``.bashrc``.)

You should now be able to successfully ``import stpsf`` in a Python session.

.. warning::

   If you have previously installed the data files for an earlier version of STPSF, and then update to a newer version, the
   software may prompt you that you must download and install a new updated version of the data files.

.. Note::

   **For STScI Users Only:** Users at STScI may access the required data files from the Central Storage network. Set the following environment variables in your ``bash`` shell. (You will probably want to add this to your ``.bashrc``.) ::
      export STPSF_PATH="/grp/jwst/ote/stpsf-data"
      export PYSYN_CDBS="/grp/hst/cdbs"

Software Requirements
---------------------


See `the requirements.txt specification file <https://github.com/spacetelescope/stpsf/blob/develop/requirements.txt>`_ for the required package dependencies.

**Required Python version**: STPSF 1.1 and above require Python 3.10 or higher.

The major dependencies are the standard `NumPy, SciPy <http://www.scipy.org/scipylib/download.html>`_, `matplotlib <http://matplotlib.org>`_ stack, and `Astropy <http://astropy.org>`_.

**Recommended Python packages**:

* `synphot <https://synphot.readthedocs.io/>`_ enables the simulation
  of PSFs with proper spectral response to realistic source spectra.  Without
  this, PSF fidelity is reduced. See above for :ref:`installation instructions
  for synphot <synphot_install>`.  Stsynphot is recommended for most users.

**Optional Python packages**:

Some calculations with POPPY can benefit from the optional packages `psutil <https://pypi.python.org/pypi/psutil>`_ and `pyFFTW <https://pypi.python.org/pypi/pyFFTW>`_, but these are not needed in general. See `the POPPY installation docs <http://poppy-optics.readthedocs.io/en/stable/installation.html>`_ for more details.
These optional packages are only worth adding for speed improvements if you are spending substantial time running calculations.

Additional packages are needed for the optional use of GPUs to accelerate calculations. See the POPPY documentation.

.. _install_dev_version:

Installing a pre-release version or contributing to STPSF development
-----------------------------------------------------------------------

The `STPSF source code repository <https://github.com/spacetelescope/stpsf>`_ is hosted at GitHub, as is the repository for `POPPY <https://github.com/spacetelescope/poppy>`_. Users may clone or fork in the usual manner. Pull requests with code enhancements welcomed.

To install the current development version of STPSF, you can use ``pip`` to install directly from a ``git`` repository. To install STPSF and POPPY from ``git``, uninstall any existing copies of STPSF and POPPY, then invoke pip as follows::

    $ pip install -e git+https://github.com/spacetelescope/poppy.git#egg=poppy \
       -e git+https://github.com/spacetelescope/stpsf.git#egg=stpsf

This will create directories ``./src/poppy`` and ``./src/stpsf`` in your current directory containing the cloned repository. If you have commit access to the repository, you may want to clone via ssh with a URL like ``git+ssh://git@github.com:spacetelescope/stpsf.git``. Documentation of the available options for installing directly from Git can be found in the `pip documentation <http://pip.readthedocs.org/en/latest/reference/pip_install.html#git>`_.

Remember to :ref:`install the required data files <data_install>`, if you have not already installed them.
