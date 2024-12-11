STPSF: Simulated Point Spread Functions for the James Webb and Nancy Grace Roman Space Telescopes
===================================================================================================

.. image:: https://github.com/spacetelescope/stpsf/blob/stable/docs/readme_fig.png?raw=true

.. image:: https://img.shields.io/pypi/v/stpsf.svg
   :target: https://pypi.python.org/pypi/stpsf
   :alt: Badge showing current released PyPI version

.. image:: https://github.com/spacetelescope/stpsf/workflows/CI/badge.svg?branch=develop
   :target: https://github.com/spacetelescope/stpsf/actions
   :alt: Github Actions CI Status

.. image:: https://codecov.io/gh/spacetelescope/stpsf/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/spacetelescope/stpsf

.. |Documentation Status| image:: https://img.shields.io/readthedocs/stpsf/latest.svg?logo=read%20the%20docs&logoColor=white&label=Docs&version=latest
   :target: https://stpsf.readthedocs.io/en/latest/
   :alt: Documentation Status

.. image:: https://img.shields.io/badge/ascl-1504.007-blue.svg?colorB=262255
   :target: http://ascl.net/1504.007


**ADVISORY:  STPSF IS CURRENTLY BEING MIGRATED FROM WEBBPSF - THIS REPOSITORY IS NOT READY FOR USE**

STPSF produces simulated PSFs for the James Webb Space Telescope, NASA's
flagship infrared space telescope. STPSF can simulate images for any of the
four science instruments plus the fine guidance sensor, including both direct
imaging, coronagraphic, and spectroscopic modes.

STPSF also supports simulating PSFs for the upcoming Nancy Grace Roman Space Telescope (formerly WFIRST),
including its Wide Field Instrument and a preliminary version of the Coronagraph Instrument.

.. note::

   The current Roman WFI optical model was provided by Goddard Space Flight Center circa 2021 (the Cycle 9 reference data); a new optical model is currently being implemented in STPSF.

Developed by Marshall Perrin, Joseph Long, Shannon Osborne, Robel Geda, Bradley Sappington, Marcio Meléndez,
Charles-Philippe Lajoie, Jarron Leisenring, Neil Zimmerman, Keira Brooks,
Justin Otor, Trey Kulp, Lauren Chambers, Alden Jurling, and collaborators, 2010-2024.

Documentation can be found online at https://stpsf.readthedocs.io

STPSF requires input data for its simulations, including optical path
difference (OPD) maps, filter transmission curves, and coronagraph Lyot mask
shapes. These data files are not included in this source distribution.
Please see the documentation to download the required data files.
