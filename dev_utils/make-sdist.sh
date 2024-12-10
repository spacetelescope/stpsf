#!/bin/sh
# Script to make a distributable version of STPSF, with various packaging tweaks

VER="0.2.6"
TAR=/usr/bin/tar  # make sure to use the BSD version, required for the -L option

# Create a source distribution
#   use python sdist, but then re-make the tar file so we can
#   include the stsci_distutils_hack and defsetup files
python setup.py sdist
cd dist
tar xvzf stpsf-$VER.tar.gz
\cp ../stsci_distutils_hack.py stpsf-$VER
\cp ../defsetup.py stpsf-$VER
$TAR -cvz --exclude .svn --exclude '*.pyc' -f stpsf-$VER.tar.gz stpsf-$VER
\cp stpsf-$VER.tar.gz ~/web/software/stpsf
mv stpsf-$VER.tar.gz ..
cd ..


