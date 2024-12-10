#!/bin/bash
# Script to make a distributable version of STPSF, with various packaging tweaks
set -e

if ! [[ $1 ]]; then
  echo "Provide a version string, e.g.:"
  echo "    ./make-data-sdist.sh 0.3.3"
  exit 1
fi

if ! [[ $DATAROOT ]]; then
  DATAROOT="/grp/jwst/ote/stpsf-data-source/"
fi
echo "Using data from $DATAROOT"

# If on Mac OS, tell tar to not include ._* files for
# HFS-specific extended attributes
export COPYFILE_DISABLE=1

# Prepare to create the data tarfile
# Also exclude various things we don't want to distribute, like .svn, the old OPDs, and the data source directories

TMPDIR="/tmp/stpsf-data"

mkdir -p "$TMPDIR"
rsync -avz --delete --exclude '._*' --exclude '_Obsolete' \
    --exclude .svn --exclude OPD_RevT --exclude TFI --exclude .DS_Store \
    --exclude sources --exclude "*py" --exclude "OTE_source" \
    --exclude "SI_WFE_source" --exclude README_DEVEL.md \
    "$DATAROOT" "$TMPDIR"

VER="$1"
echo "$VER" > $TMPDIR/version.txt
echo "Saving version number $VER to version.txt"


# Some temporary extras to support pre- and post-
# renaming of some WFIRST stuff
#ln -s $TMPDIR/WFI $TMPDIR/WFIRSTImager
#ln -s $TMPDIR/AFTA_WFC_C5_Pupil_Shortwave_Norm_2048px.fits $TMPDIR/AFTA_symmetrical.fits


# create distributable tar file
tar -cvz -C "$TMPDIR/.." \
    -f "stpsf-data-$VER.tar.gz" stpsf-data

echo "File output to:    $PWD/stpsf-data-$VER.tar.gz"
echo
echo "If that works, remove $TMPDIR"
