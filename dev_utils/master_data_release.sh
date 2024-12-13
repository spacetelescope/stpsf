#!/bin/bash
# Top-level script to make a distributable version of the data files

# See /grp/stpsf/stpsf-data-source/README_DEVEL.txt  # TODO SAPP  VERIFY THIS IS CREATED

if ! [[ $1 ]]; then
  echo "Provide a version string, e.g.:"
  echo "    ./master_data_release.sh  0.3.3"
  exit 1
fi

VER="$1"
TMPDIR="/tmp/stpsf-data"

./make-data-sdist.sh $VER

echo
echo "Copying latest data to /grp/jwst/ote for internal stsci use..."
main_directory="/grp/jwst/ote"  # TODO SAPP  - update when stpsf new data mount
new_directory="$main_directory/stpsf-data-$VER"
symlink_directory="/grp/jwst/ote/stpsf-data"
legacy_webbpsf_symlink_directory="/grp/jwst/ote/webbpsf-data"

cp "$PWD/stpsf-data-$VER.tar.gz" "$main_directory"
mkdir "$new_directory"
tar -xzf "$PWD/stpsf-data-$VER.tar.gz" -C "$new_directory"
rm "$symlink_directory"
ln -s "$new_directory/stpsf-data" "$symlink_directory"

# Allow legacy webbpsf users to continue using the stpsf data
rm "$legacy_webbpsf_symlink_directory"
ln -s "$symlink_directory" "$legacy_webbpsf_symlink_directory"

./make-minimal-datafiles.py  ${PWD}/stpsf-data-${VER}.tar.gz $VER

echo
echo "================================================="
echo "Data extracted for internal use with updated symlink  $symlink_directory -> $new_directory"
echo "Legacy WebbPSF symlink $legacy_webbpsf_symlink_directory -> $symlink_directory"
echo
echo "OUTPUT FILES:"
echo
echo ${PWD}/stpsf-data-${VER}.tar.gz
echo ~/tmp/minimal-stpsf-data-${VER}/minimal-stpsf-data-${VER}.tar.gz
echo
echo You probably want to test if those look as expected, and if so then copy into the Box folder 'stpsf_data_public'
echo "================================================="
echo

