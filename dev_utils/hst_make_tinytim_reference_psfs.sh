#!/usr/bin/env bash
# hst_make_tinytim_reference_psfs.sh: generate TinyTim WFC3 (IR and UVIS) reference PSFs non-interactively,
# for validating the STPSF HST models (see run_hst_regression_tests.py).
#
# Needs a compiled TinyTim installation (tiny1, tiny2, tiny3), pointed to by the TINYTIM environment variable.
#
# Two modes:
#
# 1) Standard suite (default when run with no camera/filter options, or with -S):
#      hst_make_tinytim_reference_psfs.sh [-S] [-o OUTDIR]
#    Generates the matrix below, at 3.0" diameter, G8V star (spectrum #13), zero despace:
#      WFC3 IR (detector 1014x1014):
#        filters F105W, F125W, F160W, plus monochromatic 1600 nm
#        positions (507,507) center, (50,50) lower-left corner, (960,960) upper-right corner
#      WFC3 UVIS (UVIS1 = upper chip, UVIS2 = lower chip; each 4096x2051):
#        filters F275W, F438W, F606W, F814W, plus monochromatic 600 nm
#        positions UVIS2 (2048,1024) chip center, UVIS2 (100,100) lower-left corner of the full detector,
#                  UVIS1 (2048,1024) chip center, UVIS1 (4000,1950) upper-right corner of the full detector
#    That is 12 IR + 20 UVIS = 32 PSFs, taking roughly 1-2 minutes. A MANIFEST.csv lists every case and the
#    field-dependent Zernikes TinyTim applied there.
#
# 2) Single PSF:
#      hst_make_tinytim_reference_psfs.sh -c ir|uvis1|uvis2 [-f FILTER | -m WAVELENGTH_NM] [-x X] [-y Y]
#                             [-s SPECTRUM_NUM] [-d DIAMETER_ARCSEC] [-z DESPACE_MICRONS]
#                             [-r ROOTNAME] [-o OUTDIR]
#    Defaults: F160W (IR) or F606W (UVIS), chip center, spectrum 13, 3.0", despace 0.
#
# Outputs per case, in OUTDIR (default /tmp/tinytim_ref):
#   ROOT.par         tiny1 parameter file (contains the full .pup/.tab; edit it to change aberrations)
#   ROOT00_psf.fits  tiny2 undistorted, finely sampled PSF; no detector kernel
#                    (pixel scale ~ lambda_min/(2D)/1.3; see PIXSCALE in the header)
#   ROOT00.fits      tiny3 detector-sampled PSF (IR 0.13"/px, UVIS 0.04"/px) with geometric distortion
#                    and the IPC (IR) or charge-diffusion (UVIS) kernel applied
#   ROOT.*.log       stdout from each stage
# Header keywords FOCUS, X_ASTIG, Y_ASTIG, X_COMA, Y_COMA hold the applied Z4-Z8 (RMS waves at 547 nm).
#
# Dumping maps: set the three "Write out ..." flags in ROOT.par to 1 and rerun tiny2 in OUTDIR with
# TINYTIM set. It writes tinytim_pupil.fits, tinytim_wave.fits and critpsf.fits at the first
# wavelength, then stops.
#
# Note: if tiny1 runs out of stdin it loops forever and floods its output, so it runs under an
# alarm and its output is capped.

set -euo pipefail

if [ -z "${TINYTIM:-}" ]; then
  echo "Set the TINYTIM environment variable to the TinyTim installation directory" >&2
  exit 1
fi
REPO="$(cd "$TINYTIM" && pwd)"
export TINYTIM="$REPO"

usage() { sed -n '2,42p' "$0" | sed 's/^# \{0,1\}//'; }

CAMERA=""; FILTER=""; MONO=""; X=""; Y=""; SPEC=13; DIAM=3.0; DESPACE=0
ROOT=""; OUTDIR="/tmp/tinytim_ref"; SUITE=0

while getopts "Sc:f:m:x:y:s:d:z:r:o:h" opt; do
  case $opt in
    S) SUITE=1 ;;
    c) CAMERA="$OPTARG" ;;
    f) FILTER="$OPTARG" ;;
    m) MONO="$OPTARG" ;;
    x) X="$OPTARG" ;;
    y) Y="$OPTARG" ;;
    s) SPEC="$OPTARG" ;;
    d) DIAM="$OPTARG" ;;
    z) DESPACE="$OPTARG" ;;
    r) ROOT="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    h|*) usage; exit 0 ;;
  esac
done
if [ -z "$CAMERA$FILTER$MONO$X$Y$ROOT" ]; then SUITE=1; fi

for exe in tiny1 tiny2 tiny3; do
  if [ ! -x "$REPO/$exe" ]; then
    echo "Missing $REPO/$exe; build with 'make macosx' (or 'make linux') in $REPO" >&2
    exit 1
  fi
done

mkdir -p "$OUTDIR"
cd "$OUTDIR"
MANIFEST="$OUTDIR/MANIFEST.csv"

# Print the value of a FITS header keyword (first HDU), or an empty string.
fits_key() {
  head -c 57600 "$1" | fold -w 80 | awk -v k="$2" '
    substr($0,1,8)=="END     " {exit}
    { key=substr($0,1,8); sub(/ +$/,"",key)
      if (key==k) { v=substr($0,11); sub(/\/.*/,"",v); gsub(/[ \x27]/,"",v); print v; exit } }'
}

# run_one CAMERA FILTER MONO X Y ROOT
#   CAMERA is ir, uvis1 or uvis2. Pass FILTER for a filter PSF, or MONO (nm) for a monochromatic PSF.
run_one() {
  local cam="$1" filt="$2" mono="$3" x="$4" y="$5" root="$6" pre answers
  case "$cam" in
    ir)    pre='23\n' ;;
    uvis1) pre='22\n1\n' ;;
    uvis2) pre='22\n2\n' ;;
    *) echo "Unknown camera '$cam' (use ir, uvis1, uvis2)" >&2; return 1 ;;
  esac
  # tiny1 prompts: camera [+ chip], position, filter or mono+wavelength,
  # spectrum (list choice 1 + number; skipped for mono), diameter, despace, rootname.
  if [ -n "$mono" ]; then
    answers=$(printf "${pre}%s %s\nmono\n%s\n%s\n%s\n%s\n" "$x" "$y" "$mono" "$DIAM" "$DESPACE" "$root")
  else
    answers=$(printf "${pre}%s %s\n%s\n1\n%s\n%s\n%s\n%s\n" "$x" "$y" "$filt" "$SPEC" "$DIAM" "$DESPACE" "$root")
  fi

  rm -f "$root.par" "${root}00.fits" "${root}00_psf.fits"
  printf '%s\n' "$answers" | perl -e 'alarm 30; exec @ARGV' "$REPO/tiny1" "$root.par" 2>&1 \
    | head -c 200000 > "$root.tiny1.log" || true
  if [ ! -s "$root.par" ]; then
    echo "tiny1 failed for $root; see $OUTDIR/$root.tiny1.log" >&2
    return 1
  fi
  "$REPO/tiny2" "$root.par" > "$root.tiny2.log" 2>&1
  "$REPO/tiny3" "$root.par" > "$root.tiny3.log" 2>&1
  if [ ! -s "${root}00.fits" ]; then
    echo "tiny2/tiny3 failed for $root; see logs in $OUTDIR" >&2
    return 1
  fi

  local t2="${root}00_psf.fits" t3="${root}00.fits" nw
  nw=$(awk '/# Number of wavelengths/ {print $1; exit}' "$root.par")
  printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
    "$root" "$cam" "${filt:-mono}" "$mono" "$x" "$y" "$([ -n "$mono" ] && echo "" || echo "$SPEC")" \
    "$DIAM" "$DESPACE" "$nw" \
    "$(fits_key "$t2" FOCUS)" "$(fits_key "$t2" X_ASTIG)" "$(fits_key "$t2" Y_ASTIG)" \
    "$(fits_key "$t2" X_COMA)" "$(fits_key "$t2" Y_COMA)" \
    "$t2" "$(fits_key "$t2" PIXSCALE)" "$t3" "$(fits_key "$t3" PIXSCALE)" >> "$MANIFEST"
  echo "  done: $root"
}

if [ ! -f "$MANIFEST" ]; then
  echo "root,camera,filter,mono_nm,x,y,spectrum_num,diameter_arcsec,despace_um,n_waves,Z4_focus,Z5_xastig,Z6_yastig,Z7_xcoma,Z8_ycoma,tiny2_file,tiny2_pixscale,tiny3_file,tiny3_pixscale" > "$MANIFEST"
fi

if [ "$SUITE" -eq 1 ]; then
  echo ">> Generating standard WFC3 reference suite in $OUTDIR"
  IR_POS=("507 507 ctr" "50 50 ll" "960 960 ur")
  for p in "${IR_POS[@]}"; do
    set -- $p
    for f in f105w f125w f160w; do run_one ir "$f" "" "$1" "$2" "ir_${f}_$3"; done
    run_one ir "" 1600 "$1" "$2" "ir_mono1600_$3"
  done
  UVIS_POS=("uvis2 2048 1024 ctr" "uvis2 100 100 ll" "uvis1 2048 1024 ctr" "uvis1 4000 1950 ur")
  for p in "${UVIS_POS[@]}"; do
    set -- $p
    for f in f275w f438w f606w f814w; do run_one "$1" "$f" "" "$2" "$3" "${1}_${f}_$4"; done
    run_one "$1" "" 600 "$2" "$3" "${1}_mono600_$4"
  done
else
  CAMERA="${CAMERA:-ir}"
  case "$CAMERA" in
    ir)  X="${X:-507}";  Y="${Y:-507}";  DEF_FILTER=f160w ;;
    uvis1|uvis2) X="${X:-2048}"; Y="${Y:-1024}"; DEF_FILTER=f606w ;;
    *) echo "Unknown camera '$CAMERA' (use ir, uvis1, uvis2)" >&2; exit 1 ;;
  esac
  if [ -z "$MONO" ]; then FILTER="${FILTER:-$DEF_FILTER}"; fi
  if [ -z "$ROOT" ]; then
    if [ -n "$MONO" ]; then ROOT="${CAMERA}_mono${MONO}_${X}_${Y}"; else ROOT="${CAMERA}_${FILTER}_${X}_${Y}"; fi
  fi
  run_one "$CAMERA" "$FILTER" "$MONO" "$X" "$Y" "$ROOT"
fi

echo ">> Manifest: $MANIFEST"
