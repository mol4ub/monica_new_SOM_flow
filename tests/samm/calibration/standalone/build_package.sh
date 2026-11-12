#!/bin/bash
# Builds the standalone SAMM package (Windows via MinGW-w64 from nixpkgs, and Linux) into <out_dir>.
# usage (repository root): tests/samm/calibration/standalone/build_package.sh <out_dir>
set -e
OUT=$1; mkdir -p "$OUT/params" "$OUT/output" "$OUT/source"
H=tests/samm/calibration/standalone
MCF=$(nix --extra-experimental-features 'nix-command flakes' build --inputs-from . 'nixpkgs#pkgsCross.mingwW64.windows.mcfgthreads^*' --no-link --print-out-paths)
INC=$(echo "$MCF" | grep -- '-dev$'); LIB=$(echo "$MCF" | grep -v -- '-dev$')
nix --extra-experimental-features 'nix-command flakes' shell --inputs-from . nixpkgs#pkgsCross.mingwW64.buildPackages.gcc -c \
  x86_64-w64-mingw32-g++ -O3 -std=c++17 -D_USE_MATH_DEFINES -static -static-libgcc -static-libstdc++ \
  -I "$INC/include" -L "$LIB/lib" -I src/core tests/samm/calibration/samm_calib.cpp src/core/samm.cpp -o "$OUT/samm.exe"
g++ -O3 -std=c++17 -pthread -I src/core tests/samm/calibration/samm_calib.cpp src/core/samm.cpp -o "$OUT/samm_linux"
cp tests/samm/data/climate_file.csv tests/samm/data/AllMeasurements.csv "$OUT/"
cp $H/README.txt $H/run_samm.bat $H/calibrate.bat $H/plot_results.R "$OUT/"
cp $H/params/*.csv "$OUT/params/"
cp src/core/samm.h src/core/samm.cpp tests/samm/calibration/samm_calib.cpp "$OUT/source/"
