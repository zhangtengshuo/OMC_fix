#!/usr/bin/env bash
set -Eeuo pipefail

BUILD_ROOT=/scratch/alpha/build-exactden-20260820_093910
SOURCE_ROOT=$BUILD_ROOT/src/OpenMolcasUP_ExactDenEmb-65a3eeb11
BUILD_DIR=$BUILD_ROOT/openmolcas-build-j4
INSTALL_PREFIX=/scratch/alpha/OpenMolcasUP_ExactDenEmb-65a3eeb11
WIGNERNJ_PREFIX=$BUILD_ROOT/wignernj-sys
OMC_PREFIX=/opt/mamba/envs/omc
GA_PREFIX=/opt/ga-5.8.2-mambaomc
EXPERIMENT_ROOT=/scratch/alpha/exactemb-ethylene-20260820_103419
LOG_DIR=$EXPERIMENT_ROOT/logs

mkdir -p "$LOG_DIR"
export PATH="$OMC_PREFIX/bin:/usr/bin:/bin:$PATH"
export LD_LIBRARY_PATH="$OMC_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export GAROOT=$GA_PREFIX
export MKLROOT=$OMC_PREFIX
unset CC CXX FC CFLAGS CXXFLAGS FCFLAGS LDFLAGS CPPFLAGS LIBS

[[ -f "$SOURCE_ROOT/src/Driver/exactemb.prgm.src" ]] || { echo "missing exactemb.prgm.src" >&2; exit 11; }
[[ -f "$GA_PREFIX/lib/libga.a" ]] || { echo "missing GA" >&2; exit 12; }
[[ -d "$WIGNERNJ_PREFIX" ]] || { echo "missing offline wignernj prefix" >&2; exit 13; }
grep -q "if (FALSE)" "$SOURCE_ROOT/CMakeLists.txt" || { echo "XMLLINT patch missing" >&2; exit 14; }

cd "$BUILD_DIR"
echo "==== repair configure $(date --iso-8601=seconds) ===="
set +e
cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -DLINALG=MKL -DHDF5=ON -DMPI=ON -DGA=ON \
  -DBUILD_SHARED_LIBS=ON \
  -DEXTERNAL_LIBXC="$OMC_PREFIX" \
  -DCMAKE_PREFIX_PATH="$WIGNERNJ_PREFIX;$GA_PREFIX;$OMC_PREFIX" \
  -DCMAKE_Fortran_COMPILER="$OMC_PREFIX/bin/mpifort" \
  -DCMAKE_C_COMPILER="$OMC_PREFIX/bin/mpicc" \
  -DCMAKE_CXX_COMPILER="$OMC_PREFIX/bin/mpicxx" \
  -DMPI_C_COMPILER="$OMC_PREFIX/bin/mpicc" \
  -DMPI_CXX_COMPILER="$OMC_PREFIX/bin/mpicxx" \
  -DMPI_Fortran_COMPILER="$OMC_PREFIX/bin/mpifort" \
  -DPYTHON_EXECUTABLE="$OMC_PREFIX/bin/python3" \
  -DEXTRA_FFLAGS="-ffree-line-length-none -Wno-error" \
  -DCMAKE_INSTALL_PREFIX="$INSTALL_PREFIX" \
  "$SOURCE_ROOT" 2>&1 | tee "$LOG_DIR/driver-fix-configure.log"
configure_rc=${PIPESTATUS[0]}
set -e
printf 'CONFIGURE_RC=%s\n' "$configure_rc" > "$LOG_DIR/driver-fix.status"
[[ "$configure_rc" -eq 0 ]] || exit 21

echo "==== repair make -j4 $(date --iso-8601=seconds) ===="
set +e
make -j4 2>&1 | tee "$LOG_DIR/driver-fix-build.log"
build_rc=${PIPESTATUS[0]}
set -e
printf 'BUILD_RC=%s\n' "$build_rc" >> "$LOG_DIR/driver-fix.status"
[[ "$build_rc" -eq 0 ]] || exit 22

echo "==== repair make install $(date --iso-8601=seconds) ===="
set +e
make install 2>&1 | tee "$LOG_DIR/driver-fix-install.log"
install_rc=${PIPESTATUS[0]}
set -e
printf 'INSTALL_RC=%s\n' "$install_rc" >> "$LOG_DIR/driver-fix.status"
[[ "$install_rc" -eq 0 ]] || exit 23

[[ -f "$INSTALL_PREFIX/data/exactemb.prgm" ]] || { echo "installed exactemb.prgm missing" >&2; exit 24; }
grep -q 'exactemb.exe' "$INSTALL_PREFIX/data/exactemb.prgm" || { echo "installed descriptor lacks executable" >&2; exit 25; }
sha256sum "$SOURCE_ROOT/src/Driver/exactemb.prgm.src" "$BUILD_DIR/data/exactemb.prgm" "$INSTALL_PREFIX/data/exactemb.prgm" | tee "$LOG_DIR/driver-fix-sha256.txt"
printf 'VERIFY_RC=0\n' >> "$LOG_DIR/driver-fix.status"
echo "==== repair complete $(date --iso-8601=seconds) ===="
