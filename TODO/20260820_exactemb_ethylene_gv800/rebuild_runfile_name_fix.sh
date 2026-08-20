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
LOG_STEM=${1:-runfile-name-fix}

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
grep -q "EMBSRC" "$SOURCE_ROOT/src/Driver/exactemb.prgm.src" || { echo "driver descriptor lacks EMBSRC" >&2; exit 15; }
grep -q "EMBTGT" "$SOURCE_ROOT/src/Driver/exactemb.prgm.src" || { echo "driver descriptor lacks EMBTGT" >&2; exit 16; }
grep -q "Exact Emb Active" "$SOURCE_ROOT/src/runfile_util/runfile_data.F90" || { echo "RunFile labels lack Exact Emb Active" >&2; exit 17; }
grep -q "Exact Emb Pot" "$SOURCE_ROOT/src/runfile_util/runfile_data.F90" || { echo "RunFile labels lack Exact Emb Pot" >&2; exit 18; }

cd "$BUILD_DIR"
echo "==== RunFile-name repair configure $(date --iso-8601=seconds) ===="
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
  "$SOURCE_ROOT" 2>&1 | tee "$LOG_DIR/$LOG_STEM-configure.log"
configure_rc=${PIPESTATUS[0]}
set -e
printf 'CONFIGURE_RC=%s\n' "$configure_rc" > "$LOG_DIR/$LOG_STEM.status"
[[ "$configure_rc" -eq 0 ]] || exit 21

echo "==== RunFile-name repair make -j4 $(date --iso-8601=seconds) ===="
set +e
make -j4 2>&1 | tee "$LOG_DIR/$LOG_STEM-build.log"
build_rc=${PIPESTATUS[0]}
set -e
printf 'BUILD_RC=%s\n' "$build_rc" >> "$LOG_DIR/$LOG_STEM.status"
[[ "$build_rc" -eq 0 ]] || exit 22

echo "==== RunFile-name repair make install $(date --iso-8601=seconds) ===="
set +e
make install 2>&1 | tee "$LOG_DIR/$LOG_STEM-install.log"
install_rc=${PIPESTATUS[0]}
set -e
printf 'INSTALL_RC=%s\n' "$install_rc" >> "$LOG_DIR/$LOG_STEM.status"
[[ "$install_rc" -eq 0 ]] || exit 23

[[ -f "$INSTALL_PREFIX/data/exactemb.prgm" ]] || { echo "installed exactemb.prgm missing" >&2; exit 24; }
grep -q 'EMBSRC' "$INSTALL_PREFIX/data/exactemb.prgm" || { echo "installed descriptor lacks EMBSRC" >&2; exit 25; }
grep -q 'EMBTGT' "$INSTALL_PREFIX/data/exactemb.prgm" || { echo "installed descriptor lacks EMBTGT" >&2; exit 26; }
sha256sum "$SOURCE_ROOT/src/exactemb/exactemb.F90" "$SOURCE_ROOT/src/exactemb/exactemb_data.F90" "$SOURCE_ROOT/src/exactemb/rdinp_exactemb.F90" "$SOURCE_ROOT/src/runfile_util/runfile_data.F90" "$SOURCE_ROOT/src/Driver/exactemb.prgm.src" "$BUILD_DIR/data/exactemb.prgm" "$INSTALL_PREFIX/data/exactemb.prgm" | tee "$LOG_DIR/$LOG_STEM-sha256.txt"
printf 'VERIFY_RC=0\n' >> "$LOG_DIR/$LOG_STEM.status"
echo "==== RunFile-name repair complete $(date --iso-8601=seconds) ===="
