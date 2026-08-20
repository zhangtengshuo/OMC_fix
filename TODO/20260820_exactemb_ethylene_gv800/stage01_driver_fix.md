# Stage 1 repair: pymolcas EXACTEMB driver registration

Status: original registration bug fixed; stopped on separate CLEAR-context blocker

Started: 2026-08-20T10:41:01+08:00

## Trigger

The preserved Stage 1 job 854 output returned `Requested module not available` and `RUN_RC=36` before launching `exactemb.exe`.

## Root cause

pymolcas resolves `&NAME` through `$MOLCAS/data/name.prgm`. CMake built and installed `exactemb.exe`, but the source tree had no `src/Driver/exactemb.prgm.src`, so the generated and installed `data/exactemb.prgm` descriptor did not exist.

## Repair design

- Add `src/Driver/exactemb.prgm.src`.
- Map the module to `$MOLCAS/bin/exactemb.exe`.
- Map the input to logical name `EXACTINP`; `SpoolInp` derives this name from the first five characters of program name `exactemb`.
- Include `aoints.inc` for ONEINT and ORDINT mappings.
- Include `choints.inc` for the Cholesky files consumed by `CHO_X_INIT` and `CHORAS_DRV`.
- Rely on `global.prgm` for the standard current RunFile mapping, as pymolcas merges every module descriptor with that global descriptor.

## Build contract

Use the existing GV800 `openmolcas-build-j4` CMake tree and original `/opt/mamba/envs/omc` environment, offline wignernj prefix, MKL, GA, HDF5, LibXC, MPI, and install prefix. Run CMake only if required to discover the new `src/Driver/*.prgm.src`, then run `make -j4` and `make install` exactly as recorded in `TODO/20260820_GV800_exactemb_build_record.md`.

## Acceptance gates

1. Preprocessing the descriptor must produce an executable line, EXACTINP, ONEINT, and Cholesky logical files without cpp residue.
2. GV800 build and installation must complete successfully.
3. The installed `data/exactemb.prgm` must exist and name `exactemb.exe`.
4. The preserved CLEAR smoke input must reach `exactemb.exe`, print the clear diagnostic, return zero, and create a modified RunFile copy.

## Error log

No repair-build error has occurred.

## Build result

- CMake completed at 2026-08-20T10:45:13+08:00 with the original Release, MKL, HDF5, MPI, GA, LibXC, offline wignernj, compiler, and installation settings.
- The CMake output again confirmed system libwignernj 0.6.0 from `/scratch/alpha/build-exactden-20260820_093910/wignernj-sys`.
- `make -j4` completed and generated `data/exactemb.prgm` through the `prgms` target.
- Existing Fortran object and executable targets were already up to date; the repair changes only driver metadata.
- `make install` completed and installed `/scratch/alpha/OpenMolcasUP_ExactDenEmb-65a3eeb11/data/exactemb.prgm`.
- Source descriptor SHA-256: `e6bee4d309f3b900ed30fa186c62069479818f0f858261be54d657c25ef18eca`.
- Both build-tree and installed descriptor SHA-256 values are `015fbb5ea91c9957aa092bf651b836965ad838a6a9e93e891b0bb2201d6d5da3`.
- Configure, build, installation, checksum, and status logs are preserved below `/scratch/alpha/exactemb-ethylene-20260820_103419/logs` without replacing the original 09:39 build logs.

## Runtime verification policy

Run one CLEAR smoke test against a fresh copy of the accepted neutral-A RunFile. If this verification fails, stop without a second driver adjustment.

## Runtime verification result

- Slurm job 855 no longer reported `Requested module not available`.
- pymolcas printed `Start Module: exactemb`, launched `exactemb.exe`, and printed the EXACTEMB module banner; this directly verifies that the original driver-registration bug is fixed.
- The job then stopped inside `Start('exactemb')` with `RunFile does not exist`, before `ExactEmb()` reached its CLEAR branch.
- Root cause of this separate failure: every module receives the global current-RunFile mapping `$WorkDir/$Project.RunFile`, but the smoke input copied only `A_neutral_clear_target.RunFile` under its alternate target name. A normal production call after full-dimer GATEWAY/SEWARD already has the required current project RunFile.
- A corrected standalone CLEAR test would also copy a compatible RunFile to `$Project.RunFile`, or would execute CLEAR in an already initialized current OpenMolcas context.
- No input correction or second run was attempted because this verification failure occurred after the one authorized driver repair.
- Installed descriptor SHA-256: `015fbb5ea91c9957aa092bf651b836965ad838a6a9e93e891b0bb2201d6d5da3`.
- Job-855 input SHA-256: `55cd9893f640b011ed45c7acf6128495c546560d19a027d6922389b96482e2e5`.
- Job-855 output SHA-256: `cf84ac73653fa71347e851acdd496ebdf60486870d6c5be7258e0d9876474df3`.
- The target RunFile remained byte-identical with SHA-256 `8320da18ca69e7ed8bc75db21232eaf8f1fdea05f0caaba94652b99441b8371a`.
- Build status is `CONFIGURE_RC=0`, `BUILD_RC=0`, `INSTALL_RC=0`, and `VERIFY_RC=0`.
- No `alpha` job or task-owned scratch remained at 2026-08-20T10:46:40+08:00.

Completed/stopped: 2026-08-20T10:46:40+08:00
