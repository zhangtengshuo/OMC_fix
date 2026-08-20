# GV800 ethylene EXACTEMB validation

Date: 2026-08-20

Remote host: `GV800`

Remote user: `alpha`

Remote root: `/scratch/alpha/exactemb-ethylene-20260820_103419`

OpenMolcas installation: `/scratch/alpha/OpenMolcasUP_ExactDenEmb-65a3eeb11`

The experiment validates the new pure-OpenMolcas exact-density electronic embedding implementation against the established Dimerge ethylene dimer references.

## Error policy

For each distinct blocking compile or runtime problem encountered after 2026-08-20T10:46:40+08:00, record the original evidence and make at most two reasoned attempts. If the second attempt fails, stop that route and report.

Prefer a complete, self-contained OpenMolcas input that starts from GATEWAY/SEWARD and generates every required current RunFile and integral context inside the same task. Do not depend on a RunFile copied from another task when the same context can be regenerated cheaply; any same-task intermediate context switch must be explicit and checksummed.

If a source repair is required, rebuild incrementally with the GV800 recipe recorded in `TODO/20260820_GV800_exactemb_build_record.md`: the `/opt/mamba/envs/omc` compiler and MPI stack, MKL, GA, HDF5, LibXC, the offline system wignernj installation, the same CMake cache, and `make -j4`.

## Stages

| Stage | Purpose | Record | Status |
|---|---|---|---|
| 0 | Environment, provenance, remote layout, and reusable-data audit | `stage00_setup.md` | passed |
| 1 | No-active regression and CLEAR/runtime smoke tests | `stage01_runtime_smoke.md`, `stage01_runfile_name_fix.md` | passed after RunFile interface repair |
| 2 | Generate consistent A/B fragment RunFiles and AB C1/RICD context | `stage02_contexts.md` | passed (jobs 860-862, zero errors) |
| 3 | Construct B to A and A to B potentials and exercise validation guards | `stage03_exactemb_construction.md` | passed (jobs 863-865, zero errors) |
| 4 | Compare the OpenMolcas target J block with a PySCF direct-J reference | `stage04_direct_j.md` | passed (elementwise rel 1.8e-08) |
| 5 | Validate J-only and full one-shot M2 RASSCF response | `stage05_oneshot_m2.md` | passed (three gates at 1e-07; one input repair) |
| 6 | Run mutual exact-density freeze-and-thaw after all preceding gates pass | `stage06_freeze_thaw.md` | passed (6 cycles, geometric convergence; zero errors) |

## Stage 1 failure history and resolution

The driver-registration repair recorded in `stage01_driver_fix.md` succeeded: CMake generated and installed `data/exactemb.prgm`, and subsequent jobs entered `Start Module: exactemb`.

The user then authorized continuation with complete inputs and up to two attempts per distinct problem. Jobs 856 and 857 each regenerated the full neutral-A GATEWAY, SEWARD, SCF, and RASSCF context in one task before CLEAR. Job 857 additionally copied the current task's RunFile immediately before EXACTEMB, proving that the target artifact existed and was internally consistent.

Both attempts failed because the code passes the user-supplied physical target filename to `NameRun`, while OpenMolcas stores `RunName` in an eight-character field intended for driver-registered logical names. The relative filename in job 856 could not be resolved, and the absolute `/scratch/alpha/...` filename in job 857 was truncated to `/scratch`, which appeared as `/SCRATCH` in the `AixOpn` failure. This is a code/API contract bug rather than a missing-input-file problem.

The two permitted attempts are exhausted. No third CLEAR task and no Stage 2 task was submitted. A future repair must add fixed source and target logical RunFile entries to the EXACTEMB driver descriptor and update the Fortran input semantics, documentation, and complete-input workflow together.

The user subsequently authorized this repair. The implementation now accepts only validated RunFile logical names of at most eight characters with no path separator, and the driver registers `EMBSRC` and `EMBTGT` in `$WorkDir`.

The first post-repair complete-input test, job 858, reached CLEAR but exposed a related omission: all eleven new `Exact Emb *` records were absent from the typed RunFile label tables and Release OpenMolcas rejected the first as a temporary field. The second and final repair registered all six integer scalars, four real scalars, and one real array together. Job 859 then completed with `RUN_RC=0`, printed the `EMBTGT` clear diagnostic and `Happy landing`, and left no job or task scratch. Stage 1 is now passed; the full record is `stage01_runfile_name_fix.md`.
