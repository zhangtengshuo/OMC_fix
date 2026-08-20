# Stage 1: runtime smoke tests

Status: stopped after two complete-input CLEAR attempts

Started: 2026-08-20T10:36:10+08:00

## Scope

Run an ordinary RASSCF with no exact-density record, verify CLEAR handling, and prove that the new executable can open compatible RunFiles without changing an unembedded baseline.

## Error and repair log

- Attempt 1, Slurm job 852, failed before GATEWAY started.
- The preserved output reports `parnell_base: cannot make directory` below `/scratch/alpha/exactemb-ethylene-20260820_103419/work/openmolcas/A_neutral_baseline.852/A_neutral_baseline` and `RUN_RC=1`.
- Root cause: the custom runner created only `MOLCAS_SCRATCH_ROOT`, while explicit `MOLCAS_WORKDIR` named a deeper parent directory that did not yet exist; Parnell attempted to create the project subdirectory below that missing parent.
- Single repair: change the runner setup from `mkdir -p "$MOLCAS_SCRATCH_ROOT"` to `mkdir -p "$MOLCAS_WORKDIR"` and submit once under a new project name so the failed output remains immutable.
- Stop condition: any failure of the repaired Stage 1 submission ends Stage 1 without another repair.

## Repaired baseline result

- Slurm job 853 completed with `RUN_RC=0` and OpenMolcas `Happy landing`.
- The final neutral-A RunFile size is 263192 bytes, and RasOrb and VecDet were preserved.
- The final RASSCF root-1 energy is `-78.06789017 Eh`, identical at the printed precision to the accepted historical Round 2 neutral-A energy.
- No `Exact Emb` activation diagnostic appears in the output.
- The task-owned `/scratch/alpha/exactemb-ethylene-20260820_103419/work/openmolcas/A_neutral_baseline_retry1.853` tree was absent after completion, confirming cleanup.
- The repaired runner is accepted; no further repair is allowed in Stage 1.

## CLEAR result and mandatory stop

- Slurm job 854 copied the byte-identical target RunFile into task scratch, parsed the input, and then stopped before launching `exactemb.exe`.
- The exact pymolcas diagnostic is `Requested module not available`, the process status is `RUN_RC=36`, and the internal return condition is `_RC_NOT_AVAILABLE_`.
- Input SHA-256: `5c809800e4f657d53a2b09390ab7bacb432cf399ddb0a901357413f61c405a99`.
- Output SHA-256: `d3fd1fdc79a2ba36c06d51c2104c8eea0758d950feae15141a2219dc6759df5a`.
- Untouched target RunFile SHA-256: `8320da18ca69e7ed8bc75db21232eaf8f1fdea05f0caaba94652b99441b8371a`.
- The task-owned scratch tree was cleaned and no `alpha` job remained queued or running at 2026-08-20T10:41:01+08:00.

## Read-only root-cause diagnosis at the job 854 stop

- `Tools/pymolcas/molcas_wrapper.py` resolves every `&NAME` module by opening `$MOLCAS/data/name.prgm`; absence of that file raises `Unknown module` and maps to `_RC_NOT_AVAILABLE_`.
- The installed tree contains `/scratch/alpha/OpenMolcasUP_ExactDenEmb-65a3eeb11/data/expbas.prgm` but does not contain `data/exactemb.prgm`.
- At the time of job 854, the local branch contained no `src/Driver/exactemb.prgm.src`.
- Top-level CMake auto-discovery built and installed `exactemb.exe`, but module executability through pymolcas requires a separate driver descriptor generated from `src/Driver/*.prgm.src`.
- The build record's executable and `ldd` checks therefore proved linking but did not prove `&EXACTEMB` driver registration.

## Repair proposed at the job 854 stop

The proposed repair was to add `src/Driver/exactemb.prgm.src` with the `exactemb.exe` executable mapping, EXACTEMB input unit, ONEINT, Cholesky files, and standard global RunFile mapping, then regenerate the CMake `prgms` target, build with `make -j4`, install, and verify that `data/exactemb.prgm` exists. It was not attempted until the user explicitly authorized continuation; its subsequent successful execution is recorded in `stage01_driver_fix.md`.

Stopped: 2026-08-20T10:41:01+08:00

## Authorized continuation with complete input

The experiment resumed after the user required complete inputs and allowed two attempts per distinct future problem. The new CLEAR validation regenerates the neutral-A GATEWAY, SEWARD, SCF, RASSCF, and current project RunFile in one job before invoking EXACTEMB, so it has no RunFile dependency on job 853 or job 855.

### Complete-input attempt 1

- Slurm job 856 regenerated neutral-A GATEWAY, SEWARD, SCF, and CAS(2,2) RASSCF successfully and again obtained 16 electrons.
- EXACTEMB was recognized and launched, but `NameRun` could not open the relative target name `A_neutral_clear_full.RunFile`; the job ended with `RunFile does not exist` and `RUN_RC=250`.
- The current project RunFile was sufficient for module startup; the failure occurs when the CLEAR branch opens the separately named target argument.

### Complete-input attempt 2

- Slurm job 857 regenerated the complete neutral-A GATEWAY, SEWARD, SCF, and CAS(2,2) RASSCF context and again obtained the root-1 total energy `-78.06789017 Eh`.
- The input then used `>>COPY $Project.RunFile $CurrDir/A_neutral_clear_target.RunFile` to create the CLEAR target from the same task immediately before EXACTEMB; no artifact from another task was used.
- The snapshot exists, is 263192 bytes, and has SHA-256 `d6a4be3c03400459c848ce1d6ab44e65e0c8a712671148bf657b8bd1be67eb98`.
- EXACTEMB was recognized and launched, but the absolute target path was reported by `AixOpn` as `/SCRATCH`; the job ended with `Premature abort while opening file No such file or directory` and `RUN_RC=250`.
- Input SHA-256: `7cb2699eabe6733e7131842a27a72cb872f37614c6133ae31146fdaed052c327`.
- Output SHA-256: `7d8c0b34cb4ae9a9ed1854ae4d2e614e62e6ac19713b327e8ca42b2fc7b746df`.
- No `alpha` job remained queued or running, and no task-owned OpenMolcas scratch directory remained at 2026-08-20T10:53:06+08:00.

## Root cause of the two complete-input failures

- `src/runfile_util/runfile_data.F90` declares both `RunName` and its stack as `character(len=8)`.
- `NameRun` merely assigns its argument to this 8-character variable and clears the RunFile cache; it does not retain an arbitrary physical filename.
- RunFile readers call `f_inquire(RunName)`, which passes this short name through the current module's driver-file translation table.
- Consequently, attempt 1 treated `A_neutral_clear_full.RunFile` as an unregistered, truncated logical name and could not find a RunFile; attempt 2 truncated the absolute `/scratch/alpha/...` path to `/scratch`, exactly matching the `/SCRATCH` path printed by `AixOpn`.
- Existing OpenMolcas callers consistently use registered logical names such as `RUNFILE`, `AUXRFIL`, and `RUNOLD` with `NameRun`.
- The new EXACTEMB implementation and its input syntax are therefore inconsistent with the RunFile API: `SOURCE`, `TARGET`, and `CLEAR` currently accept 256-character physical filenames, but `exactemb.F90` passes those strings directly to an 8-character logical-name interface.

## Required code repair for a future continuation

Add two distinct, at-most-eight-character logical RunFile entries to `src/Driver/exactemb.prgm.src`, switch `exactemb.F90` to those fixed logical names, and revise the input/workflow so complete jobs populate the corresponding physical files explicitly. The driver mapping, parser contract, examples, and implementation documentation must be changed together; merely shortening or relocating the current filename is not a correct repair.

The two permitted attempts for this distinct problem are exhausted. No third CLEAR task and no Stage 2 task was submitted.

Stopped: 2026-08-20T10:53:06+08:00

## Subsequent authorized repair

The user authorized a code repair after this stop. Short RunFile logical names, driver mappings, parser guards, and all eleven typed Exact Emb RunFile labels were added; the complete job 859 CLEAR validation passed. See `stage01_runfile_name_fix.md` for the repair and rebuild evidence.
