# Stage 1 repair: RunFile short-name contract

Status: passed after the second and final runtime repair

Started: 2026-08-20T11:00:18+08:00

## Confirmed root cause

OpenMolcas declares `RunFile_data::RunName` as `character(len=8)`. `NameRun` therefore accepts a short RunFile logical name, normally resolved in the current module's work directory or through its driver descriptor; it is not an arbitrary physical-path interface.

The original EXACTEMB parser stored up to 256 characters and passed SOURCE, TARGET, and CLEAR values directly to `NameRun`. All four calls in `exactemb.F90` had the same truncation bug.

## Same-type audit

- The four affected EXACTEMB calls were the only `NameRun(trim(variable))` calls in the source tree.
- GENANO uses a six-character `RUNnnn` value.
- AVERD uses a seven-character `RUNnnn` value.
- SLAPAF `OldFcm` receives the registered names `RUNOLD` or `RUNFILE` from its callers.
- No second long-name or physical-path call was found outside EXACTEMB.

## Local repair

- Replace the 256-character source and target fields with eight-character logical-name fields.
- Read each name through a dedicated validator before assignment, rejecting values longer than eight characters and values containing `/` or a backslash.
- Register the recommended work-directory files `EMBSRC` and `EMBTGT` in `src/Driver/exactemb.prgm.src`.
- Keep SOURCE, TARGET, and CLEAR explicit in the input, but require short names such as `EMBSRC` and `EMBTGT` rather than paths or long physical filenames.
- Update the implementation document and prepare a complete GATEWAY/SEWARD/SCF/RASSCF/CLEAR input that copies its same-task current RunFile to `EMBTGT` before invoking EXACTEMB.

## Validation gates

1. Static preprocessing of `exactemb.prgm.src` must expose both `EMBSRC` and `EMBTGT` below `$WorkDir`.
2. The local diff must contain no whitespace errors or remaining long-name EXACTEMB call.
3. Rebuild and install on GV800 with the recorded `/opt/mamba/envs/omc`, MPI, GA, MKL, HDF5, LibXC, offline wignernj, CMake, and `make -j4` recipe.
4. Run the complete self-contained CLEAR input and require `RUN_RC=0`, the EXACTEMB clear diagnostic naming `EMBTGT`, and a preserved cleared RunFile.

## Rebuild result

- The GV800 incremental rebuild started at 2026-08-20T11:04:34+08:00 and completed at 2026-08-20T11:04:41+08:00.
- CMake configure, `make -j4`, install, and descriptor verification all returned zero.
- The installed driver descriptor contains `EMBSRC` and `EMBTGT`, and its SHA-256 is `ad563e9d20de1df29e3526993a6f020be31f6f67ec26d4581c101542e04d7b7e`.

## Complete-input validation attempt 1

- Slurm job 858 regenerated GATEWAY, SEWARD, SCF, and CAS(2,2) RASSCF from the complete input and reproduced the root-1 total energy `-78.06789017 Eh`.
- The EMIL step reported `COPY A_neutral_clear_shortname_fix.RunFile EMBTGT`, and EXACTEMB launched, proving that the long-path bug and module mapping blocker were bypassed.
- CLEAR then stopped at its first write with `Warning, writing temporary iScalar field`, `Field: Exact Emb Active`, `_INTERNAL_ERROR_`, and `RUN_RC=250`.
- Root cause: the implementation added eleven `Exact Emb *` records but did not register any of them in `RunFile_data`'s known-label tables; Release OpenMolcas deliberately aborts on an unsupported temporary field.

## Second and final repair for this runtime problem

Register all eleven records together rather than patching only the first failing field: six entries in `LabelsIS`, four in `LabelsDS`, and `Exact Emb Pot` in `LabelsDA`. Rebuild with the same recipe and rerun the complete input once. If that second validation fails, stop and report without a third attempt.

## Second rebuild result

- The second GV800 rebuild started at 2026-08-20T11:08:36+08:00 and completed at 2026-08-20T11:08:49+08:00.
- CMake configure, `make -j4`, install, and verification all returned zero.
- The rebuilt `runfile_data.F90` SHA-256 is `95eae85f4e7433f5a908993555e3476e407c624fb67b035cdb8ba22528db3b6e`.
- All eleven RunFile labels used by EXACTEMB and RASSCF are now present in the correct typed tables.

## Complete-input validation attempt 2

- Slurm job 859 completed with `RUN_RC=0`, scheduler state `COMPLETED`, and exit code `0:0`.
- The complete input regenerated GATEWAY, SEWARD, SCF, and CAS(2,2) RASSCF and reproduced the root-1 total energy `-78.06789017 Eh`.
- EXACTEMB printed `embedding flag cleared in target RunFile: EMBTGT`, followed by OpenMolcas `Happy landing`.
- The output contains no temporary-field warning.
- The preserved cleared RunFile is 263192 bytes with SHA-256 `6b9803b308db8cdcf5c59589c3ab0996ff17c362134ce21fa37a574f15a21937` and contains the registered Exact Emb integer-scalar, real-scalar, and `Exact Emb Pot` array label tables.
- Input SHA-256: `d861904649300d776ec302c77dc761d9e0271ee1f54b7ef2394e66a9a2463eb4`.
- Output SHA-256: `8aa2125a4bc5956750453e8b1e110208d573b0df943879f07de837fdcdaf9671`.
- No `alpha` job and no task-owned OpenMolcas scratch directory remained at 2026-08-20T11:09:27+08:00.

Completed: 2026-08-20T11:09:27+08:00
