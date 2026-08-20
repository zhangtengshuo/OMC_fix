# Stage 2: fragment and dimer contexts

Status: passed (zero errors, zero repair attempts)

Started: 2026-08-20T11:11:54+08:00; executed and completed: 2026-08-20T11:22 by ZCode (takeover after the operational pause)

## Scope

Generate final state-specific A and B fragment RunFiles and a matching AB C1/RICD RunFile, ONEINT, and Cholesky context from the fixed cc-pVDZ geometry.

## Error and repair log

No runtime or build error occurred; no repair attempt was used. The three Slurm jobs were submitted after the approval-service window elapsed, using the established `run_exactemb_job.sh` contract on partition `molcas` with one MPI process and 12000 MiB.

## Execution record (2026-08-20T11:19-11:22 +08:00)

- Job 860 `omc_S2_A_neutral`, job 861 `omc_S2_B_neutral`, job 862 `omc_S2_AB_ricd`: all `COMPLETED`, all `RUN_RC=0`, elapsed 4/4/3 s, all three `slurm-*.err` files empty (0 bytes).
- All three task-owned OpenMolcas scratch directories were removed by the job script cleanup trap; `work/openmolcas/` is empty.

## Invariant verification

- S2_A_neutral: 48 basis functions, 14 closed-shell + 2 active electrons, RASSCF root-1 total energy `-78.06789017 Eh` — matches the expected reference exactly.
- S2_B_neutral: identical counts and the identical root-1 energy `-78.06789017 Eh`.
- A and B preserved RunFiles are not byte-identical (SHA-256 differ, as expected because the stored geometry differs by the rigid z translation); their scientific records agree at printed precision, which satisfies the documented fallback contract.
- S2_AB_ricd: 96 basis functions in one symmetry; preserved nonempty context files `RunFile` (358784 B), `OneInt` (2501024 B), `ChVec1` (24488736 B), `ChRed` (4157136 B), `ChRst` (31256 B), `ChMap` (3248 B), plus `RICDLib`; SEWARD reported `Cholesky decomposed two-electron repulsion integrals` under `Medium Cholesky`.
- Provisional AO starts remain A=1 and B=49 (each fragment 48 AOs, dimer 96); final confirmation is deferred to the Stage 3 EXACTEMB electron-count invariant as planned.

## Fixed context contract

- Fragment A is neutral ethylene at z=0 Angstrom; fragment B is the identical rigid ethylene translated to z=3.5 Angstrom.
- Both fragments use cc-pVDZ, `Group=NoSymm`, RHF orbitals, and state-specific CAS(2e,2o) RASSCF with `FileOrb=SCFORB`.
- The full dimer lists all A atoms before all B atoms, uses cc-pVDZ and `Group=NoSymm`, and requests GATEWAY `RICD` followed by SEWARD `Medium Cholesky`.
- Every context is regenerated from GATEWAY in its own complete task under the repaired OpenMolcas installation; no RunFile from an earlier development build is reused.
- All contexts use one MPI process, matching the RunFile parallel-environment contract.

## Expected invariants

- A and B must each report 48 basis functions, 16 electrons, and the accepted neutral root-1 energy `-78.06789017 Eh` at printed precision.
- A and B preserved RunFiles must be byte-identical if translation invariance and output metadata are both unaffected by the rigid z translation; if not byte-identical, their scientific records and energies must still agree.
- AB must report 96 basis functions in one symmetry and generate a nonempty RunFile, ONEINT, ChVec1, ChRed, ChRst, and ChMap context.
- The provisional AO starts are A=1 and B=49; Stage 3 must still verify them through the EXACTEMB electron-count invariant.
