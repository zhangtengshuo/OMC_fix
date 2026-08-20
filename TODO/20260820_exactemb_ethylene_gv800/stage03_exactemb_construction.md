# Stage 3: EXACTEMB construction

Status: passed (zero errors, zero repair attempts)

Started/executed: 2026-08-20T11:24-11:27 +08:00 by ZCode

## Scope

Construct B to A and A to B electronic Coulomb potentials, validate mapped electron counts and record metadata, and deliberately exercise one AO-mapping failure guard.

## Error and repair log

No runtime error occurred in the two construction jobs, and the guard job aborted exactly as designed; no repair attempt was used.

## Execution record

Three complete self-contained inputs (regenerated full-dimer GATEWAY RICD + SEWARD Medium Cholesky context in the same task, then `>>COPY` of the Stage 2 fragment RunFiles to `EMBSRC`/`EMBTGT`, then `&EXACTEMB`), submitted on partition `molcas`, one MPI process, 12000 MiB:

- Job 863 `omc_S3_BtoA` (SOURCE EMBSRC 49, TARGET EMBTGT 1): `COMPLETED`, `RUN_RC=0`, elapsed 3 s.
- Job 864 `omc_S3_AtoB` (SOURCE EMBSRC 1, TARGET EMBTGT 49): `COMPLETED`, `RUN_RC=0`, elapsed 3 s.
- Job 865 `omc_S3_guard_srcrange` (SOURCE start deliberately 73, outside the 96-AO dimer): expected hard abort; Slurm `FAILED` with `RUN_RC=250` is the designed outcome.

All three task-owned scratch directories were removed by the job script cleanup trap.

## Construction verification (B -> A, job 863)

Verbatim summary block:

```text
  Source RunFile        : EMBSRC
  Target RunFile        : EMBTGT
  Source relaxation root:                     1
  Dimer AO count        :                    96
  Source AO block       :                    49                   96
  Target AO block       :                     1                   48
 Expected electrons    :        16.0000000000
 Tr[P S]               :        16.0000000000
 |Delta N|             :   1.705303E-13
 packed ||J||_2        :   2.087182E+01
 full-matrix ||J||_F   :   2.505784E+01
 max |J_mn|            :   2.448000E+00
  Exchange contribution : DISABLED
  Partner nuclei        : supply with standard SEWARD XField
  Target RunFile record : Exact Emb Pot
```

The electron-count invariant `|Tr[P S] - N_e| = 1.7e-13` passes the `TOLERANCE=1.0d-6` gate with nine orders of magnitude of margin, confirming the Stage 2 AO starts A=1, B=49.

## Construction verification (A -> B, job 864)

Same block with mirrored ranges (`Source AO block 1 48`, `Target AO block 49 96`); `|Delta N| = 1.776357E-13`. Because A and B are identical rigid fragments, the J statistics must match the B -> A case exactly, and they do: `packed ||J||_2 = 2.087182E+01`, `full-matrix ||J||_F = 2.505784E+01`, `max |J_mn| = 2.448000E+00`. This mirror equality is an additional physical consistency check.

## Guard verification (job 865)

Verbatim abort section:

```text
 EXACTEMB: SOURCE AO block lies outside the dimer AO range.
[ process      0]: xquit (rc =    128): _INTERNAL_ERROR_
 SOURCE =                    73                  120  DIMER = 1,                   96
Program aborted. Backtrace: ...
```

The guard names the violation, prints the offending block 73..120 against the dimer range 1..96, and hard-aborts (`rc=-6` SIGABRT via standard OpenMolcas Abend; pymolcas `RUN_RC=250`). No `Exact Emb Pot` was written for the guard case.

## Preserved artifacts

- `work/stage03/S3_BtoA/S3_EMBTGT_BtoA.RunFile` and `work/stage03/S3_AtoB/S3_EMBTGT_AtoB.RunFile`: 273088 bytes each, versus 263192 bytes for the unmodified fragment RunFile seed. The 9896-byte increment equals 9408 bytes (48*49/2 real*8 elements = packed 48-AO lower-triangular J) plus record overhead, confirming `Exact Emb Pot` was physically written.
- Full job outputs, provenance, and status files remain in the respective `work/stage03/*/` run directories; `.err` files are empty for jobs 863/864.

## Conclusion

Stage 3 passes: both directed constructions satisfy the electron-count invariant at machine precision with mirror-consistent J statistics, and the AO-range guard fires with the correct diagnostic and a clean hard abort. The next gate is Stage 4 (PySCF direct-J numerical comparison of the target J block).
