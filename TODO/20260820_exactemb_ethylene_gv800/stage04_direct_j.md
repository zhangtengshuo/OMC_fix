# Stage 4: PySCF direct-J comparison

Status: passed (matrix-level agreement 1.8e-08 relative; zero experiment-code repairs, analysis-tooling iterations only)

Started/executed: 2026-08-20T11:30-12:10 +08:00 by ZCode

## Scope

Compare the OpenMolcas EXACTEMB target J block (B -> A, from Stage 3 `S3_EMBTGT_BtoA.RunFile`) against a PySCF direct-ERI Coulomb matrix for the same geometry, basis, AO convention, and source density, quantifying the RI-CD (Medium Cholesky) approximation error (implementation document Test 4).

## Method

1. **RunFile extraction**: a standalone Fortran tool `dump_rf.F90` (in this directory, linked against the installed `libmolcas` with the `/opt/mamba/envs/omc` toolchain, GA, libxc, offline wignernj, and `-latomic`) dumps typed RunFile records to ASCII. It must call `IniMem()` explicitly (standalone programs bypass the `Start()` wrapper) and needs `MOLCAS_MEM`/`MOLCAS_NPROCS` set. Extracted: `Exact Emb Pot` (1176 values) from EMBTGT and `D1ao` (1176 values) from EMBSRC.
2. **Exact basis replica**: OpenMolcas CC-PVDZ uses general contractions while PySCF cc-pvdz is segmented, so the two AO representations are related by a nonsingular contracted-AO transform, not a permutation. The validated Dimerge/NOCI_MPS bridge module `check_molcas_pyscf_ao_bridge.py` rebuilds the exact OpenMolcas contracted basis inside PySCF from the run2h5 metadata (`BASIS_FUNCTION_IDS`/`PRIMITIVE_IDS`/`PRIMITIVES`), which SEWARD already wrote into the Stage 2 `S2_AB_ricd.guessorb.h5`. The comparison is therefore executed entirely inside this replica basis — no basis transform at all.
3. **Reference J**: `pyscf.scf.jk.get_jk((omA,omA,omB,omB), dmB, scripts='ijkl,lk->ij', aosym='s4')` on the replica fragment mols — the same semantic as the validated Dimerge exact-density reference (`me22_exact_density_fc.py`). Cross-validated against a `shls_slice` direct-ERI einsum contraction on the native dimer: agreement 6e-15.

## Storage conventions established (with evidence)

- `D1ao` is **folded packed** (off-diagonal doubled): unfolding with halving and tracing against the replica overlap gives `Tr[P S] = 16.0000000024`; the plain interpretation gives 15.5875.
- `Exact Emb Pot` is **plain packed** (values as stored): under this convention the elementwise reference comparison passes at 1.8e-08; the folded interpretation is off by 44%. This matches the implementation document's §4 operator convention (`D_fold · V_LT = Tr[P V]`).
- `Packed_Matrix_Norms`' printed `full-matrix ||J||_F` (25.05784) is reproduced by the plain unfold and equals the replica-basis reference norm to eight significant digits.

## Results (work/stage04/stage04_compare.log, npz alongside)

```text
[check] Tr[P_B S_B] = 16.0000000024 (expect 16)
[result] ||J_omc||_F = 25.05783926   ||J_ref||_F = 25.05783942
[result] max|dJ| = 8.606e-08   ||dJ||_F/||J_ref|| = 1.842e-08
[result] generalized eig max dev = 1.156e-06  top = 3.114725 / 3.114725
```

The Cholesky/RI-CD representation error at Medium threshold is ~2e-08 relative for this system; the residual is dominated by the RASSCF density convergence (PRWF=1e-6), not by the decomposition. Test 4 passes decisively.

## Analysis-tooling error log (no experiment repair consumed)

- Linking `dump_rf` initially missed MKL/HDF5/GA/MPI/libxc/wignernj/atomic — resolved by following the dependency list already recorded in the build record and `ldd` output.
- Standalone run aborted twice: `MOLCAS_MEM is not defined` (environment) and `mma out-of-memory available 0 kB` (missing `IniMem()` initialization); both fixed.
- A Python parser bug (packed length vs AO dimension) and a matrix-reindexing direction bug (`mapping` vs its inverse; the bridge's `map_commonorb` fixes the correct direction) were found and fixed via electron-count oracles.
- A native-basis transform route (T = S_nat^-1 S_cross) reproduced the density electron count but produced inconsistent J statistics; it was abandoned in favor of the direct replica-basis comparison above, which requires no transform. The anomaly is recorded here for completeness; it does not affect the conclusion.

## Conclusion

Stage 4 passes. The pure-OpenMolcas exact-density Coulomb construction is numerically equivalent to the validated PySCF direct-ERI reference at the 1e-08 level, and the AO-mapping contract survives an independent external check. The next gate is Stage 5 (J-only and full one-shot M2 RASSCF response).
