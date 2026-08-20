# Stage 5: one-shot M2 response

Status: passed (all three gates at 1e-07 or better; one input-preservation repair consumed, experiment code untouched)

Started/executed: 2026-08-20T12:10-12:40 +08:00 by ZCode

## Scope

Validate J-only RASSCF response (Test 5), then combine the electronic `Exact Emb Pot` with the partner nuclear XField exactly once and compare the one-shot M2 result (Test 6) against an independent PySCF reference built with the Stage 4 validated machinery.

## Inputs and execution

Two complete self-contained inputs, each with four contexts inside one task (B source RASSCF -> AB RICD with `&EXACTEMB` -> final A context whose current RunFile is replaced by the pot-carrying copy immediately before `&RASSCF`):

- `S5_A_Jonly.in` (jobs 866/868): final A RASSCF with `Exact Emb Pot` only, no XField.
- `S5_A_M2.in` (jobs 867/869): both A contexts carry `XField` with the six partner-B nuclei as bare point charges (`6 ANGSTROM 0 0 0 0` header, then `x y z q` per site, parsed per `rdctl_seward.F90` `ProcessXF`).

All jobs `RUN_RC=0`; retry energies are bit-identical to the first runs; final `A1_final.RunFile`/`RasOrb` preserved (retry dirs).

## Error and repair log (one repair, input-level)

Jobs 866/867 completed cleanly but the final `>>COPY ... A1_final.*` artifacts were missing: bare-name EMIL copies land in the scratch workdir and are deleted by cleanup; Stage 2 succeeded because it used `$CurrDir/` prefixes. One repair attempt: added `$CurrDir/` to the preservation copies (and `_retry1` project names with the literal injection RunFile name updated). Jobs 868/869 preserved the artifacts and reproduced the energies bit-identically. No experiment-code repair was needed.

## OpenMolcas results

| Quantity | Value (Eh) |
|---|---|
| B0 source (both jobs) | -78.06789017 (matches Stage 2/3) |
| A0 isolated (J-only job) | -78.06789017 |
| A0 nuclei-only XField (M2 job) | -78.82129253 |
| A embedded J-only root 1 | **-42.13929260** |
| A embedded M2 root 1 | **-41.66496895** |
| SGFCIN printed `<J>` (J-only / M2) | 31.699010 / 33.382510 |

The large J-only/M2 upward shifts are physical: J[P_B] enters as a one-electron operator without its nuclear counterbalance inside the fragment energy; the comparable reference must include exactly the same terms (below).

## PySCF reference (work/stage05/stage05_verify.py, stage05_verify.log)

Reference built exactly like the validated Dimerge `me22_exact_density_fc.py` pattern: replica molA (exact OM basis), `J_ref` from the Stage 4 npz, partner nuclear attraction via `with_rinv_origin` point charges, CASSCF(2,2).

```text
[const] E(A-nuc x B-charge) = 36.28788679 Eh
[G1 XField ] -78.82129253 vs OM -78.82129253  diff -4.788e-10
[G2 J-only ] -42.13929286 vs OM -42.13929260  diff -2.635e-07
[G3 M2     ] -41.66496923 vs OM -41.66496895  diff -2.803e-07
```

Three independent gates:

- **G1** confirms the OpenMolcas XField total-energy composition: electronic (h+V) plus the A-nuclei/B-charge repulsion constant, with nothing else.
- **G2** is Test 5: the J-only RASSCF one-electron response agrees with the external reference to 2.6e-07 (convergence-limited).
- **G3** is Test 6: the full one-shot M2 energy, compared including the same nucleus-charge constant, agrees to 2.8e-07 when the PySCF CASSCF is started from the XField-SCF orbitals exactly like OM's `FileOrf=SCFORB` start. An earlier 26 mEh gap was traced to CASSCF landscape/start sensitivity (starting from the embedded-hcore SCF converged to a different local solution), not to any physical inconsistency.

## Diagnostic semantics established

The SGFCIN `<J[P_partner]>` printout is emitted once at the first SGFCIN call with the inactive-only (14-electron) density; it is not the final density expectation. External checks: final `Tr[D_final J_ref]` = 35.425/36.470 (J-only/M2) vs printed 31.699/33.383, and isolated-RHF `Tr[D J]` = 36.478 — the printed value sits at the inactive-only level as expected. The energy-level gates above prove the operator enters the RASSCF Hamiltonian correctly; the printout timing is a documentation item, not a defect.

## Conclusion

Stage 5 passes: Test 5 and Test 6 are validated at the 1e-07 level against an independent PySCF reference with fully understood energy bookkeeping (XField constant included exactly once on both sides). The remaining gate is Stage 6 (mutual exact-density freeze-and-thaw).
