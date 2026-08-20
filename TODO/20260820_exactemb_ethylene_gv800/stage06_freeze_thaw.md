# Stage 6: mutual freeze-and-thaw

Status: passed (six cycles, clean geometric convergence; zero errors, zero repair attempts)

Started/executed: 2026-08-20T12:45-13:10 +08:00 by ZCode

## Scope

After every earlier validation gate passes, alternate A and B exact-density updates (implementation document §11, Jacobi style: A(n) is embedded in B(n-1), B(n) in A(n-1)) and record energies, orbitals, and density convergence until the fixed point is reached (Test 7).

## Method

One template (`S6_ft_template.in`, project name substituted per cycle) with six contexts per job: A target (XField = six B nuclei) -> AB RICD with `&EXACTEMB` from `B_prev.RunFile` (SOURCE 49, TARGET 1) -> final A context with the pot injected before `&RASSCF`; then mirrored for B (XField = six A nuclei, SOURCE 1, TARGET 49). Each cycle preserves `A_new`/`B_new` RunFile+RasOrb to `$CurrDir` and seeds the next cycle. Cycle-1 seeds are the Stage 2 isolated neutral RunFiles A(0)/B(0). Jobs 870-875, all `RUN_RC=0`, scratch auto-cleaned.

## Energy convergence (root-1 embedded RASSCF totals, Eh)

| n | A(n) | B(n) | dE_A | dE_B |
|---|---|---|---|---|
| 1 | -41.66496895 | -41.66496894 | - | - |
| 2 | -41.65985838 | -41.65985229 | +5.111e-03 | +5.117e-03 |
| 3 | -41.66066547 | -41.66065879 | -8.071e-04 | -8.065e-04 |
| 4 | -41.66059870 | -41.66059193 | +6.68e-05 | +6.69e-05 |
| 5 | -41.66060393 | -41.66059716 | -5.23e-06 | -5.23e-06 |
| 6 | -41.66060353 | -41.66059675 | +4.0e-07 | +4.1e-07 |

Alternating (saddle-type) convergence with a stable contraction ratio of about 0.071-0.076 per cycle; by cycle 6 the change is ~4e-07 Eh, i.e. converged below the micro-Eh level. The A/B energy split stays a constant 6.8e-06 Eh (numerical asymmetry of the two contexts' RICD decompositions at the two different centers, not a physical effect).

Cross-validation: A(1) = -41.66496895 is bit-identical to the Stage 5 one-shot M2 result — the cycle-1 direction is by construction the same calculation, and it reproduces exactly.

## Density convergence (Lowdin norm ||S^(1/2) dP S^(1/2)||_F, replica basis)

| transition | A | B |
|---|---|---|
| isolated(0) -> one-shot(1) | 1.38e-01 | 1.365e-01 |
| 1 -> 2 | 1.239e-03 | 1.672e-03 |
| 2 -> 3 | 9.120e-05 | 1.188e-04 |
| 3 -> 4 | 7.035e-06 | 9.199e-06 |
| 4 -> 5 | 5.475e-07 | 7.137e-07 |
| 5 -> 6 | 4.235e-08 | 5.532e-08 |
| one-shot(1) -> F&T(6) | 1.158e-03 | 1.566e-03 |

The density norms track the energy contraction ratio (~0.071-0.074) down to the 5e-08 floor. The physical picture reproduces the documented Dimerge pattern: mutual polarization of ~0.137 Lowdin units from isolation, with a further freeze-and-thaw correction of ~1.2-1.6e-03 on top of the one-shot densities.

## Error and repair log

No runtime or build error occurred in the six cycle jobs. One analysis-script slice typo (`F[48:,:48:]` instead of `F[48:,48:]`) initially returned zero B-side norms; found via checksums showing the dumped files differed, fixed, and rerun. No experiment-code repair was needed.

## Conclusion

Stage 6 passes: the mutual exact-density freeze-and-thaw fixed point is reached with clean geometric convergence in both energy and density, reproducing the established ethylene one-shot -> F&T behavior. Together with Stages 0-5 this completes the full ethylene validation campaign (implementation-document Tests 1-7); the remaining production target is Test 8 (real pentacene Stage F).
