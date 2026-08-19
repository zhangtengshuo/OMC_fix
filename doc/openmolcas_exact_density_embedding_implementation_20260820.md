# Pure OpenMolcas exact-density M2 embedding implementation

Date: 2026-08-20  
Branch: `agent/exact-density-embedding`

## 1. Scope

This branch implements the electronic exact-density part of the fragment electrostatic embedding model used by Dimerge/NOCI-FR:

\[
\hat h_A^{\rm emb}
=
\hat h_A
+
\hat V_{\rm nuc}^{B}
+
\hat J[P_B].
\]

For target fragment A and source/environment fragment B,

\[
J^{B\to A}_{\mu\nu}
=
\sum_{\kappa\lambda\in B}
(\mu_A\nu_A|\kappa_B\lambda_B)
P^B_{\lambda\kappa}.
\]

The new `&EXACTEMB` module constructs `J[P_B]` entirely inside OpenMolcas from the final RASSCF AO density and the full-dimer Cholesky vectors. No PySCF object, AO bridge, multipole fit, OFE functional, exchange operator, or Pauli projector is used in the production path.

The partner nuclear attraction is kept on the already validated standard OpenMolcas route: target-fragment SEWARD includes the partner nuclei through `XField`. Therefore the full M2 preparation operator is

\[
V_{\rm emb}=V_{\rm nuc}^{\rm partner}+J[P_{\rm partner}].
\]

`EXACTEMB` itself stores only the electronic Coulomb matrix. This separation makes double counting explicit: partner nuclei belong in XField exactly once, and the partner electronic density belongs in `Exact Emb Pot` exactly once.

This embedding is for fragment-state preparation only. The final NOCI Hamiltonian must remain the bare physical full-dimer Hamiltonian; do not carry XField or `Exact Emb Pot` into the final NOCI B3 integral generation.

---

## 2. Implemented source changes

### 2.1 New OpenMolcas program

New directory:

```text
src/exactemb/
```

Files:

```text
CMakeLists.txt
main.F90
exactemb_data.F90
rdinp_exactemb.F90
exactemb.F90
```

OpenMolcas' top-level CMake logic auto-discovers lowercase program directories under `src`, so `src/exactemb/CMakeLists.txt` produces `exactemb.exe` through the standard `prog_template.cmake` path.

### 2.2 RASSCF one-electron Hamiltonian hook

Modified file:

```text
src/rasscf/sgfcin.F90
```

`SGFCIN` reads the normal `OneHam` and then checks the target RunFile for:

```text
Exact Emb Active
Exact Emb Pot
```

If active, it performs

\[
H^{(1)}_{\rm RASSCF}
\leftarrow
H^{(1)}_{\rm RASSCF}+J[P_{\rm partner}].
\]

The addition is made to `Tmp1`, the same one-electron Hamiltonian used to construct the core energy, inactive Fock matrix, active-space one-electron Hamiltonian, orbital gradient, and CI Hamiltonian. It is therefore not an energy-only correction.

`SGFCIN` explicitly aborts if stock `OFEMbedding` and this exact-density M2 hook are active at the same time.

---

## 3. Source density: final state-specific full RASSCF `D1ao`

No new density-export code is required.

At the end of RASSCF, the existing `PutRlx`/`Export1` path reconstructs and writes the selected relaxation-root density to the RunFile as:

```text
D1ao
```

The density contains frozen/inactive plus active contributions. `EXACTEMB` also reads:

```text
nFro
nIsh
nActel
Relax CASSCF root
```

and computes the expected source electron number

\[
N_e^{\rm expected}
=
2(n_{\rm Fro}+n_{\rm Ish})+n_{\rm ActEl}.
\]

The source RunFile must therefore be the **final RunFile from the accepted RASSCF state**, with `RLXROOT`/state tracking chosen consistently with the intended physical state. Do not use an intermediate RunFile from an unfinished root-tracking calculation.

---

## 4. Packed-density convention and why direct `D1ao` mapping is correct

OpenMolcas `Fold` converts a symmetric square density to lower-triangular expectation-value storage by keeping diagonal elements unchanged and multiplying off-diagonal elements by two:

\[
D^{\rm fold}_{ii}=D_{ii},
\qquad
D^{\rm fold}_{ij}=2D_{ij}\quad(i>j).
\]

Consequently, for a lower-triangular one-electron matrix `V`,

\[
D^{\rm fold}\cdot V^{\rm LT}
=
\operatorname{Tr}(PV).
\]

The existing Cholesky Fock code uses this packed density convention. `EXACTEMB` therefore maps the source RunFile's `D1ao` directly into the corresponding full-dimer lower-triangular AO block. No additional factor of two is applied.

---

## 5. AO mapping contract

Version 1 deliberately imposes a strict, transparent mapping contract:

1. full dimer, source fragment, and target fragment must all use `Group=NoSymm` / C1;
2. the source fragment AO basis must be identical to one contiguous AO block of the full dimer;
3. the target fragment AO basis must be identical to another contiguous AO block of the full dimer;
4. fragment geometry, atom order, basis assignments, spherical/Cartesian convention, and global orientation must match the corresponding full-dimer fragment;
5. source and target blocks must not overlap.

The input specifies only the first AO of each fragment in the dimer. The fragment AO counts are read from their RunFiles automatically.

Example for a dimer ordered A then B, each with 120 AOs:

```text
A: dimer AO 1...120
B: dimer AO 121...240
```

For A <- B:

```text
SOURCE first AO = 121
TARGET first AO = 1
```

For B <- A, swap them.

The code performs range and non-overlap checks and then a numerical electron-count invariant using the full-dimer overlap matrix:

\[
N_e^{\rm mapped}=D_{AB}^{\rm fold}\cdot S_{AB}^{\rm LT}.
\]

It hard-aborts unless

\[
|N_e^{\rm mapped}-N_e^{\rm expected}|<\texttt{TOLERANCE}.
\]

This catches many AO-order/basis mismatches before a Coulomb matrix is generated. It does not replace the requirement that the fragment and dimer geometries/basis conventions actually be identical.

---

## 6. Full-dimer Cholesky contraction

`&EXACTEMB` must be called while the **current OpenMolcas RunFile/ONEINT/Cholesky context belongs to the full AB dimer**.

The source folded density is embedded into a zero full-dimer density:

\[
D^{B,AB}_{\kappa\lambda}
=
\begin{cases}
D^B_{\kappa\lambda}, & \kappa,\lambda\in B,\\
0, & \text{otherwise}.
\end{cases}
\]

The implementation reuses the existing OpenMolcas `CHORAS_DRV` path with

```text
ALGO = 1
Deco = false
REORD = false
ExFac = 0
```

`ExFac=0` makes `CHORAS_DRV` Coulomb-only. The existing reduced-vector kernel then evaluates the Cholesky contraction equivalent to

\[
d_Q^B
=
\sum_{\kappa\lambda\in B}
L^Q_{\kappa\lambda}P^B_{\lambda\kappa},
\]

\[
J_{\mu\nu}^B
=
\sum_Q L^Q_{\mu\nu}d_Q^B.
\]

After the full-dimer `J` is built, `EXACTEMB` extracts only the target-target block and writes it to the target RunFile.

No exchange matrix is constructed. The source spin density is not used; open-shell T1/D+/D- states enter through their spin-summed `D1ao`.

---

## 7. RunFile records

`EXACTEMB` writes the following target RunFile fields:

| Record | Meaning |
|---|---|
| `Exact Emb Active` | nonzero => RASSCF must add the potential |
| `Exact Emb Ver` | format/version, currently 1 |
| `Exact Emb Pot` | packed target AO Coulomb matrix `J[P_partner]` |
| `Exact Emb Root` | source `Relax CASSCF root`, or -1 if absent |
| `Exact Emb SrcAO` | source first AO in full dimer |
| `Exact Emb TgtAO` | target first AO in full dimer |
| `Exact Emb DimAO` | full-dimer AO dimension |
| `Exact Emb Nelec` | mapped source `Tr[P S]` |
| `Exact Emb JFrob` | full symmetric Frobenius norm of target J |
| `Exact Emb JMax` | maximum absolute target J element |
| `Exact Emb Energy` | current RASSCF expectation value of the electronic embedding operator, written by SGFCIN |

The source density itself is not copied into the target RunFile.

---

## 8. `&EXACTEMB` input

### 8.1 Construct B -> A electronic Coulomb embedding

```text
&EXACTEMB
SOURCE
B.RunFile
121
TARGET
A.RunFile
1
TOLERANCE
1.0d-6
PRINT
1
END
```

Semantics:

```text
SOURCE
<source fragment RunFile>
<first source AO in current full-dimer AO order>

TARGET
<target fragment RunFile>
<first target AO in current full-dimer AO order>
```

`SOURCE` and `TARGET` file names are passed directly to OpenMolcas `NameRun`; they must be accessible in the current work directory/environment.

`PRINT` currently controls the summary diagnostics; nonnegative values print the construction summary.

### 8.2 Clear a target RunFile

```text
&EXACTEMB
CLEAR
A.RunFile
END
```

This sets `Exact Emb Active=0`. A stale potential can therefore remain physically present in the RunFile without being consumed by RASSCF.

---

## 9. Required file/context separation

There are two integral contexts and they must not be confused.

### 9.1 Full-dimer EXACTEMB context

Used only to build the Coulomb operator:

```text
AB.RunFile
AB.ONEINT
AB Cholesky/RI-CD files
B.RunFile   (source density)
A.RunFile   (target metadata / output record)
```

The current RunFile when `&EXACTEMB` starts must be the full-dimer AB RunFile, and the current Cholesky files must be those generated for the same dimer geometry/basis.

### 9.2 Target-fragment RASSCF context

Used to optimize A in B's electrostatic environment:

```text
A.RunFile containing Exact Emb Pot
A fragment ONEINT generated by SEWARD
partner-B nuclei represented by standard XField
A starting RasOrb/FileOrb
```

Do **not** run the target fragment RASSCF against the full-dimer ONEINT. Copy the modified target RunFile back to the target-fragment job/work directory, where its own target-fragment `ONEINT` contains the nuclear XField.

---

## 10. One-shot M2 workflow

For A <- B:

1. Run final isolated/source RASSCF for B at the exact dimer geometry/orientation and save `B0.RunFile`.
2. Prepare the target A integral context with standard SEWARD and partner-B nuclear XField.
3. Prepare a full-dimer AB C1 SEWARD Cholesky/RI-CD context with the same basis and AO order.
4. In the full-dimer context, make `B0.RunFile` and the target `A.RunFile` accessible.
5. Run `&EXACTEMB` with B as SOURCE and A as TARGET.
6. Return the modified `A.RunFile` to A's target-fragment work directory.
7. Run A RASSCF. SGFCIN reads `Exact Emb Pot` and optimizes with

\[
h_A+V_{\rm nuc}^B+J[P_B^{(0)}].
\]

8. Preserve the final `A1.RunFile`; its `D1ao` is the polarized source density for a subsequent freeze-thaw cycle.

Repeat symmetrically for B <- A.

---

## 11. Mutual freeze-thaw workflow

Define

\[
P_A^{(n+1)}=\mathcal F_A[P_B^{(n)}],
\qquad
P_B^{(n+1)}=\mathcal F_B[P_A^{(n)}].
\]

Each outer iteration is:

```text
B(n).RunFile --EXACTEMB in AB Cholesky context--> Exact Emb Pot in A target RunFile
A target SEWARD/XField + RASSCF ----------------> A(n+1).RunFile

A(n).RunFile --EXACTEMB in AB Cholesky context--> Exact Emb Pot in B target RunFile
B target SEWARD/XField + RASSCF ----------------> B(n+1).RunFile
```

The partner density is frozen for the entire inner RASSCF optimization. `EXACTEMB` is **not** called from a RASSCF microiteration. This preserves a clean outer freeze-thaw fixed-point problem.

For the first production implementation, keep the AB Cholesky context cached for a fixed geometry and regenerate only `Exact Emb Pot` from the updated source `D1ao`.

---

## 12. State-specific use

The embedding is state-conditioned. For example, the potential used for A:S1 in a B:S0 environment is generated from the final B:S0 density, whereas A:D+ in a B:D- environment uses the B:D- density.

The source RunFile must therefore be associated with the correct physical state and root. Recommended production metadata outside OpenMolcas should retain at least:

```text
fragment role
physical state label
charge
spin/multiplicity
accepted RASSCF root
geometry id
basis id
F&T cycle
```

`EXACTEMB` records the source `Relax CASSCF root`, but it cannot determine whether that root has the desired S1/T1/D+/D- character. Existing root tracking remains mandatory.

---

## 13. Meaning of "exact-density"

The environment is represented by its complete AO 1-RDM, not by fitted point charges, LoProp multipoles, or truncated transition-density descriptors. In this sense the density representation is exact within the chosen fragment wavefunction/basis.

The two-electron integrals are evaluated through the existing OpenMolcas Cholesky/RI-CD representation. Therefore the most precise name is:

> **Cholesky-accelerated exact-density electrostatic embedding**

or

> **OpenMolcas Cholesky exact-density M2**.

The residual integral error is controlled by the Cholesky/RI-CD approximation and its threshold. A small direct-ERI/PySCF calculation should remain the numerical reference for quantifying this error.

---

## 14. Hard exclusions in version 1

This implementation intentionally does **not** add:

- partner exchange `K[P_partner]`;
- Pauli projector / Huzinaga term;
- OFE nonadditive kinetic potential;
- OFE exchange-correlation potential;
- dispersion;
- intermolecular CT mixing during fragment orbital optimization;
- analytic embedding gradients;
- automatic fragment AO matching across arbitrary atom reorderings or rotations.

These are distinct physical/model extensions and should not be silently mixed into the M2 validation benchmark.

---

## 15. Compile/runtime constraints to check first

This branch was intentionally not compiled or numerically tested before handoff. The first build should therefore check the following interfaces first:

1. `exactemb.exe` is discovered by CMake and linked against `libmolcas`.
2. `Fock_util_interface::CHORAS_DRV` is visible to the new program.
3. the current OpenMolcas compiler accepts the `DSBA_Type` allocation/reference calls used in `exactemb.F90`.
4. `NameRun` can open the supplied source/target RunFile names in the chosen work directory.
5. source/target RunFiles were created with a parallel environment compatible with the `EXACTEMB` run; stock OpenMolcas RunFile opening can reject a changed process count.
6. the full-dimer Cholesky files remain available when `&EXACTEMB` executes.

No numerical conclusion should be drawn before these basic runtime contracts pass.

---

## 16. Recommended numerical validation order

### Test 1: no active record regression

Run an ordinary RASSCF without `Exact Emb Active`. Energy, orbitals, CI and output must be unchanged from `main`.

### Test 2: `CLEAR`

Create/retain an `Exact Emb Pot`, run `CLEAR`, and verify target RASSCF follows the ordinary Hamiltonian.

### Test 3: source electron-count invariant

For a known ethylene source state, verify

\[
\operatorname{Tr}(PS)=N_e
\]

to near machine precision. Deliberately use a wrong `SOURCE` AO start and verify that `EXACTEMB` aborts.

### Test 4: Coulomb matrix against PySCF direct J

For the same geometry, basis, AO convention, and source density, compare the OpenMolcas target block against the PySCF reference semantic

```python
get_jk((molA, molA, molB, molB), dmB,
       scripts='ijkl,lk->ij', aosym='s4')
```

The comparison isolates the Cholesky/RI-CD approximation and any AO-mapping errors.

### Test 5: J-only RASSCF response

Use an electronic potential with no partner nuclear XField and verify the RASSCF one-electron energy/orbital response is consistent with the stored `Exact Emb Energy` diagnostic.

### Test 6: full one-shot M2

Enable partner nuclear XField plus `Exact Emb Pot` and compare against the existing validated PySCF exact-density/OpenMolcas hybrid ethylene result.

Compare at minimum:

```text
CASSCF energy
D_P density distance
frontier / singly occupied NO overlap
state/root identity
```

### Test 7: mutual freeze-thaw

Reproduce the previous ethylene one-shot -> F&T density changes and convergence pattern.

### Test 8: real pentacene Stage F

Only after the ethylene matrix-level comparison passes, run the real pentacene comparison:

```text
isolated
LoProp multipoles
exact frozen-density electrostatics
mutual exact-density electrostatic polarization
```

and propagate the resulting fragment states into the unchanged bare-Hamiltonian NOCI calculation.

---

## 17. Failure diagnostics already implemented

`EXACTEMB` hard-aborts on:

- non-C1 dimer/source/target;
- non-Cholesky full-dimer context;
- missing or wrong-size source `D1ao`;
- source/target AO blocks outside the dimer AO range;
- overlapping source and target AO blocks;
- electron-count mismatch beyond `TOLERANCE`.

RASSCF hard-aborts on:

- active `Exact Emb Pot` with wrong packed dimension;
- exact-density embedding combined with OFEMbedding;
- non-C1 use of exact-density version 1.

Construction diagnostics report:

```text
source relaxation root
full-dimer/source/target AO ranges
expected electron count
Tr[P S]
|Delta N|
packed J norm
full symmetric Frobenius J norm
max |J_mn|
```

RASSCF reports once per execution that exact-density electronic embedding is active, together with source root/electron count, J norm/max element, and the current `D . J` expectation value.

---

## 18. Production rule for NOCI-FR

The resulting embedded RASSCF wavefunctions/orbitals may be used as fragment basis states in CommonOrb/NOCI-FR.

However:

\[
H_{IJ}^{\rm final}
=
\langle\Phi_I^{\rm emb}|H_{AB}|\Phi_J^{\rm emb}\rangle
\]

must still use the physical bare dimer Hamiltonian. Do not generate the final NOCI B3 integrals with partner XField or the `Exact Emb Pot` added to `OneHam`.

This is the boundary that prevents Hamiltonian-level double counting while retaining environment-induced orbital relaxation in the fragment states.
