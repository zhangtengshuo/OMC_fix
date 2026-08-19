# EXACTEMB example workflow

These templates exercise the pure OpenMolcas exact-density electronic embedding implementation on branch `agent/exact-density-embedding`.

## 1. Two distinct OpenMolcas contexts are required

### Full-dimer AB context: construct `J[P_partner]`

The `&EXACTEMB` call must run with the current OpenMolcas files belonging to the full dimer:

```text
AB.RunFile
AB.ONEINT
AB Cholesky/RI-CD vector files
```

The full dimer must be C1 (`Group=NoSymm`) in version 1.

The final source-fragment RASSCF RunFile and target-fragment RunFile must also be accessible in this work directory. `&EXACTEMB` reads `D1ao` from the source RunFile, contracts it with the AB Cholesky vectors, extracts the target-target AO block, and stores that matrix as `Exact Emb Pot` in the target RunFile.

### Target-fragment context: consume the operator in RASSCF

After `&EXACTEMB`, return/copy the modified target RunFile to the target fragment's own work directory. Its SEWARD/ONEINT must describe only the target fragment plus the partner **nuclei** through the standard OpenMolcas XField route.

RASSCF then sees

\[
h^{emb}=h_{fragment}+V_{nuc}^{partner}+J[P_{partner}].
\]

Do not run the target RASSCF against the full-dimer ONEINT.

## 2. AO-order contract

Version 1 assumes source and target are non-overlapping contiguous AO blocks of the full dimer. The fragment coordinates/orientation, atom order, basis assignment and AO convention must be identical to the corresponding dimer blocks.

For example, if

```text
A = dimer AO   1...120
B = dimer AO 121...240
```

then:

- `exactemb_A_from_B.input`: B -> A uses SOURCE first AO 121, TARGET first AO 1.
- `exactemb_B_from_A.input`: A -> B uses SOURCE first AO 1, TARGET first AO 121.

The program validates the mapped source density by requiring

\[
|\operatorname{Tr}(PS)-N_e| < \texttt{TOLERANCE}.
\]

A wrong AO start should therefore normally fail before the Cholesky contraction.

## 3. Small-system matrix regression

For ethylene/debugging, set:

```text
PRINT
2
```

`EXACTEMB` then prints the complete target lower-triangular `J[P_source]` matrix. Compare it element-by-element with the PySCF direct-ERI reference semantic:

```python
get_jk((molA, molA, molB, molB), dmB,
       scripts='ijkl,lk->ij', aosym='s4')
```

Do not use `PRINT 2` for pentacene-sized production calculations.

## 4. One-shot sequence

For one simultaneous/Jacobi-style update from isolated densities:

```text
A0.RunFile, B0.RunFile
       |
       +-- EXACTEMB B0 -> A target RunFile -- target A RASSCF --> A1.RunFile
       |
       +-- EXACTEMB A0 -> B target RunFile -- target B RASSCF --> B1.RunFile
```

The partner density is frozen during each complete RASSCF optimization.

## 5. Mutual freeze-thaw sequence

Repeat at the outer level:

```text
B(n) -- EXACTEMB --> A potential -- RASSCF --> A(n+1)
A(n) -- EXACTEMB --> B potential -- RASSCF --> B(n+1)
```

For a strict Jacobi iteration, build both new potentials from cycle `n` densities before accepting `A(n+1), B(n+1)`. A sequential/Gauss-Seidel variant can instead use the newly accepted A density immediately for B, but do not mix the two conventions in a benchmark.

`EXACTEMB` is never called inside a RASSCF microiteration.

## 6. Disable a stored potential

Use:

```text
exactemb_clear.input
```

This sets `Exact Emb Active=0`. The old `Exact Emb Pot` record may remain in the RunFile but RASSCF will not consume it.

## 7. Required first validation order

1. ordinary RASSCF on this branch with no `Exact Emb Active` record: must reproduce `main`;
2. `CLEAR` regression;
3. `Tr[P S]` electron-count check;
4. ethylene `PRINT 2` OpenMolcas J vs PySCF direct-J matrix;
5. J-only RASSCF response;
6. full `XField + Exact Emb Pot` one-shot M2 vs the existing ethylene reference;
7. mutual freeze-thaw;
8. only then real pentacene Stage F.

The branch was intentionally handed off without compilation or numerical execution; the user-side build is the first compiler/runtime validation.

## 8. NOCI boundary

The embedded fragment orbitals/wavefunctions may feed CommonOrb/NOCI-FR, but the final dimer NOCI integrals must be generated with the bare physical dimer Hamiltonian:

\[
H_{IJ}=\langle\Phi_I^{emb}|H_{AB}|\Phi_J^{emb}\rangle.
\]

Do not leave either fragment XField or `Exact Emb Pot` in the final B3/NOCI Hamiltonian.
