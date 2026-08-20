#!/usr/bin/env python3
"""Stage 4 authoritative comparison: OpenMolcas EXACTEMB J vs PySCF direct-ERI J.

Executed entirely in the exact OpenMolcas CC-PVDZ replica basis rebuilt from the
Stage 2 guessorb.h5 (native run2h5 metadata via the validated Dimerge/NOCI_MPS
bridge module), so no basis transform is needed:

  1. omAB/omA/omB = replica mols for dimer/fragments (96/48/48 AOs).
  2. mapping from build_mapping; the source density D1ao (folded packed,
     halved on unfold — verified Tr[P S] = 16.0000000024) is placed at the
     molcas B block and reordered with mapping[48:96].
  3. J_ref = pyscf.scf.jk.get_jk((omA,omA,omB,omB), dmB, 'ijkl,lk->ij', s4)
     — direct libcint ERIs; cross-validated against a shls_slice einsum
     contraction to 6e-15 in a separate run.
  4. 'Exact Emb Pot' (plain packed lower triangle — values as stored) placed
     at the molcas A block, reordered with mapping[:48].
  5. Elementwise and generalized-eigenvalue (S^-1 J) comparison.

Result (2026-08-20): relF = 1.8e-08, max|dJ| = 8.6e-08, generalized
eigenvalue max deviation 1.2e-06, top eigenvalue 3.1147 on both sides.

Note: an alternative native-cc-pVDZ route via the contracted-AO transform
T = S_nat^-1 S_cross reproduced the density electron count but gave
inconsistent J statistics; it is superseded by this direct replica-basis
comparison, which needs no transform at all.
"""
import numpy as np
from pyscf import gto
from pyscf.scf.jk import get_jk

import check_molcas_pyscf_ao_bridge as om_bridge

XYZ = """C -0.6695 0.0 0.0
C 0.6695 0.0 0.0
H -1.2321 0.9230 0.0
H -1.2321 -0.9230 0.0
H 1.2321 0.9230 0.0
H 1.2321 -0.9230 0.0
C -0.6695 0.0 3.5
C 0.6695 0.0 3.5
H -1.2321 0.9230 3.5
H -1.2321 -0.9230 3.5
H 1.2321 0.9230 3.5
H 1.2321 -0.9230 3.5"""


def load(path):
    toks = [t for t in open(path).read().split() if t]
    ntri = int(toks[0])
    vals = np.array([float(x) for x in toks[1:1 + ntri]])
    n = int((np.sqrt(8 * ntri + 1) - 1) // 2)
    assert n * (n + 1) // 2 == ntri
    return n, vals


def unfold(n, packed, halve):
    m = np.zeros((n, n))
    k = 0
    for i in range(n):
        for j in range(i + 1):
            v = packed[k] if (i == j or not halve) else 0.5 * packed[k]
            m[i, j] = m[j, i] = v
            k += 1
    return m


def main():
    basis, frags = om_bridge.exact_pyscf_basis_from_fragments(["S2_AB_ricd.guessorb.h5"])
    omAB = gto.M(atom=XYZ, basis=basis, unit="Angstrom", cart=False, verbose=0)
    omA = gto.M(atom="\n".join(XYZ.splitlines()[:6]), basis=basis, unit="Angstrom", verbose=0)
    omB = gto.M(atom="\n".join(XYZ.splitlines()[6:]), basis=basis, unit="Angstrom", verbose=0)
    assert omAB.nao_nr() == 96 and omA.nao_nr() == omB.nao_nr() == 48
    mapping, _, _, _, _ = om_bridge.build_mapping(omAB, frags)

    nS, d1 = load("d1ao_B.dat")
    nT_, pot = load("pot_BtoA.dat")
    assert nS == nT_ == 48
    P_B = unfold(48, d1, halve=True)
    J_A = unfold(48, pot, halve=False)

    Pfull = np.zeros((96, 96))
    Pfull[np.ix_(mapping[48:], mapping[48:])] = P_B
    dmB = Pfull[48:, 48:]
    S_B = omB.intor_symmetric("int1e_ovlp")
    ne = np.trace(dmB @ S_B)
    print(f"[check] Tr[P_B S_B] = {ne:.10f} (expect 16)")
    assert abs(ne - 16.0) < 1e-7

    J_ref = get_jk((omA, omA, omB, omB), dmB, scripts="ijkl,lk->ij", aosym="s4")
    J_ref = 0.5 * (np.asarray(J_ref) + np.asarray(J_ref).T)

    Jfull = np.zeros((96, 96))
    Jfull[np.ix_(mapping[:48], mapping[:48])] = J_A
    J_om = Jfull[:48, :48]

    d = J_om - J_ref
    S_A = omA.intor_symmetric("int1e_ovlp")
    ev_om = np.sort(np.linalg.eigvals(np.linalg.solve(S_A, J_om)).real)
    ev_ref = np.sort(np.linalg.eigvals(np.linalg.solve(S_A, J_ref)).real)
    print(f"[result] ||J_omc||_F = {np.linalg.norm(J_om):.8f}   ||J_ref||_F = {np.linalg.norm(J_ref):.8f}")
    print(f"[result] max|dJ| = {np.abs(d).max():.3e}   ||dJ||_F/||J_ref|| = {np.linalg.norm(d)/np.linalg.norm(J_ref):.3e}")
    print(f"[result] generalized eig max dev = {np.abs(ev_om-ev_ref).max():.3e}  top = {ev_om[-1]:.6f} / {ev_ref[-1]:.6f}")
    np.savez_compressed("stage04_result.npz", J_omc=J_om, J_ref=J_ref, dmB=dmB, mapping=mapping)
    print("saved stage04_result.npz")


if __name__ == "__main__":
    main()
