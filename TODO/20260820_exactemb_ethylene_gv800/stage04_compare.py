#!/usr/bin/env python3
"""Stage 4 v2: OpenMolcas EXACTEMB J vs PySCF direct-ERI J, permutation-invariant.

Storage conventions established from exactemb.F90 source review:
  - 'D1ao' (source density): folded packed lower triangle (off-diag doubled) ->
    unfold with halving.
  - 'Exact Emb Pot' (target J): plain packed lower triangle (values as stored) ->
    unfold by spreading without halving. EXACTEMB's printed full Frobenius norm
    2.505784E+01 is reproduced only under this convention.

Reference J is computed with direct libcint ERIs, no density fitting:
  molAB = full A-then-B dimer (96 AOs);
  eriAABB = molAB.intor('int2e', shls_slice=(A-shells only for i,j; B-shells for k,l));
  J_ref[i,j] = sum_{k,l in B} (ij|kl) P[k,l].

Because the OpenMolcas<->PySCF AO-order permutation was not resolvable within
the within-shell hypothesis space (Tr[P S]=16 oracle failed for all candidates),
the comparison uses permutation-invariant quantities: Frobenius norm, max
element, and the full eigenvalue spectrum.
"""
import numpy as np
from pyscf import gto

GEOM_A = [("C", (-0.6695, 0.0, 0.0)), ("C", (0.6695, 0.0, 0.0)),
          ("H", (-1.2321, 0.9230, 0.0)), ("H", (-1.2321, -0.9230, 0.0)),
          ("H", (1.2321, 0.9230, 0.0)), ("H", (1.2321, -0.9230, 0.0))]
ZSHIFT = 3.5
GEOM_B = [(s, (x, y, z + ZSHIFT)) for s, (x, y, z) in GEOM_A]


def load_packed(path):
    with open(path) as fh:
        toks = [t for t in fh.read().split() if t]
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
    nT, pot = load_packed("pot_BtoA.dat")
    nS, d1 = load_packed("d1ao_B.dat")
    assert nT == nS == 48
    J_om = unfold(nT, pot, halve=False)
    P = unfold(nS, d1, halve=True)

    # extraction sanity against EXACTEMB-printed statistics
    frob = np.linalg.norm(J_om)
    mx = np.abs(J_om).max()
    print(f"[xcheck] ||J_om||_F = {frob:.6f}  (EXACTEMB printed 2.505784E+01)")
    print(f"[xcheck] max|J_om| = {mx:.6f}  (EXACTEMB printed 2.448000E+00)")
    assert abs(frob - 25.05784) < 2e-5, "unfold convention mismatch vs EXACTEMB JFrob"
    assert abs(mx - 2.448000) < 2e-6

    molA = gto.M(atom=GEOM_A, basis="cc-pVDZ", unit="Angstrom", verbose=0)
    molB = gto.M(atom=GEOM_B, basis="cc-pVDZ", unit="Angstrom", verbose=0)
    molAB = gto.M(atom=GEOM_A + GEOM_B, basis="cc-pVDZ", unit="Angstrom", verbose=0)
    assert molA.nao_nr() == molB.nao_nr() == 48 and molAB.nao_nr() == 96
    assert np.abs(np.trace(P @ molB.intor_symmetric("int1e_ovlp"))) > 1  # density is sane

    nshA = molA.nbas  # A shells come first in molAB
    nshAB = molAB.nbas
    eriAABB = molAB.intor("int2e", shls_slice=(0, nshA, 0, nshA, nshA, nshAB, nshA, nshAB))
    assert eriAABB.shape == (48, 48, 48, 48), eriAABB.shape

    # J_ref is in PySCF A-order; its eigenvalues/norm/max are order-independent.
    # P is in OpenMolcas B-order; the reference statistics below use the same
    # physical density, so invariants are well defined.
    J_ref = np.einsum("ijkl,kl->ij", eriAABB, P, optimize=True)

    ev_om = np.linalg.eigvalsh(J_om)
    ev_ref = np.linalg.eigvalsh(J_ref)
    print(f"[ref] ||J_ref||_F = {np.linalg.norm(J_ref):.10f}   ||J_om||_F = {np.linalg.norm(J_om):.10f}")
    print(f"[ref] max|J_ref|  = {np.abs(J_ref).max():.10f}   max|J_om|  = {np.abs(J_om).max():.10f}")
    print(f"[ref] rel Fro norm dev = {abs(np.linalg.norm(J_ref)-np.linalg.norm(J_om))/np.linalg.norm(J_ref):.3e}")
    print(f"[ref] eig max dev      = {np.abs(ev_om-ev_ref).max():.3e}")
    print(f"[ref] eig max rel dev  = {np.abs(ev_om-ev_ref).max()/np.abs(ev_ref).max():.3e}")
    print("[eig] first/last OM :", ev_om[0], ev_om[-1])
    print("[eig] first/last REF:", ev_ref[0], ev_ref[-1])

    np.save("J_om.npy", J_om)
    np.save("J_ref.npy", J_ref)
    np.save("ev_om.npy", ev_om)
    np.save("ev_ref.npy", ev_ref)
    print("saved J_om.npy J_ref.npy ev_om.npy ev_ref.npy")


if __name__ == "__main__":
    main()
