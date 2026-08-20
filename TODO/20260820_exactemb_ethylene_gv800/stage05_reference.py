#!/usr/bin/env python3
"""Stage 5 reference and verification.

1. PySCF CASSCF(2e,2o) on the replica molA with embedded hcores
   (hcore + J_ref[P_B]) and (hcore + V_nuc^B + J_ref[P_B]),
   compared against the OpenMolcas embedded RASSCF root-1 energies
   (-42.13929260 and -41.66496895).
2. External recomputation of the SGFCIN <J[P_partner]> diagnostic from the
   preserved final densities (Tr[D_final . J_ref] vs the printed values
   3.169901E+01 and 3.338251E+01).
"""
import numpy as np
from pyscf import gto, scf
from types import MethodType

import check_molcas_pyscf_ao_bridge as om_bridge

XYZ_A = """C -0.6695 0.0 0.0
C 0.6695 0.0 0.0
H -1.2321 0.9230 0.0
H -1.2321 -0.9230 0.0
H 1.2321 0.9230 0.0
H 1.2321 -0.9230 0.0"""
XYZ_B = XYZ_A.replace("0.0\n", "3.5\n").replace(" 0.0", "")  # placeholder, built explicitly below
XYZ_B = """C -0.6695 0.0 3.5
C 0.6695 0.0 3.5
H -1.2321 0.9230 3.5
H -1.2321 -0.9230 3.5
H 1.2321 0.9230 3.5
H 1.2321 -0.9230 3.5"""
XYZ_AB = XYZ_A + "\n" + XYZ_B


def load(path):
    toks = [t for t in open(path).read().split() if t]
    ntri = int(toks[0])
    vals = np.array([float(x) for x in toks[1:1 + ntri]])
    return vals


def unfold_fold(n, packed):
    m = np.zeros((n, n))
    k = 0
    for i in range(n):
        for j in range(i + 1):
            v = packed[k] if i == j else 0.5 * packed[k]
            m[i, j] = m[j, i] = v
            k += 1
    return m


def casscf_in_potential(mol, potential, label, e_om):
    mf = mol.RHF()
    hcore = np.asarray(mf.get_hcore(), dtype=float) + potential
    hcore = 0.5 * (hcore + hcore.T)

    def get_hcore(this, mol=None):
        return hcore
    mf.get_hcore = MethodType(get_hcore, mf)
    mf.conv_tol = 1e-12
    mf.kernel()
    mc = mf.CASCI(2, 2)
    mc.kernel()
    ncas_root = mc.e_tot
    mcs = mf.CASSCF(2, 2)
    mcs.conv_tol = 1e-10
    mcs.conv_tol_grad = 1e-8
    mcs.kernel()
    print(f"[{label}] pyscf CASCI={ncas_root:.8f} CASSCF={mcs.e_tot:.8f} "
          f"OM RASSCF={e_om:.8f} d(CASSCF)={mcs.e_tot - e_om:+.3e}")
    return mcs.e_tot


def main():
    basis, frags = om_bridge.exact_pyscf_basis_from_fragments(["S2_AB_ricd.guessorb.h5"])
    omAB = gto.M(atom=XYZ_AB, basis=basis, unit="Angstrom", cart=False, verbose=0)
    omA = gto.M(atom=XYZ_A, basis=basis, unit="Angstrom", verbose=0)
    omB = gto.M(atom=XYZ_B, basis=basis, unit="Angstrom", verbose=0)
    mapping, _, _, _, _ = om_bridge.build_mapping(omAB, frags)
    assert np.array_equal(np.sort(mapping[:48]), np.arange(48))
    assert np.array_equal(np.sort(mapping[48:]), np.arange(48, 96))

    z = np.load("stage04_result.npz")
    J_ref = z["J_ref"]  # replica basis, A block 48x48

    # partner nuclear attraction on A from B nuclei (point charges, Bohr coords)
    V = np.zeros((48, 48))
    for atom in range(omB.natm):
        q = omB.atom_charge(atom)
        with omA.with_rinv_origin(omB.atom_coord(atom)):
            V -= q * omA.intor("int1e_rinv")
    V = 0.5 * (V + V.T)

    casscf_in_potential(omA, J_ref, "J-only", -42.13929260)
    casscf_in_potential(omA, V + J_ref, "M2    ", -41.66496895)

    # external <J> diagnostic check from preserved final densities
    for label, path, printed in (("J-only", "d1ao_A1_Jonly.dat", 31.69901),
                                 ("M2", "d1ao_A1_M2.dat", 33.38251)):
        P48 = unfold_fold(48, load(path))
        Pfull = np.zeros((96, 96))
        Pfull[np.ix_(mapping[:48], mapping[:48])] = P48
        P_rep = Pfull[:48, :48]
        dj = np.trace(P_rep @ J_ref)
        print(f"[<J> {label}] external Tr[D J_ref] = {dj:.6f}   OM printed = {printed:.6f}   "
              f"diff = {dj - printed:+.2e}")


if __name__ == "__main__":
    main()
