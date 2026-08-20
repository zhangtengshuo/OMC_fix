#!/usr/bin/env python3
"""Stage 5 verification (final consolidated form).

Three gates, all in the exact OpenMolcas CC-PVDZ replica basis:

  G1  XField composition: pyscf CASSCF(h + V_nuc^B) + (A-nuc x B-charge)
      constant == OpenMolcas A0 nuclei-only RASSCF total (-78.82129253).
  G2  J-only (Test 5): pyscf CASSCF(h + J_ref) == OM -42.13929260.
  G3  Full M2 (Test 6): pyscf CASSCF(h + V + J_ref, started from XField-SCF
      orbitals exactly like OM's FileOrb=SCFORB start) + constant
      == OM -41.66496895.

The originally observed 26 mEh M2 gap was CASSCF start/landscape sensitivity:
starting from the embedded-hcore SCF lands on a different local solution;
starting from the XField SCF orbitals (the OM start) reproduces OM to 2.8e-07.

Also documented: the SGFCIN <J[P_partner]> printout is emitted once at the
first SGFCIN call with the inactive-only (14e) density; it is not the final
density expectation (final Tr[D J_ref] = 35.43/36.47 vs printed 31.70/33.38).
"""
import numpy as np
from pyscf import gto
from types import MethodType

import check_molcas_pyscf_ao_bridge as om_bridge

XYZ_A = """C -0.6695 0.0 0.0
C 0.6695 0.0 0.0
H -1.2321 0.9230 0.0
H -1.2321 -0.9230 0.0
H 1.2321 0.9230 0.0
H 1.2321 -0.9230 0.0"""
XYZ_B = """C -0.6695 0.0 3.5
C 0.6695 0.0 3.5
H -1.2321 0.9230 3.5
H -1.2321 -0.9230 3.5
H 1.2321 0.9230 3.5
H 1.2321 -0.9230 3.5"""

OM = {"A0_xfield": -78.82129253, "Jonly": -42.13929260, "M2": -41.66496895}


def main():
    basis, _ = om_bridge.exact_pyscf_basis_from_fragments(["S2_AB_ricd.guessorb.h5"])
    omA = gto.M(atom=XYZ_A, basis=basis, unit="Angstrom", verbose=0)
    omB = gto.M(atom=XYZ_B, basis=basis, unit="Angstrom", verbose=0)
    J = np.load("stage04_result.npz")["J_ref"]
    Enn = sum(omA.atom_charge(a) * omB.atom_charge(b) /
              np.linalg.norm(omA.atom_coord(a) - omB.atom_coord(b))
              for a in range(omA.natm) for b in range(omB.natm))
    V = np.zeros((48, 48))
    for b in range(omB.natm):
        with omA.with_rinv_origin(omB.atom_coord(b)):
            V -= omB.atom_charge(b) * omA.intor("int1e_rinv")
    V = 0.5 * (V + V.T)
    print(f"[const] E(A-nuc x B-charge) = {Enn:.8f} Eh")

    def hf_with(pot):
        mf = omA.RHF()
        hc = mf.get_hcore() + pot
        hc = 0.5 * (hc + hc.T)
        mf.get_hcore = MethodType(lambda self, mol=None: hc, mf)
        mf.conv_tol = 1e-12
        return mf

    def cas(mf, mo=None):
        mc = mf.CASSCF(2, 2)
        mc.conv_tol = 1e-11
        mc.conv_tol_grad = 1e-9
        mc.kernel(mo_coeff=mo)
        return mc.e_tot

    mf_v = hf_with(V); mf_v.kernel()
    e1 = cas(mf_v) + Enn
    print(f"[G1 XField ] {e1:.8f} vs OM {OM['A0_xfield']:.8f}  diff {e1 - OM['A0_xfield']:+.3e}")

    mf_j = hf_with(J); mf_j.kernel()
    e2 = cas(mf_j)
    print(f"[G2 J-only ] {e2:.8f} vs OM {OM['Jonly']:.8f}  diff {e2 - OM['Jonly']:+.3e}")

    mo_xf = mf_v.mo_coeff
    mf_m2 = hf_with(V + J); mf_m2.kernel(mo_coeff=mo_xf)
    e3 = cas(mf_m2, mo=mo_xf) + Enn
    print(f"[G3 M2     ] {e3:.8f} vs OM {OM['M2']:.8f}  diff {e3 - OM['M2']:+.3e}")


if __name__ == "__main__":
    main()
