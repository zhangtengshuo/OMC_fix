!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
!***********************************************************************
!
! EXACTEMB constructs the electronic part of the M2 exact-density
! electrostatic embedding operator,
!
!   J^B_mn = sum_kl (mn|kl) P^B_lk,
!
! using the full-dimer OpenMolcas Cholesky vectors. The source density
! is the final state-specific full AO/SO density D1ao exported by
! RASSCF. The target-fragment block is written to the target RunFile as
! "Exact Emb Pot" and is consumed by RASSCF/SGFCIN.
!
! The partner nuclear attraction is deliberately NOT generated here.
! It is supplied by the standard OpenMolcas XField route in the target
! SEWARD calculation. The complete M2 preparation Hamiltonian is thus
!
!   h_emb = h_fragment + V_nuc(partner) [XField] + J[P_partner].
!
! No exchange, non-additive kinetic, XC, Pauli-projector, or dispersion
! contribution is included by this module.
!***********************************************************************

subroutine ExactEmb(iReturn)

use ExactEmb_Data, only: SourceRunName, TargetRunName, SourceFirstAO, TargetFirstAO, PrintLevel, ElectronTolerance, ClearMode
use Index_Functions, only: iTri
use OneDat, only: sNoNuc, sNoOri
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: iReturn

integer(kind=iwp) :: i, j, ii, jj, iRc, iComp, iOpt, iSyLbl
integer(kind=iwp) :: nSymD, nSymS, nSymT, nD, nS, nT, nTriD, nTriS, nTriT
integer(kind=iwp) :: nData, nActElS, SourceRoot, ExactEmbVersion
integer(kind=iwp) :: nBasD(8), nBasS(8), nBasT(8), nFroS(8), nIshS(8)
integer(kind=iwp) :: SrcLast, TgtLast
real(kind=wp) :: ExpectedElectrons, MeasuredElectrons, DiffElectrons
real(kind=wp) :: JFrob, JMax, JPackedNorm
logical(kind=iwp) :: Found, DoCholesky
character(len=8) :: OneLbl
real(kind=wp), allocatable :: DSource(:), DDimer(:), JFull(:), JTarget(:), Overlap(:)
real(kind=wp), external :: dDot_

call RdInp_ExactEmb()
iReturn = 0
ExactEmbVersion = 1

! ----------------------------------------------------------------------
! CLEAR is intentionally only a switch. A stale matrix may remain on the
! RunFile but RASSCF will not consume it while Exact Emb Active == 0.
! ----------------------------------------------------------------------
if (ClearMode) then
  call NameRun(TargetRunName)
  call Put_iScalar('Exact Emb Active',0_iwp)
  call Put_iScalar('Exact Emb Ver',ExactEmbVersion)
  call Put_dScalar('Exact Emb Energy',Zero)
  call NameRun('#Pop')
  write(u6,*)
  write(u6,*) 'EXACTEMB: embedding flag cleared in target RunFile: ',trim(TargetRunName)
  return
end if

! ----------------------------------------------------------------------
! Full-dimer metadata. EXACTEMB must be called while the current
! RunFile/ONEINT/Cholesky context belongs to the full AB dimer.
! ----------------------------------------------------------------------
nBasD(:) = 0
call Get_iScalar('nSym',nSymD)
if (nSymD /= 1) then
  write(u6,*) 'EXACTEMB: version 1 requires full-dimer Group=NoSymm (nSym=1).'
  write(u6,*) 'EXACTEMB: dimer nSym = ',nSymD
  call Abend()
end if
call Get_iArray('nBas',nBasD,nSymD)
nD = nBasD(1)
nTriD = nD*(nD+1)/2

call DecideOnCholesky(DoCholesky)
if (.not. DoCholesky) then
  write(u6,*) 'EXACTEMB: current full-dimer integral context is not Cholesky/RI-CD.'
  write(u6,*) 'EXACTEMB: run full-dimer SEWARD with Cholesky/RICD before EXACTEMB.'
  call Abend()
end if

! ----------------------------------------------------------------------
! Read the final state-specific full density from the source fragment.
! Existing RASSCF PutRlx/Export1 writes the final relaxation-root D1ao.
! ----------------------------------------------------------------------
nBasS(:) = 0
nFroS(:) = 0
nIshS(:) = 0
SourceRoot = -1
nData = 0

call NameRun(SourceRunName)
call Get_iScalar('nSym',nSymS)
if (nSymS /= 1) then
  write(u6,*) 'EXACTEMB: source fragment must use Group=NoSymm (nSym=1).'
  call Abend()
end if
call Get_iArray('nBas',nBasS,nSymS)
nS = nBasS(1)
nTriS = nS*(nS+1)/2

call Qpg_dArray('D1ao',Found,nData)
if ((.not. Found) .or. (nData /= nTriS)) then
  write(u6,*) 'EXACTEMB: source RunFile lacks a compatible final D1ao.'
  write(u6,*) 'EXACTEMB: expected packed size = ',nTriS,' found = ',nData
  write(u6,*) 'EXACTEMB: run RASSCF to completion and preserve its final RunFile.'
  call Abend()
end if
call mma_allocate(DSource,nTriS,Label='ExactEmb_DSource')
call Get_dArray('D1ao',DSource,nTriS)

call Get_iArray('nFro',nFroS,nSymS)
call Get_iArray('nIsh',nIshS,nSymS)
call Qpg_iScalar('nActel',Found)
if (.not. Found) then
  write(u6,*) 'EXACTEMB: source RunFile lacks nActel; cannot validate full electron count.'
  call Abend()
end if
call Get_iScalar('nActel',nActElS)
ExpectedElectrons = real(2*(nFroS(1)+nIshS(1))+nActElS,kind=wp)

call Qpg_iScalar('Relax CASSCF root',Found)
if (Found) call Get_iScalar('Relax CASSCF root',SourceRoot)
call NameRun('#Pop')

! ----------------------------------------------------------------------
! Read target dimensions. The target RunFile receives the J operator.
! ----------------------------------------------------------------------
nBasT(:) = 0
call NameRun(TargetRunName)
call Get_iScalar('nSym',nSymT)
if (nSymT /= 1) then
  write(u6,*) 'EXACTEMB: target fragment must use Group=NoSymm (nSym=1).'
  call Abend()
end if
call Get_iArray('nBas',nBasT,nSymT)
nT = nBasT(1)
nTriT = nT*(nT+1)/2
call NameRun('#Pop')

! ----------------------------------------------------------------------
! Version 1 has an explicit contiguous AO-block mapping contract.
! Fragment geometry/orientation/basis/AO convention must match the
! corresponding full-dimer block exactly.
! ----------------------------------------------------------------------
SrcLast = SourceFirstAO+nS-1
TgtLast = TargetFirstAO+nT-1
if ((SourceFirstAO < 1) .or. (SrcLast > nD)) then
  write(u6,*) 'EXACTEMB: SOURCE AO block lies outside the dimer AO range.'
  write(u6,*) 'SOURCE = ',SourceFirstAO,SrcLast,' DIMER = 1,',nD
  call Abend()
end if
if ((TargetFirstAO < 1) .or. (TgtLast > nD)) then
  write(u6,*) 'EXACTEMB: TARGET AO block lies outside the dimer AO range.'
  write(u6,*) 'TARGET = ',TargetFirstAO,TgtLast,' DIMER = 1,',nD
  call Abend()
end if
if (.not. ((SrcLast < TargetFirstAO) .or. (TgtLast < SourceFirstAO))) then
  write(u6,*) 'EXACTEMB: SOURCE and TARGET AO blocks overlap.'
  call Abend()
end if

call mma_allocate(DDimer,nTriD,Label='ExactEmb_DDimer')
DDimer(:) = Zero

do i=1,nS
  ii = SourceFirstAO+i-1
  do j=1,i
    jj = SourceFirstAO+j-1
    DDimer(iTri(ii,jj)) = DSource(iTri(i,j))
  end do
end do

! ----------------------------------------------------------------------
! Electron-count invariant in full-dimer representation. D1ao is in
! OpenMolcas Fold convention, so D_fold . S_LT = Tr[P S].
! ----------------------------------------------------------------------
call mma_allocate(Overlap,nTriD+4,Label='ExactEmb_Overlap')
iRc = -1
iOpt = ibset(ibset(0,sNoOri),sNoNuc)
iComp = 1
iSyLbl = 1
OneLbl = 'Mltpl  0'
call RdOne(iRc,iOpt,OneLbl,iComp,Overlap,iSyLbl)
if (iRc /= 0) then
  write(u6,*) 'EXACTEMB: failed to read full-dimer overlap matrix from ONEINT.'
  call Abend()
end if
MeasuredElectrons = dDot_(nTriD,DDimer,1,Overlap,1)
DiffElectrons = abs(MeasuredElectrons-ExpectedElectrons)

if (DiffElectrons > ElectronTolerance) then
  write(u6,*)
  write(u6,*) 'EXACTEMB: electron-count invariant FAILED.'
  write(u6,'(A,F20.10)') '  expected electrons = ',ExpectedElectrons
  write(u6,'(A,F20.10)') '  Tr[P S]             = ',MeasuredElectrons
  write(u6,'(A,ES14.6)') '  absolute difference = ',DiffElectrons
  write(u6,*) 'Check geometry/orientation, basis/AO ordering, SOURCE first AO, and final RASSCF D1ao.'
  call Abend()
end if

! ----------------------------------------------------------------------
! Existing fock_util wrapper performs a Coulomb-only Cholesky build:
! ExFac=0, ALGO=1, Deco=false, REORD=false. No K[P] is evaluated.
! ----------------------------------------------------------------------
call mma_allocate(JFull,nTriD,Label='ExactEmb_JFull')
call Cho_ExactEmb_Coulomb(nSymD,nBasD,DDimer,JFull)

! ----------------------------------------------------------------------
! Extract J_target,target[P_source].
! ----------------------------------------------------------------------
call mma_allocate(JTarget,nTriT,Label='ExactEmb_JTarget')
JTarget(:) = Zero

do i=1,nT
  ii = TargetFirstAO+i-1
  do j=1,i
    jj = TargetFirstAO+j-1
    JTarget(iTri(i,j)) = JFull(iTri(ii,jj))
  end do
end do

call Packed_Matrix_Norms(JTarget,nT,JPackedNorm,JFrob,JMax)

call NameRun(TargetRunName)
call Put_iScalar('Exact Emb Active',1_iwp)
call Put_iScalar('Exact Emb Ver',ExactEmbVersion)
call Put_iScalar('Exact Emb Root',SourceRoot)
call Put_iScalar('Exact Emb SrcAO',SourceFirstAO)
call Put_iScalar('Exact Emb TgtAO',TargetFirstAO)
call Put_iScalar('Exact Emb DimAO',nD)
call Put_dScalar('Exact Emb Nelec',MeasuredElectrons)
call Put_dScalar('Exact Emb JFrob',JFrob)
call Put_dScalar('Exact Emb JMax',JMax)
call Put_dArray('Exact Emb Pot',JTarget,nTriT)
call NameRun('#Pop')

if (PrintLevel >= 0) then
  write(u6,*)
  write(u6,*) '============================================================'
  write(u6,*) ' OpenMolcas EXACTEMB: exact-density Coulomb embedding built'
  write(u6,*) '============================================================'
  write(u6,*) ' Source RunFile        : ',trim(SourceRunName)
  write(u6,*) ' Target RunFile        : ',trim(TargetRunName)
  write(u6,*) ' Source relaxation root: ',SourceRoot
  write(u6,*) ' Dimer AO count        : ',nD
  write(u6,*) ' Source AO block       : ',SourceFirstAO,SrcLast
  write(u6,*) ' Target AO block       : ',TargetFirstAO,TgtLast
  write(u6,'(A,F20.10)') ' Expected electrons    : ',ExpectedElectrons
  write(u6,'(A,F20.10)') ' Tr[P S]               : ',MeasuredElectrons
  write(u6,'(A,ES14.6)') ' |Delta N|             : ',DiffElectrons
  write(u6,'(A,ES14.6)') ' packed ||J||_2        : ',JPackedNorm
  write(u6,'(A,ES14.6)') ' full-matrix ||J||_F   : ',JFrob
  write(u6,'(A,ES14.6)') ' max |J_mn|            : ',JMax
  write(u6,*) ' Exchange contribution : DISABLED'
  write(u6,*) ' Partner nuclei        : supply with standard SEWARD XField'
  write(u6,*) ' Target RunFile record : Exact Emb Pot'
  write(u6,*) '============================================================'
  write(u6,*)
end if

! PRINT >= 2 is intended only for small-system matrix-level regression
! against a direct-ERI/PySCF J reference. Do not enable for pentacene.
if (PrintLevel >= 2) then
  call TriPrt('EXACTEMB target J[P_source]','(5G20.12)',JTarget,nT)
end if

call mma_deallocate(JTarget)
call mma_deallocate(JFull)
call mma_deallocate(Overlap)
call mma_deallocate(DDimer)
call mma_deallocate(DSource)

contains

subroutine Packed_Matrix_Norms(AP,n,PNorm,FNorm,AMax)
  real(kind=wp), intent(in) :: AP(*)
  integer(kind=iwp), intent(in) :: n
  real(kind=wp), intent(out) :: PNorm, FNorm, AMax
  integer(kind=iwp) :: ia, ja, ij
  real(kind=wp) :: psq, fsq_local

  psq = Zero
  fsq_local = Zero
  AMax = Zero
  do ia=1,n
    do ja=1,ia
      ij = iTri(ia,ja)
      psq = psq+AP(ij)*AP(ij)
      if (ia == ja) then
        fsq_local = fsq_local+AP(ij)*AP(ij)
      else
        fsq_local = fsq_local+2.0_wp*AP(ij)*AP(ij)
      end if
      AMax = max(AMax,abs(AP(ij)))
    end do
  end do
  PNorm = sqrt(psq)
  FNorm = sqrt(fsq_local)
end subroutine Packed_Matrix_Norms

end subroutine ExactEmb
