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
!   J^B_{mn} = sum_kl (mn|kl) P^B_lk,
!
! by reusing the full-dimer OpenMolcas Cholesky vectors.  The source
! density is the final state-specific full AO/SO density D1ao exported
! by RASSCF.  The resulting target-fragment packed AO matrix is written
! to the target RunFile under the record "Exact Emb Pot".  RASSCF's
! SGFCIN adds this record to OneHam when "Exact Emb Active" is nonzero.
!
! The partner nuclear attraction is intentionally NOT generated here.
! It is supplied by the standard OpenMolcas XField route in the target
! SEWARD calculation.  Thus the complete M2 preparation Hamiltonian is
!
!   h_emb = h_fragment + V_nuc(partner) [XField] + J[P_partner].
!
! No exchange, non-additive kinetic, or XC embedding term is included.
!***********************************************************************

subroutine ExactEmb(iReturn)

use ExactEmb_Data, only: SourceRunFile, TargetRunFile, SourceFirstAO, TargetFirstAO, PrintLevel, ElectronTolerance, ClearMode
use Fock_util_interface, only: CHOras_drv
use Fock_util_global, only: ALGO, Deco, REORD
use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
use Index_Functions, only: iTri
use OneDat, only: sNoNuc, sNoOri
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: iReturn

integer(kind=iwp) :: i, j, ii, jj, iRc, iRcFinal, iComp, iOpt, iSyLbl
integer(kind=iwp) :: nSymD, nSymS, nSymT, nD, nS, nT, nTriD, nTriS, nTriT, nSqD
integer(kind=iwp) :: nData, nActElS, SourceRoot, ExactEmbVersion
integer(kind=iwp) :: nBasD(8), nBasS(8), nBasT(8), nFroS(8), nIshS(8), nOcc(8)
integer(kind=iwp) :: OldAlgo
integer(kind=iwp) :: SrcLast, TgtLast
real(kind=wp) :: ExpectedElectrons, MeasuredElectrons, DiffElectrons
real(kind=wp) :: JFrob, JMax, JPackedNorm
real(kind=wp) :: Dummy(1)
logical(kind=iwp) :: Found, DoCholesky, OldDeco, OldReord
character(len=8) :: OneLbl
real(kind=wp), allocatable :: DSource(:), DDimer(:), DSquare(:), JFull(:), JTarget(:), Overlap(:), CMODummy(:)
type(DSBA_Type) :: FSQ(1)
real(kind=wp), external :: dDot_

call RdInp_ExactEmb()
iReturn = 0
Dummy(1) = Zero
ExactEmbVersion = 1

if (ClearMode) then
  call NameRun(trim(TargetRunFile))
  call Put_iScalar('Exact Emb Active',0_iwp)
  call Put_iScalar('Exact Emb Ver',ExactEmbVersion)
  call Put_dArray('Exact Emb Pot',Dummy,0_iwp)
  call Put_dArray('Exact Emb Energy',Dummy,0_iwp)
  call NameRun('#Pop')
  write(u6,*)
  write(u6,*) 'EXACTEMB: embedding flag cleared in target RunFile: ',trim(TargetRunFile)
  return
end if

! ----------------------------------------------------------------------
! Full-dimer metadata.  EXACTEMB must be called while the current
! RunFile/Cholesky context belongs to the full AB dimer.
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
nSqD = nD*nD

call DecideOnCholesky(DoCholesky)
if (.not. DoCholesky) then
  write(u6,*) 'EXACTEMB: the current full-dimer integral context is not Cholesky/RI-CD.'
  write(u6,*) 'EXACTEMB: run full-dimer SEWARD with Cholesky/RICD before EXACTEMB.'
  call Abend()
end if

! ----------------------------------------------------------------------
! Read the final state-specific full density from the source fragment.
! PutRlx + Export1 guarantee that RASSCF D1ao contains frozen/inactive
! plus active density for the selected relaxation root.
! ----------------------------------------------------------------------
nBasS(:) = 0
nFroS(:) = 0
nIshS(:) = 0
SourceRoot = -1

call NameRun(trim(SourceRunFile))
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
  write(u6,*) 'EXACTEMB: source RunFile does not contain a compatible final D1ao.'
  write(u6,*) 'EXACTEMB: expected packed size = ',nTriS,' found = ',nData
  write(u6,*) 'EXACTEMB: run RASSCF to completion and preserve its final RunFile.'
  call Abend()
end if
call mma_allocate(DSource,nTriS,Label='ExactEmb_DSource')
call Get_dArray('D1ao',DSource,nTriS)

call Get_iArray('nFro',nFroS,nSymS)
call Get_iArray('nIsh',nIshS,nSymS)
call Get_iScalar('nActel',nActElS)
ExpectedElectrons = real(2*(nFroS(1)+nIshS(1))+nActElS,kind=wp)

call Qpg_iScalar('Relax CASSCF root',Found)
if (Found) call Get_iScalar('Relax CASSCF root',SourceRoot)
call NameRun('#Pop')

! ----------------------------------------------------------------------
! Read target fragment dimensions.  The target RunFile receives J.
! ----------------------------------------------------------------------
nBasT(:) = 0
call NameRun(trim(TargetRunFile))
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
! Validate explicit contiguous AO-block mapping into the full dimer.
! The fragment basis/order must be identical to the corresponding dimer
! AO block.  Electron-number validation below is an independent numerical
! guard against a wrong source block.
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
  write(u6,*) 'EXACTEMB: SOURCE and TARGET AO blocks overlap; this is not a fragment embedding map.'
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
! Electron-count invariant in the full-dimer AO representation.
! D1ao is already in OpenMolcas folded (packed expectation-value)
! convention, hence dot(D_packed,S_packed) = Tr[P S].
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
  write(u6,*) 'Check fragment/dimer AO ordering, basis, geometry, SOURCE first AO, and final RASSCF D1ao.'
  call Abend()
end if

! ----------------------------------------------------------------------
! Build J[P_source] in the full dimer AO space using the existing
! OpenMolcas Coulomb-only Cholesky Fock machinery.  CHORAS_DRV with
! ExFac=0 disables exchange.  For Coulomb, CHO_FOCKTWO_RED consumes only
! the packed DLT density; DSquare/CMODummy are deliberately zero dummies.
! ----------------------------------------------------------------------
call mma_allocate(DSquare,nSqD,Label='ExactEmb_DSquare')
call mma_allocate(CMODummy,nSqD,Label='ExactEmb_CMODummy')
call mma_allocate(JFull,nTriD,Label='ExactEmb_JFull')
DSquare(:) = Zero
CMODummy(:) = Zero
JFull(:) = Zero
nOcc(:) = 0
nOcc(1) = 1

call Allocate_DT(FSQ(1),nBasD,nBasD,nSymD)
FSQ(1)%A0(:) = Zero

OldAlgo = ALGO
OldDeco = Deco
OldReord = REORD
ALGO = 1
Deco = .false.
REORD = .false.

iRc = 0
call CHO_X_INIT(iRc,Zero)
if (iRc /= 0) then
  write(u6,*) 'EXACTEMB: CHO_X_INIT returned rc = ',iRc
  call Abend()
end if

call CHOras_drv(nSymD,nBasD,nOcc,DSquare,DDimer,JFull,Zero,FSQ,CMODummy)

iRcFinal = 0
call CHO_X_FINAL(iRcFinal)
if (iRcFinal /= 0) then
  write(u6,*) 'EXACTEMB: warning: CHO_X_FINAL returned rc = ',iRcFinal
end if

ALGO = OldAlgo
Deco = OldDeco
REORD = OldReord
call Deallocate_DT(FSQ(1))

! ----------------------------------------------------------------------
! Extract target-target block J_AA[P_B] and store it in target RunFile.
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

call NameRun(trim(TargetRunFile))
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
  write(u6,*) ' Source RunFile       : ',trim(SourceRunFile)
  write(u6,*) ' Target RunFile       : ',trim(TargetRunFile)
  write(u6,*) ' Source relaxation root: ',SourceRoot
  write(u6,*) ' Dimer AO count       : ',nD
  write(u6,*) ' Source AO block      : ',SourceFirstAO,SrcLast
  write(u6,*) ' Target AO block      : ',TargetFirstAO,TgtLast
  write(u6,'(A,F20.10)') ' Expected electrons   : ',ExpectedElectrons
  write(u6,'(A,F20.10)') ' Tr[P S]              : ',MeasuredElectrons
  write(u6,'(A,ES14.6)') ' |Delta N|            : ',DiffElectrons
  write(u6,'(A,ES14.6)') ' packed ||J||_2       : ',JPackedNorm
  write(u6,'(A,ES14.6)') ' full-matrix ||J||_F  : ',JFrob
  write(u6,'(A,ES14.6)') ' max |J_mn|           : ',JMax
  write(u6,*) ' Exchange contribution: DISABLED (ExFac = 0)'
  write(u6,*) ' Partner nuclei       : supply with standard SEWARD XField'
  write(u6,*) ' Target RunFile record: Exact Emb Pot'
  write(u6,*) '============================================================'
  write(u6,*)
end if

call mma_deallocate(JTarget)
call mma_deallocate(JFull)
call mma_deallocate(CMODummy)
call mma_deallocate(DSquare)
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
