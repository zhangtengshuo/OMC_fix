!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
!***********************************************************************
!
! Coulomb-only wrapper used by the standalone EXACTEMB program.
!
! The wrapper deliberately lives inside fock_util so that the program
! does not depend directly on the internal Fock_util_interface and
! Fock_util_global module files.  DLT uses the standard OpenMolcas Fold
! convention (diagonal once, off-diagonal twice).  FLT is returned in
! lower-triangular storage and contains only J[DLT].
!***********************************************************************

subroutine Cho_ExactEmb_Coulomb(nSym,nBas,W_DLT,W_FLT)

use Fock_util_interface, only: CHOras_drv
use Fock_util_global, only: ALGO, Deco, REORD
use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(8)
real(kind=wp), intent(in) :: W_DLT(*)
real(kind=wp), intent(out) :: W_FLT(*)

integer(kind=iwp) :: iRc, iRcFinal, nSq, OldAlgo, nOcc(8)
logical(kind=iwp) :: OldDeco, OldReord
type(DSBA_Type) :: FSQ(1)
real(kind=wp), allocatable :: W_DSQ(:), W_CMO(:)

if (nSym /= 1) then
  write(u6,*) 'Cho_ExactEmb_Coulomb: version 1 requires nSym=1.'
  call Abend()
end if

nSq = sum(nBas(1:nSym)*nBas(1:nSym))
call mma_allocate(W_DSQ,nSq,Label='ExactEmb_DSQ')
call mma_allocate(W_CMO,nSq,Label='ExactEmb_CMO')
W_DSQ(:) = Zero
W_CMO(:) = Zero
W_FLT(1:nBas(1)*(nBas(1)+1)/2) = Zero

nOcc(:) = 0
! CHORAS_DRV uses nOcc to decide whether a symmetry block can be skipped.
! For Coulomb-only contraction the numerical density is entirely W_DLT;
! a positive sentinel keeps the single C1 block active.
nOcc(1) = 1

call Allocate_DT(FSQ(1),nBas,nBas,nSym)
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
  write(u6,*) 'Cho_ExactEmb_Coulomb: CHO_X_INIT returned rc = ',iRc
  call Abend()
end if

! ExFac=0 => DoExchange=.false.; FactC remains exactly one.
call CHOras_drv(nSym,nBas,nOcc,W_DSQ,W_DLT,W_FLT,Zero,FSQ,W_CMO)

iRcFinal = 0
call CHO_X_FINAL(iRcFinal)
if (iRcFinal /= 0) then
  write(u6,*) 'Cho_ExactEmb_Coulomb: warning: CHO_X_FINAL returned rc = ',iRcFinal
end if

ALGO = OldAlgo
Deco = OldDeco
REORD = OldReord

call Deallocate_DT(FSQ(1))
call mma_deallocate(W_CMO)
call mma_deallocate(W_DSQ)

end subroutine Cho_ExactEmb_Coulomb
