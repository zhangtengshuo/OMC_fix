!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
!***********************************************************************

module ExactEmb_Data

use Definitions, only: wp, iwp

implicit none
private

character(len=8) :: SourceRunName = ''
character(len=8) :: TargetRunName = ''
integer(kind=iwp) :: SourceFirstAO = 0_iwp
integer(kind=iwp) :: TargetFirstAO = 0_iwp
integer(kind=iwp) :: PrintLevel = 1_iwp
real(kind=wp) :: ElectronTolerance = 1.0e-6_wp
logical(kind=iwp) :: ClearMode = .false.

public :: SourceRunName, TargetRunName, SourceFirstAO, TargetFirstAO, PrintLevel, ElectronTolerance, ClearMode
public :: Reset_ExactEmb_Input

contains

subroutine Reset_ExactEmb_Input()
  SourceRunName = ''
  TargetRunName = ''
  SourceFirstAO = 0_iwp
  TargetFirstAO = 0_iwp
  PrintLevel = 1_iwp
  ElectronTolerance = 1.0e-6_wp
  ClearMode = .false.
end subroutine Reset_ExactEmb_Input

end module ExactEmb_Data
