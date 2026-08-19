!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
!***********************************************************************

program Main

#ifdef _FPE_TRAP_
use, intrinsic :: IEEE_Exceptions, only: IEEE_Set_Halting_Mode, IEEE_Usual
use Definitions, only: DefInt
#endif
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: rc

#ifdef _FPE_TRAP_
call IEEE_Set_Halting_Mode(IEEE_Usual,.true._DefInt)
#endif

call Start('exactemb')
call ExactEmb(rc)
call Finish(rc)

end program Main
