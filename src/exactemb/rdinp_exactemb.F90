!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
!***********************************************************************

subroutine RdInp_ExactEmb()

use ExactEmb_Data, only: SourceRunFile, TargetRunFile, SourceFirstAO, TargetFirstAO, PrintLevel, ElectronTolerance, ClearMode, &
                         Reset_ExactEmb_Input
use spool, only: SpoolInp
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: LuSpool, ios
character(len=256) :: Line, KeyWord
logical(kind=iwp) :: Done

call Reset_ExactEmb_Input()

LuSpool = 17
call SpoolInp(LuSpool)
rewind(LuSpool)
call RdNLst(LuSpool,'EXACTEMB')

Done = .false.
do while (.not. Done)
  read(LuSpool,'(A)',iostat=ios) Line
  if (ios /= 0) then
    write(u6,*) 'EXACTEMB: unexpected end of input; END keyword is required.'
    call Abend()
  end if
  if ((len_trim(Line) == 0) .or. (Line(1:1) == '*')) cycle

  KeyWord = adjustl(Line)
  call UpCase(KeyWord)

  select case (KeyWord(1:min(4,len_trim(KeyWord))))
  case ('SOUR')
    call ReadTextLine(LuSpool,SourceRunFile,'SOURCE RunFile')
    call ReadIntLine(LuSpool,SourceFirstAO,'SOURCE first AO in dimer')
  case ('TARG')
    call ReadTextLine(LuSpool,TargetRunFile,'TARGET RunFile')
    call ReadIntLine(LuSpool,TargetFirstAO,'TARGET first AO in dimer')
  case ('TOLE')
    call ReadRealLine(LuSpool,ElectronTolerance,'TOLERANCE')
  case ('PRIN')
    call ReadIntLine(LuSpool,PrintLevel,'PRINT')
  case ('CLEA')
    ClearMode = .true.
    call ReadTextLine(LuSpool,TargetRunFile,'CLEAR target RunFile')
  case ('END ')
    Done = .true.
  case default
    write(u6,*) 'EXACTEMB: unknown input keyword: ',trim(Line)
    call Abend()
  end select
end do

if (len_trim(TargetRunFile) == 0) then
  write(u6,*) 'EXACTEMB: TARGET RunFile was not specified.'
  call Abend()
end if

if (ClearMode) return

if (len_trim(SourceRunFile) == 0) then
  write(u6,*) 'EXACTEMB: SOURCE RunFile was not specified.'
  call Abend()
end if
if (SourceFirstAO < 1) then
  write(u6,*) 'EXACTEMB: SOURCE first AO must be >= 1.'
  call Abend()
end if
if (TargetFirstAO < 1) then
  write(u6,*) 'EXACTEMB: TARGET first AO must be >= 1.'
  call Abend()
end if
if (ElectronTolerance <= 0.0_wp) then
  write(u6,*) 'EXACTEMB: TOLERANCE must be positive.'
  call Abend()
end if

contains

subroutine ReadTextLine(LU,Value,What)
  integer(kind=iwp), intent(in) :: LU
  character(len=*), intent(out) :: Value
  character(len=*), intent(in) :: What
  character(len=256) :: LocalLine
  integer(kind=iwp) :: ierr

  do
    read(LU,'(A)',iostat=ierr) LocalLine
    if (ierr /= 0) then
      write(u6,*) 'EXACTEMB: missing value for ',trim(What)
      call Abend()
    end if
    if ((len_trim(LocalLine) == 0) .or. (LocalLine(1:1) == '*')) cycle
    Value = trim(adjustl(LocalLine))
    exit
  end do
end subroutine ReadTextLine

subroutine ReadIntLine(LU,Value,What)
  integer(kind=iwp), intent(in) :: LU
  integer(kind=iwp), intent(out) :: Value
  character(len=*), intent(in) :: What
  character(len=256) :: LocalLine
  integer(kind=iwp) :: ierr

  call ReadTextLine(LU,LocalLine,What)
  read(LocalLine,*,iostat=ierr) Value
  if (ierr /= 0) then
    write(u6,*) 'EXACTEMB: invalid integer for ',trim(What),': ',trim(LocalLine)
    call Abend()
  end if
end subroutine ReadIntLine

subroutine ReadRealLine(LU,Value,What)
  integer(kind=iwp), intent(in) :: LU
  real(kind=wp), intent(out) :: Value
  character(len=*), intent(in) :: What
  character(len=256) :: LocalLine
  integer(kind=iwp) :: ierr

  call ReadTextLine(LU,LocalLine,What)
  read(LocalLine,*,iostat=ierr) Value
  if (ierr /= 0) then
    write(u6,*) 'EXACTEMB: invalid real value for ',trim(What),': ',trim(LocalLine)
    call Abend()
  end if
end subroutine ReadRealLine

end subroutine RdInp_ExactEmb
