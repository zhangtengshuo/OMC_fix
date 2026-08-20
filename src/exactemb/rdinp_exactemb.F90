!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
!***********************************************************************

subroutine RdInp_ExactEmb()

use ExactEmb_Data, only: SourceRunName, TargetRunName, SourceFirstAO, TargetFirstAO, PrintLevel, ElectronTolerance, ClearMode, &
                         Reset_ExactEmb_Input
use spool, only: SpoolInp
use Definitions, only: wp, iwp, u6

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
    call ReadRunNameLine(LuSpool,SourceRunName,'SOURCE RunFile logical name')
    call ReadIntLine(LuSpool,SourceFirstAO,'SOURCE first AO in dimer')
  case ('TARG')
    call ReadRunNameLine(LuSpool,TargetRunName,'TARGET RunFile logical name')
    call ReadIntLine(LuSpool,TargetFirstAO,'TARGET first AO in dimer')
  case ('TOLE')
    call ReadRealLine(LuSpool,ElectronTolerance,'TOLERANCE')
  case ('PRIN')
    call ReadIntLine(LuSpool,PrintLevel,'PRINT')
  case ('CLEA')
    ClearMode = .true.
    call ReadRunNameLine(LuSpool,TargetRunName,'CLEAR target RunFile logical name')
  case ('END ')
    Done = .true.
  case default
    write(u6,*) 'EXACTEMB: unknown input keyword: ',trim(Line)
    call Abend()
  end select
end do

if (len_trim(TargetRunName) == 0) then
  write(u6,*) 'EXACTEMB: TARGET RunFile logical name was not specified.'
  call Abend()
end if

if (ClearMode) return

if (len_trim(SourceRunName) == 0) then
  write(u6,*) 'EXACTEMB: SOURCE RunFile logical name was not specified.'
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

subroutine ReadRunNameLine(LU,Value,What)
  integer(kind=iwp), intent(in) :: LU
  character(len=8), intent(out) :: Value
  character(len=*), intent(in) :: What
  character(len=256) :: LocalLine

  call ReadTextLine(LU,LocalLine,What)
  LocalLine = trim(adjustl(LocalLine))
  if (len_trim(LocalLine) > len(Value)) then
    write(u6,*) 'EXACTEMB: ',trim(What),' must contain at most 8 characters: ',trim(LocalLine)
    call Abend()
  end if
  if ((index(LocalLine,'/') > 0) .or. (index(LocalLine,achar(92)) > 0)) then
    write(u6,*) 'EXACTEMB: ',trim(What),' must be a short name in the OpenMolcas work directory, not a path: ',trim(LocalLine)
    call Abend()
  end if
  Value = trim(LocalLine)
end subroutine ReadRunNameLine

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
