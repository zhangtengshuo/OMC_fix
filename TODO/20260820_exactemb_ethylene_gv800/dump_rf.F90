program Dump_RF
use Definitions, only: wp, iwp
use OneDat, only: sNoNuc, sNoOri
implicit none
character(len=8) :: RFName
character(len=64) :: Label, OutFile
integer(kind=iwp) :: nSym, nBas(8), nFro(8), nIsh(8), nActel, nData, iu, i, iRc, iOpt
logical :: Found
real(kind=wp), allocatable :: Arr(:)

call getarg(1, RFName)
call getarg(2, Label)
call getarg(3, OutFile)
if (len_trim(RFName) == 0 .or. len_trim(Label) == 0 .or. len_trim(OutFile) == 0) then
  write(*,*) 'usage: dump_rf <RunFileShortName> <Label> <OutFile>'
  stop 2
end if

call IniMem()

write(*,*) 'RunFile = ', trim(RFName)
call NameRun(RFName)
call Get_iScalar('nSym', nSym)
write(*,*) 'nSym = ', nSym
call Get_iArray('nBas', nBas, nSym)
write(*,*) 'nBas = ', nBas(1:nSym)
call Get_iArray('nFro', nFro, nSym)
call Get_iArray('nIsh', nIsh, nSym)
write(*,*) 'nFro/nIsh = ', nFro(1), nIsh(1)
call Qpg_iScalar('nActel', Found)
if (Found) then
  call Get_iScalar('nActel', nActel)
  write(*,*) 'nActel = ', nActel
end if

call Qpg_dArray(trim(Label), Found, nData)
write(*,*) 'Qpg_dArray ', trim(Label), ' Found=', Found, ' nData=', nData
if (.not. Found .and. trim(Label) /= 'OVLPC') then
  write(*,*) 'record not found'
  stop 1
end if

if (trim(Label) == 'OVLPC') then
  ! read the packed overlap from the ONEINT attached to the current RunFile
  nData = nBas(1)*(nBas(1)+1)/2
  allocate(Arr(nData+4))
  call RdOne(iRc,ibset(ibset(0,sNoOri),sNoNuc),'Mltpl  0',1,Arr,1)
  write(*,*) 'RdOne overlap iRc=', iRc, ' nData=', nData
  if (iRc /= 0) then
    write(*,*) 'RdOne failed'
    stop 1
  end if
else
  allocate(Arr(nData))
  call Get_dArray(trim(Label), Arr, nData)
end if

open(newunit=iu, file=trim(OutFile), status='unknown', action='write')
write(iu,'(I0)') nData
do i = 1, nData
  write(iu,'(ES25.17)') Arr(i)
end do
close(iu)
write(*,*) 'DUMPED ', trim(Label), ' -> ', trim(OutFile), ' nData=', nData
stop 0
end program Dump_RF
