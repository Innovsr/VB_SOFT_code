!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fort_py_connect(flgst, nao, nae, mult)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
use commondat1
implicit none

integer::wig1,nnae,STDOUT,argnum,nstr7,str5(2000,20)
!integer::nflgst, nnao, nnae, nmult
logical :: fileexists
character(len=35)::inputfilename
common/str/str5,nstr7

noq0=100
noq1=100
noq2=100
noq3=100
symm=1
ovval=1.000
ovopt=0
nlpset=0
nfset=0
flg1=1
mult=1
flgst=1
key_frag=0

print*,flgst, nao, nae, mult


open(unit=7,file='structures.dat',status='unknown')
open(unit=10,file='structure_set_1.dat',status='unknown')
open(unit=23,file='Rumer_Sets_all.dat',status='unknown')
open(unit=35,file='quality_str.dat',status='unknown')

serial=0
input_flg=0
nstr7=0
Rumwrite=0

!argnum=iargc()
!if(argnum.eq.0)then
!inputfilename='input.dat'
!input_flg=1
!INQUIRE(FILE=TRIM(inputfilename),EXIST=fileexists)
!IF (fileexists) THEN
!print*,'Reading from input.dat file.'
!!call input
!ELSE
!PRINT*,'SORRY input.dat file does not exist and you also did not put .xmi file'
!stop
!ENDIF
!endif

if(input_flg.eq.0)call read_xmi


if(flgst.eq.1.or.flgst.eq.2) then
call cov_struc
endif
!if(flgst.eq.1.or.flgst.eq.3) then
!!call ion_struc
!endif

call close_file

return
end subroutine fort_py_connect
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
