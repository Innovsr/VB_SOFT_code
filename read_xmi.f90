!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_xmi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
use commondat1
implicit none

integer::argnum,l,MDP,sl(50),k,ll,j,i,i1,io,endflg(50),l1,l2,gap(10,100),keyn(15),least&
,l3,l4,k1,bfimat(5,20),gap1(20),dash(20),bfiln(5),&
lastlbfi,ufrzorbn,orbn1,orbn2,gap2(200),k2,k3,k4,k5,k6,k7,nfrag
integer::chinst,orbsl1(50),k8,k9,sl2(10),MDP1
logical :: fileexists
character(len=35)::inputfilename,keywd(15,50),line7(50),line
character(len=4)::lname
character(len=55)::sttr,charst(500),string(4,50)
character(len=90)::line3,line2,line4,line5,line6
character(len=200)::lowercase
integer::atsymset(20,20),nsym,syn(50),at_sym(50),fragorb,ifrag,iorb,ibfi

common/ats/atsymset,nsym,syn,at_sym
print*,'enter read_xmi'

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


!if(argnum.eq.1)then
!endif

if(input_flg.eq.1)goto 500

call getarg(1,inputfilename)
l=len(TRIM(inputfilename))
lname=inputfilename(l-3:l)

if(lname.ne.'.xmi') then
print*,'please put the .xmi file as input'
stop
endif
INQUIRE(FILE=TRIM(inputfilename),EXIST=fileexists)
IF (fileexists) THEN
open(unit=21,file=TRIM(inputfilename),status='old')
ELSE
PRINT*,'SORRY This input file does not exist or you may not provide the &
filename at all'
stop
ENDIF
write(9,*)inputfilename
write(23,*)inputfilename

MDP=0
do
read(21,'(a)',iostat=io)

if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(21)
!print*,'MDP',MDP

do i=1,MDP-1
read(21,'(a)')charst(i)
if(index(charst(i),'$end').ne.0)then
MDP1=i
endif
enddo

rewind(21)

do i=1,10
sl1(i)=0
enddo
l=0
l1=0

do i=1,MDP1

if(index(charst(i),'#').eq.1.or.index(charst(i),';').eq.1)goto 789
if(index(charst(i),'#').gt.1)then
ll=index(charst(i),'#')
do k=1,ll
line=charst(i)(k:k)
if(line.eq.'')then
if(trim(charst(i)(1:k-1)).ne.'')then
charst(i)=charst(i)(1:ll-1)
goto 111
endif
endif
enddo
endif

if(index(charst(i),';').gt.1)then
ll=index(charst(i),';')
do k=1,ll
line=charst(i)(k:k)
if(line.eq.'')then
if(trim(charst(i)(1:k-1)).ne.'')then
charst(i)=charst(i)(1:ll-1)
goto 111
endif
endif
enddo
endif


111 charst(i)=lowercase(charst(i))

if(index(charst(i),'$ctrl').ne.0.or.index(charst(i),'$end').ne.0&
.or.index(charst(i),'$orb').ne.0.or.index(charst(i),'$frag').ne.0&
.or.index(charst(i),'$bfi').ne.0) then
l=l+1
if(index(charst(i),'$end').ne.0)then
endflg(l)=1
else
endflg(l)=0
endif
if(endflg(l).ne.endflg(l-1).or.l.eq.1)then
l1=l1+1
sl(l1)=i
if(index(charst(i),'$ctrl').ne.0)sl1(l1)=1
if(index(charst(i),'$frag').ne.0)sl1(l1)=2
if(index(charst(i),'$orb').ne.0)sl1(l1)=3
if(index(charst(i),'$bfi').ne.0)sl1(l1)=4
endif
if(endflg(l).eq.0.and.endflg(l-1).eq.0)then
print*,'You might forget to put $end somewhere in the .xmi file'
stop
endif
endif
789 enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! if fragment section present then orbital section should be read before
!!!!!!!  fragment section

fragorb=0
ifrag=2
iorb=3
k1=0
k2=0
do i=1,10
sl2(i)=0
enddo
do i=1,l1
if(sl1(i).eq.2)then
k1=1
k3=i
endif
enddo
do i=1,l1
if(sl1(i).eq.3)then
k2=1
k4=i
endif
enddo
fragorb=k1*k2
if(fragorb.eq.1)then
sl2(k3)=sl(k4)
sl2(k3+1)=sl(k4+1)
sl2(k4)=sl(k3)
sl2(k4+1)=sl(k3+1)
ifrag=3
iorb=2
do i=1,l1
if(sl2(i).ne.0)then
sl(i)=sl2(i)
endif
enddo
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
l=0
l4=0
do i=1,l1,2
l4=l4+1
l3=0
do j=sl(i)+1,sl(i+1)-1
if(sl1(i).eq.1)then
l=l+1
l2=1
gap(1,1)=0
ll=len(trim(charst(j)))
do k=1,ll+1
line5=charst(j)(k:k)
if(line5.eq.'#')goto 110
if(line5.eq.'')then
l2=l2+1
gap(l,l2)=k
endif
enddo
110 do k=1,l2
if((gap(l,k+1)-gap(l,k)).ne.1)then
l3=l3+1
keywd(l4,l3)=trim(charst(j)(gap(l,k):gap(l,k+1)-1))
endif
enddo
endif
enddo

keyn(l4)=l3
do k=1,l3
if(trim(keywd(l4,k)(2:7)).eq.'frgtyp')key_frag=1
if(trim(keywd(l4,k)(2:4)).eq.'nae')then
ll=len(trim(keywd(l4,k)))
do k1=1,ll+1
line5=keywd(l4,k)(k1:k1)
if(line5.eq.'=')then
line3=keywd(l4,k)(k1+1:ll)
read(line3,'(I2)')nae
endif
enddo
endif
enddo



do k=1,l3
if(trim(keywd(l4,k)(2:4)).eq.'nao')then
ll=len(trim(keywd(l4,k)))
do k1=1,ll+1
line5=keywd(l4,k)(k1:k1)
if(line5.eq.'=')then
line3=keywd(l4,k)(k1+1:ll)
read(line3,'(I2)')nao
endif
enddo
endif
enddo


do k=1,l3
if(trim(keywd(l4,k)(2:5)).eq.'nmul')then
ll=len(trim(keywd(l4,k)))
do k1=1,ll+1
line5=keywd(l4,k)(k1:k1)
if(line5.eq.'=')then
line3=keywd(l4,k)(k1+1:ll)
read(line3,'(I2)')mult
endif
enddo
endif
enddo

do k=1,l3
if(trim(keywd(l4,k)(2:7)).eq.'chinst')then
ll=len(trim(keywd(l4,k)))
do k1=1,ll+1
line5=keywd(l4,k)(k1:k1)
if(line5.eq.'=')then
line3=keywd(l4,k)(k1+1:ll)
read(line3,'(I2)')chinst
endif
enddo
if(chinst.eq.1)then
flg1=0
itb=1
syb=3
nnb=2
radical=0
mnbond=0
endif
endif
enddo


do k=1,l3
if(trim(keywd(l4,k)(2:4)).eq.'str')then
ll=len(trim(keywd(l4,k)))
do k1=1,ll+1
line5=keywd(l4,k)(k1:k1)
if(line5.eq.'=')then
line3=keywd(l4,k)(k1+1:ll)
if(trim(line3).eq.'full')flgst=1
if(trim(line3).eq.'cov')flgst=2
if(trim(line3).eq.'ion')flgst=3

endif
enddo
endif
enddo
rewind(21)

if(sl1(i).eq.iorb)then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! reading of orbital section !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


ll=0
l=0
do j=sl(i)+1,sl(i+1)-1
ll=len(trim(charst(j)))
line3=trim(charst(j)(1:1))
if(ll.eq.0)goto 300
if(trim(line3).eq.'#')goto 300
l=l+1
if(l.eq.1)goto 300
k2=0
do k=1,ll+1
line3=trim(charst(j)(k:k))
if(line3.eq.';')goto 301
if(line3.eq.'')then
k2=k2+1
gap2(k2)=k
endif
enddo
301 if(gap2(1).eq.1)then
k5=0
do k=2,k2
if(abs(gap2(k)-gap2(k-1)).ne.1)then
line3=trim(charst(j)(gap2(k-1)+1:gap2(k)-1))
ll=len(trim(line3))
do k3=1,ll
line2=trim(line3(k3:k3))
if(line2.eq.'-')then
line5=trim(charst(j)(gap2(k-1)+1:gap2(k-1)+k3-1))
line6=trim(charst(j)(gap2(k-1)+k3+1:gap2(k)-1))
read(line5,'(I10)')orbn1
read(line6,'(I10)')orbn2
do k4=orbn1,orbn2
k5=k5+1
orbs(l,k5)=k4
enddo
goto 444
endif
enddo
k5=k5+1
read(line3,'(I10)')orbs(l,k5)
endif
444 enddo
endif
if(gap2(k).ne.1)then
k5=0
do k=1,k2
if(k.eq.1)then
line3=trim(charst(j)(1:gap2(k)-1))
ll=len(trim(line3))
do k3=1,ll
line2=trim(line3(k3:k3))
if(line2.eq.'-')then
line5=trim(charst(j)(1:k3-1))
line6=trim(charst(j)(k3+1:gap2(k)-1))
read(line5,'(I10)')orbn1
read(line6,'(I10)')orbn2
do k4=orbn1,orbn2
k5=k5+1
orbs(l,k5)=k4
enddo
goto 445
endif
 enddo
k5=k5+1
read(line3,'(I10)')orbs(l,k5)
endif
if(k.ne.1)then
if(abs(gap2(k)-gap2(k-1)).ne.1)then
line3=trim(charst(j)(gap2(k-1)+1:gap2(k)-1))
ll=len(trim(line3))
do k3=1,ll
line2=trim(line3(k3:k3))
if(line2.eq.'-')then
line5=trim(charst(j)(gap2(k-1)+1:gap2(k-1)+k3-1))
line6=trim(charst(j)(gap2(k-1)+k3+1:gap2(k)-1))
read(line5,'(I10)')orbn1
read(line6,'(I10)')orbn2
do k4=orbn1,orbn2
k5=k5+1
orbs(l,k5)=k4
enddo
goto 445
endif
 enddo
k5=k5+1
read(line3,'(I10)')orbs(l,k5)
endif
endif
445 enddo
endif

orbn(l)=k5
300 enddo
rewind(21)
orbsl=l
niao=l-1-nao
endif

do l=1,orbsl
print*,'orbitals',orbsl,niao,orbn(l),(orbs(l,k),k=1,4)
enddo

if(sl1(i).eq.4)then

ll=0
l3=0
do j=sl(i)+1,sl(i+1)-1
l3=l3+1
gap1(1)=0
ll=len(trim(charst(j)))
l=1
l1=0
do k=1,ll+1
line3=trim(charst(j)(k:k))
if(line3.eq.''.or.line3.eq.'-')then
l=l+1
gap1(l)=k
if(line3.eq.'-')dash(l)=l
endif
enddo
do k=2,l
if(gap1(k)-gap1(k-1).gt.1)then
l1=l1+1
line3=trim(charst(j)(gap1(k-1)+1:gap1(k)-1))
read(line3,'(I3)')bfimat(l3,l1)
if(dash(k).eq.k)then
l1=l1+1
bfimat(l3,l1)=0
endif
endif
enddo
bfiln(l3)=l1
enddo
lastlbfi=l3
l2=0
j=lastlbfi
do k=1,bfiln(j)-1
if(k.eq.1.and.bfimat(j,k).gt.1)then
do k1=1,bfimat(j,k)-1
l2=l2+1
freezorb(l2)=k1
enddo
endif
if(bfimat(j,k+1).ne.0.or.k.ne.1)then
if(abs(bfimat(j,k+1)-bfimat(j,k)).ne.1)then
if(abs(bfimat(j,k+1)-bfimat(j,k)).ne.bfimat(j,k+1))then
do k1=bfimat(j,k)+1,bfimat(j,k+1)-1
l2=l2+1
freezorb(l2)=k1
enddo
endif
endif
endif
enddo

frzn=l2

endif
if(sl1(i).eq.ifrag)then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! reading of fragment section !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
k6=0
k7=0
k8=0
k9=0
do j=1,20
gap1(j)=0
enddo
ll=0
l=0
do j=sl(i)+1,sl(i+1)-1
ll=len(trim(charst(j)))
line3=trim(charst(j)(1:1))
if(ll.eq.0)goto 310
if(trim(line3).eq.'#')goto 310
l=l+1
if(l.eq.1)goto 310
k2=0

do k=1,ll+1
line3=trim(charst(j)(k:k))
if(line3.eq.';')goto 311
if(line3.eq.'')then
k2=k2+1
gap2(k2)=k
endif
enddo
311 if(gap2(1).eq.1)then

if(key_frag.eq.1)then
do k=2,k2

if(abs(gap2(k)-gap2(k-1)).ne.1)then
line3=trim(charst(j)(gap2(k-1)+1:gap2(k)-1))

do l1=orbsl-nao+1,orbsl
if(orbs(l1,1).eq.l-1)then
ll=len(trim(line3))
do k3=1,ll
line2=trim(line3(k3:k3))
if(line2.eq.'s')then
k6=k6+1
orbsym(4,k6)=l1-1
goto 208
endif
enddo
do k3=1,ll
line2=trim(line3(k3:k3+1))
if(line2.eq.'px')then
k7=k7+1
orbsym(1,k7)=l1-1
goto 208
endif
enddo
do k3=1,ll
line2=trim(line3(k3:k3+1))
if(line2.eq.'py')then
k8=k8+1
orbsym(2,k8)=l1-1
goto 208
endif
enddo
do k3=1,ll
line2=trim(line3(k3:k3+1))
if(line2.eq.'pz')then
k9=k9+1
orbsym(3,k9)=l1-1
goto 208
endif
enddo
endif
208 enddo
do l1=orbsl-nao+1,orbsl
if(orbs(l1,1).eq.l-1)then
l3=0
ll=len(trim(line3))
do k3=1,ll
line2=trim(line3(k3:k3))
if(line2.eq.'s')then
l3=l3+1
norbsym1(l1-1,l3)='S'
goto 218
endif
enddo
218 do k3=1,ll
line2=trim(line3(k3:k3+1))
if(line2.eq.'px')then
l3=l3+1
norbsym1(l1-1,l3)='X'
goto 219
endif
enddo
219 do k3=1,ll
line2=trim(line3(k3:k3+1))
if(line2.eq.'py')then
l3=l3+1
norbsym1(l1-1,l3)='Y'
goto 220
endif
enddo
220 do k3=1,ll
line2=trim(line3(k3:k3+1))
if(line2.eq.'pz')then
l3=l3+1
norbsym1(l1-1,l3)='Z'
goto 221
endif
enddo
endif
221 num_norbsym1(l1-1)=l3
enddo
goto 210
endif
enddo
endif

210 k5=0
do k=2,k2-1
if(abs(gap2(k)-gap2(k-1)).ne.1)then
line3=trim(charst(j)(gap2(k)+1:gap2(k+1)-1))
ll=len(trim(line3))
do k3=1,ll
line2=trim(line3(k3:k3))
if(line2.eq.'-')then
line5=trim(charst(j)(gap2(k)+1:gap2(k)+k3-1))
line6=trim(charst(j)(gap2(k)+k3+1:gap2(k+1)-1))
read(line5,'(I10)')orbn1
read(line6,'(I10)')orbn2
do k4=orbn1,orbn2
k5=k5+1
frag_at(l,k5)=k4
enddo
goto 454
endif
enddo
k5=k5+1
read(line3,'(I10)')frag_at(l,k5)
endif
454 enddo
frag_atn(l)=k5
endif
if(gap2(1).ne.1)then
if(key_frag.eq.1)then
line3=trim(charst(j)(1:gap2(1)-1))

do l1=orbsl-nao+1,orbsl
if(orbs(l1,1).eq.l-1)then
ll=len(trim(line3))
do k3=1,ll
line2=trim(line3(k3:k3))
if(line2.eq.'s')then
k6=k6+1
orbsym(4,k6)=l1-1
goto 209
endif
enddo
do k3=1,ll
line2=trim(line3(k3:k3+1))
if(line2.eq.'px')then
k7=k7+1
orbsym(1,k7)=l1-1
goto 209
endif
enddo
do k3=1,ll
line2=trim(line3(k3:k3+1))
if(line2.eq.'py')then
k8=k8+1
orbsym(2,k8)=l1-1
goto 209
endif
enddo
do k3=1,ll
line2=trim(line3(k3:k3+1))
if(line2.eq.'pz')then
k9=k9+1
orbsym(3,k9)=l1-1
goto 209
endif
enddo
endif
209 enddo
do l1=orbsl-nao+1,orbsl
if(orbs(l1,1).eq.l-1)then
l3=0
ll=len(trim(line3))
do k3=1,ll
line2=trim(line3(k3:k3))
if(line2.eq.'s')then
l3=l3+1
norbsym1(l1-1,l3)='S'
goto 318
endif
enddo
318 do k3=1,ll
line2=trim(line3(k3:k3+1))
if(line2.eq.'px')then
l3=l3+1
norbsym1(l1-1,l3)='X'
print*,'norbsym1(l1-1,l3)X',l1-1,norbsym1(l1-1,l3)
goto 319
endif
enddo
319 do k3=1,ll
line2=trim(line3(k3:k3+1))
if(line2.eq.'py')then
l3=l3+1
norbsym1(l1-1,l3)='Y'
goto 320
endif
enddo
320 do k3=1,ll
line2=trim(line3(k3:k3+1))
if(line2.eq.'pz')then
l3=l3+1
norbsym1(l1-1,l3)='Z'
goto 321
endif
enddo
endif
321 num_norbsym1(l1-1)=l3
enddo
endif

k5=0
do k=2,k2
line3=trim(charst(j)(gap2(k-1)+1:gap2(k)-1))
if(trim(line3).eq.'') goto 404
ll=len(trim(line3))
do k3=1,ll
line2=trim(line3(k3:k3))
if(line2.eq.'-')then
line5=trim(charst(j)(gap2(k-1)+1:gap2(k-1)+k3-1))
line6=trim(charst(j)(gap2(k-1)+k3+1:gap2(k)-1))
read(line5,'(I10)')orbn1
read(line6,'(I10)')orbn2
do k4=orbn1,orbn2
k5=k5+1
frag_at(l,k5)=k4
enddo
goto 455
endif
 enddo
k5=k5+1
404 read(line3,'(I10)')frag_at(l,k5)
455 enddo
frag_atn(l)=k5
endif
310 enddo
rewind(21)
nfrag=l

norbsym(4)=k6
norbsym(1)=k7
norbsym(2)=k8
norbsym(3)=k9


endif

enddo


vacorb=nao-nae
nlp=nae-nao
if(vacorb.gt.1)nlp=0
nlast=(mult-1)

call read_info

500 return
end subroutine read_xmi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
