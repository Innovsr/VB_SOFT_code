!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine symmetry_cal_sig(nl,str1,ncqs,stsymsc,symq,nssym)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
use commondat1
implicit none

integer::m19,m18,i,i1,i2,i3,i4,i5,i6,i7,i8,ii,iii,iiii,nl,ncqs,k,k1,k2,n1,n2,n3,n4,n5,n6,l1,l3,l6,&
nd,fullgrp,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,&
k25,k26,k27,n,j1,j2,j3,j4,j5,j6,j7,j8,j9,j10,nssym,j,jj,nn1,nn2,kkk,n17,n18,n9,ii2,nnat_bond_new_1(100,2)
integer::nnat_bond_new(100,2),tot_orb(100),nn_count(50),sig_orb(100),sig_orb_1(100),nsig,nnatom,natom&
,nn4,nn5,nna
real*8::tscore,score(100,2),new_score(100,2),stsymsc(15000),ssym(15000),order(15000),lpna,&
losc,lone_score(15000),symq1(5000),cot(1000),sig_sym_sc(15000),sym_sets(15000),n7,n8,bd_dist
integer::kk,nscore(100),nnscore,ipnum,loop_score,loopsc,ncrow
!double precision::rscore,strscore
real*8::rscore,strscore,scr,bndscore(20),coord_score
!real::nnscore,rscore,strscore
integer::nn(10),str1(15000,20),mbondq(15000),nn_group(50,10),nelimt(50),&
full_nn_group(1000),sl_group(50,10),symq(15000),connectivity_row(20,20)
!integer::kval(1000),deadend(20)
integer::atsymset(20,20),nsym,syn(50),at_sym(50),numbond,numbond1,loop_score_row(20),coordination_mat(20,20)
real*8::atoset_symscr(1000),piscr
real::stsymsc1(15000)

common/ats/atsymset,nsym,syn,at_sym
common/loops/full_nn_group,fullgrp,natom,nelimt,sl_group,tot_orb,nn_group
print*,'enter symmetry_cal_sig'
!! iabd--> inactive bonds associated with the atoms
!! ialp--> inactive lone pairs associated with the atoms
!! score()--> scoring of the atoms depending on the the 'iabd', 'ialp' and atomic number (at_number)
open(unit=15,file='symsig.temp',status='unknown')

do k1=1,nsym
print*,'at_sym(k1)',at_sym(k1)
enddo



!do i3=1,atom
!print*,'atoset_symm_sig',(atoset(active_atoms(i3),i4),i4=1,atn(active_atoms(i3)))
!print*,'active_atoms',active_atoms(i3)
!enddo


do i=1,100
do i1=1,2
score(i,i1)=0.0
enddo
enddo

print*,'prime_num',prime_num(2)
i4=0
do i3=1,atom
!lpna=prime_num(1)
lpna=0.0
k=i3
!k=active_atom(i3)
if(input_flg.eq.1)then
do i1=1,niabd
if(k.eq.iabd(i1))then
lpna=lpna+1.0/prime_num(2)
print*,'lpna',lpna,prime_num(2)
endif
enddo
!print*,'rajat roy',nialp
do i2=1,nialp
if(k.eq.ialp(i2))then
!print*,'rajat'
lpna=lpna+1/prime_num(1)
endif
enddo
do i2=1,niach
if(k.eq.iach(i2))then
!print*,'rajat'
lpna=lpna+1/prime_num(3)
endif
enddo
else

!lpna=lpna+biasval(i3)+dist_nnat(i3)+atn(active_atoms(i3))/prime_num(1)
lpna=lpna+biasval(i3)+dist_nnat(i3)
endif
!score(i3,2)=lpna+1.0/prime_num(at_num(i3))
!score(i3,2)=lpna+1.0/at_num(i3)
score(i3,2)=lpna+at_num(i3)
print*,'score:symme',score(i3,2),at_num(i3),biasval(i3),dist_nnat(i3)
enddo



print*,'symmetry_cal',nnnatom
n6=0
n7=0
n8=0
do i2=1,nnnatom
do i1=1,2
nnat_bond_new(i2,i1)=nnat_bond(i2,i1)+nao+niao
enddo
print*,(nnat_bond_new(i2,i1),i1=1,2),'|',(nnat_bond(i2,i1),i1=1,2)
enddo
n5=0
do i2=1,nnnatom
do i1=1,2
do i3=1,n5
if(nnat_bond_new(i2,i1).eq.tot_orb(i3)) goto 537
enddo
n5=n5+1
tot_orb(n5)=nnat_bond_new(i2,i1)
537 enddo
enddo
!n5=0
!do i1=1,atom
!n5=n5+1
!tot_orb(n5)=i1+nao+niao
!enddo

do i1=1,atom
n1=0
do i3=1,nnnatom
do i4=1,2
if(active_atoms(i1).eq.nnat_bond(i3,i4))then
!if(i1.eq.nnat_bond(i3,i4))then
n1=n1+1
endif
enddo
enddo
nn_count(i1)=n1
enddo
do i3=1,atom
print*,'nn_count',nn_count(i3)
enddo



! do i3=1,atom
!print*,'active_atom',i3,active_atoms(i3)
! do i4=1,atn(active_atoms(i3))
!print*,'atn',i4,atn(active_atoms(i3))
!print*,'atn(i3),syn(2)',syn(nsym),(atsymset(nsym,k2),k2=1,syn(nsym))
!enddo
!enddo

do i3=1,atom
print*,'atoset',(atoset(i3,i4),i4=1,atn(i3))
enddo


print*,'atom:atom',atom
n2=nnnatom
n4=0
n3=0
n6=0
 do i3=1,atom
print*,'iiiiiiiiiiiiiiiiiiii33333333333',i3
n1=0
do i4=1,atn(active_atoms(i3))
! do i4=1,atn(i3)
!print*,'atn**',i4,atn(active_atoms(i3))
 do k2=1,syn(nsym)
!print*,'atn(i3),syn(2)',atn(active_atoms(i3)),syn(nsym),atsymset(nsym,k2)
 if(atoset(active_atoms(i3),i4).eq.atsymset(nsym,k2))then
! if(atoset(i3,i4).eq.atsymset(nsym,k2))then
n1=n1+1
n6=n6+1
!sig_orb(n1)=atoset(i3,i4)
sig_orb(n1)=atoset(active_atoms(i3),i4)
sig_orb_1(n6)=sig_orb(n1)
print*,'sig_orb(n1)',sig_orb(n1),n1
 endif
 enddo
 enddo
!!nsig=n1
 if(n1.gt.0)then
n3=n3+1
k=0
!if(nn_count(i3).eq.1)k=1
!print*,'k',k,i3,nn_count(i3)
 do i2=1,nnnatom
 do i1=1,2
 if(i3.eq.nnat_bond(i2,i1))then
k=k+1
if(k.gt.n1)goto 113
n4=n4+1
!sig_orb_1(n4)=sig_orb(k)
nnat_bond_new(i2,i1)=sig_orb(k)
print*,'nnat_bond_new(i2,i1)***',i2,i1,nnat_bond(i2,i1),nnat_bond_new(i2,i1),sig_orb(k)
 endif
 enddo
113 enddo
!print*,'n1',n1,niao+nao+n3,n3
 do i1=1,n1
n2=n2+1
n5=n5+1
nnat_bond_new(n2,1)=niao+nao+i3
nnat_bond_new(n2,2)=sig_orb(i1)
tot_orb(n5)=sig_orb(i1)
 enddo
if(n1.gt.2)then
 do i1=1,n1
if(i1.ne.n1)ii2=i1+1
if(i1.eq.n1)ii2=1
n2=n2+1
nnat_bond_new(n2,1)=sig_orb(i1)
nnat_bond_new(n2,2)=sig_orb(ii2)
 enddo
endif
endif
 enddo
!!!!!!!!!!!!
nnatom=n2
nsig=n6
!!!! nsig=number of sigma orbital
!!!! nnatom= number of connected pair of orbitals in the extended CM
print*,'nsig',nsig,nnatom

natom=n5
!!!! natom=number of extended orbitals

do i2=1,nnnatom
print*,'nnat_bond',(nnat_bond(i2,i3),i3=1,2)
enddo
do i2=1,nnatom
print*,'nnat_bond_new',(nnat_bond_new(i2,i3),i3=1,2),nnatom
enddo
!print*,'sig_orb_1(n4)',(sig_orb_1(i2),i2=1,n4)
do i2=1,n5
print*,'tot_orb',tot_orb(i2),n5,nnatom
enddo

call nnat_bond_sig(nnat_bond_new,nnatom,sig_orb_1,nsig,nnat_bond_new_1,nna)
do i2=1,nna
print*,'nnat_bond_new_new',(nnat_bond_new_1(i2,i3),i3=1,2),nna
enddo

!231 format(20I3)
m18=0
m19=0
do i3=1,n5
k=tot_orb(i3)
m19=0
!do i2=1,nnatom
do i2=1,nna
do i1=1,2
!print*,'atom',k,nnat_bond_new(i2,i1)
!if(k.eq.nnat_bond_new(i2,i1))then
if(k.eq.nnat_bond_new_1(i2,i1))then
if(i1.eq.1)i5=2
if(i1.eq.2)i5=1
m19=m19+1
m18=m18+1
nn_group(i3,m19)=nnat_bond_new_1(i2,i5)
!print*,'nnat_bond(i2,i5)',i3,k,i2,i1,i5,m19,nnat_bond_new(i2,i5),nn_group(i3,m19)
sl_group(i3,m19)=m18
endif
enddo
enddo
nelimt(i3)=m19
!print*,'nelimt(i3)',nelimt(i3)
enddo


!!! 'nn_group' specify orbitals in tot_orb array are connected with which other
!!! orbitals

do i3=1,n5
print*,'*nn_group*',nelimt(i3),(nn_group(i3,i4),i4=1,nelimt(i3))
enddo
do i3=1,n5
print*,'*sl_group*',(sl_group(i3,i4),i4=1,nelimt(i3))
enddo

m19=0
do i3=1,n5
do i4=1,nelimt(i3)
m19=m19+1
full_nn_group(m19)=nn_group(i3,i4)
enddo
enddo
fullgrp=m19

write(*,*)'full_nn_g',m19,(full_nn_group(i3),i3=1,m19)
!goto 200
!endif
!!!!! Scoring of the structures has been started from here !!!

do i=1,ncqs
sig_sym_sc(i)=1.0
enddo


do i=1,ncqs
loopsc=0
do i3=1,atom
new_score(i3,2)=score(i3,2)
enddo

!! if the structures have active lone pairs and or radicals the associated atoms are being scored in 'new_score()'.

!print*,'symmetry_cal_sym_cross'
write(*,231)(str1(i,i1),i1=1,nae)
!print*,'tscore**,nl',nl
!!if(nl.ne.0)then
!!!print*,'rajat',nl
!!do i2=1,nl*2,2
!!!print*,i2
!!do i3=1,atom
!!do i4=1,atn(active_atoms(i3))
!!!print*,'atoset(active_atoms(i3),i4)',atoset(active_atoms(i3),i4)
!!if(str1(i,i2).eq.atoset(active_atoms(i3),i4))then
!!
!!new_score(i3,2)=new_score(i3,2)+1.0/prime_num(2)
!!!print*,'new score2',new_score(i3,2)
!!endif
!!enddo
!!enddo
!!enddo
!!endif
!!
!!if(nlast.ne.0)then
!!!print*,'rajat',nlast,nae-nlast+1,atom
!!do i2=nae-nlast+1,nae
!!do i3=1,atom
!!do i4=1,atn(active_atoms(i3))
!!if(str1(i,i2).eq.atoset(active_atoms(i3),i4))then
!!
!!new_score(i3,2)=new_score(i3,2)+1.0/prime_num(3)
!!endif
!!enddo
!!!print*,'new score3',i3,new_score(i3,2)
!!enddo
!!enddo
!!endif


!! Scoring of the bonds of the structures started from here
!! the score of the bond is given by the score of the atoms involeved in the


!bond
!do i3=1,atom
!do i4=1,atn(active_atoms(i3))
!!print*,'atoset(active_atoms(i3),i4)',atoset(active_atoms(i3),i4)
!enddo
!enddo
! do i3=1,nsig
!  print*,'sig_orb_1(i3)',sig_orb_1(i3),nsig
!enddo

strscore=0.0
ipnum=3
!print*,'strset'
numbond=0
numbond1=0
do i1=1+nl*2,nae-nlast,2
tscore=0.0
rscore=0.0
n1=0
n2=0
n3=0
n4=0
do i2=i1,i1+1
do i3=1,atom
do i4=1,atn(active_atoms(i3))
!do i4=1,atn(i3)
if(str1(i,i2).eq.atoset(active_atoms(i3),i4))then
!if(str1(i,i2).eq.atoset(i3,i4))then
if(i2.eq.i1) then
n1=active_atoms(i3)
!n1=i3
endif
if(i2.eq.i1+1) then
n2=active_atoms(i3)
!n2=i3

endif
endif
enddo
enddo
!print*,'bond***',str1(i,i2),n1,n2
enddo
n4=str1(i,i1)
n5=str1(i,i1+1)
do i2=1,nsym
do i3=1,syn(i2)
!print*,'atsymset',i2,i3,n4,n4,atsymset(i2,i3)
if(n4.eq.atsymset(i2,i3))nn4=i2
if(n5.eq.atsymset(i2,i3))nn5=i2
enddo
enddo
!print*,'nn4,nn5',nn4,nn5,n4,n5
if(nsym.eq.1)then
if(symtype.eq.'tight'.or.symtype.eq.'loose')then
if(nn4.eq.1.and.nn5.eq.1)scr=1.0
endif
endif
if(nsym.eq.2)then
if(symtype.eq.'tight'.or.symtype.eq.'loose')then
if(nn4.eq.2.and.nn5.eq.2)scr=1.0
if(nn4.eq.1.and.nn5.eq.1)scr=2.0
if(nn4.eq.1.and.nn5.eq.2)scr=4.0
if(nn4.eq.2.and.nn5.eq.1)scr=4.0
endif
endif
if(nsym.eq.3)then
if(symtype.eq.'tight')then
if(nn4.eq.3.and.nn5.eq.3)scr=1.0
if(nn4.eq.1.and.nn5.eq.1)scr=2.0
if(nn4.eq.2.and.nn5.eq.2)scr=2.0
if(nn4.eq.1.and.nn5.eq.2)scr=3.0
if(nn4.eq.2.and.nn5.eq.1)scr=3.0
if(nn4.eq.1.and.nn5.eq.3)scr=4.0
if(nn4.eq.3.and.nn5.eq.1)scr=4.0
if(nn4.eq.2.and.nn5.eq.3)scr=5.0
if(nn4.eq.3.and.nn5.eq.2)scr=5.0
endif
if(symtype.eq.'loose')then
if(nn4.eq.3.and.nn5.eq.3)scr=1.0
if(nn4.eq.1.and.nn5.eq.1)scr=2.0
if(nn4.eq.2.and.nn5.eq.2)scr=2.0
if(nn4.eq.1.and.nn5.eq.2)scr=2.0
if(nn4.eq.2.and.nn5.eq.1)scr=2.0
if(nn4.eq.1.and.nn5.eq.3)scr=3.0
if(nn4.eq.3.and.nn5.eq.1)scr=3.0
if(nn4.eq.2.and.nn5.eq.3)scr=3.0
if(nn4.eq.3.and.nn5.eq.2)scr=3.0
endif
endif
!tscore=n4*new_score(n1,2)+n5*new_score(n2,2)
tscore=new_score(n1,2)+new_score(n2,2)+scr/(new_score(n1,2)+new_score(n2,2))
!print*,'tscore,scr,OSS,ABS(topology)',tscore,scr,scr/(new_score(n1,2)+new_score(n2,2)),new_score(n1,2)+new_score(n2,2)
nnscore=0.0
nd=0
n1=0
n2=0
n3=0
n4=0
!do i3=1,nsig
!print*,'sig_orb_1',sig_orb_1(i3)
!enddo
 do i2=i1,i1+1
!print*,'i2',i2,str1(i,i2)
 do i3=1,nsig
  if(sig_orb_1(i3).eq.str1(i,i2))then
n3=n3+1
if(i2.eq.i1)n1=str1(i,i2)
if(i2.eq.i1+1)n2=str1(i,i2)
!print*,'n1,n2',n1,n2,n3,str1(i,i2)
goto 232
  endif
 enddo
 do i3=1,atom
 do i4=1,atn(active_atoms(i3))
! do i4=1,atn(i3)
  if(str1(i,i2).eq.atoset(active_atoms(i3),i4))then
!  if(str1(i,i2).eq.atoset(i3,i4))then
   if(i2.eq.i1) then
n1=i3+nao+niao
   endif
    if(i2.eq.i1+1) then
n2=i3+nao+niao
    endif
  endif
 enddo
 enddo
 232 enddo

!print*,'n3,bons',n3,n1,n2,n5,n18
!print*,'tscore',tscore

!loop_score=1
if(n1.gt.n2)then
k1=n2
k2=n1
else
k1=n1
k2=n2
endif
print*,'k1,k2,bond',k1,k2,natom
do i2=1,natom
if(k1.eq.tot_orb(i2))l1=i2
enddo
!print*,'l1',l1
do i2=1,nelimt(l1)
 if(k2.eq.nn_group(l1,i2))then
n=1
!nscore(n)=prime_num(9)
!nscore(n)=0
!print*,'nscore before loop start',nscore(n)
!goto 101
ncrow=1
loop_score=1
numbond=numbond+1
loop_score_row(numbond)=loop_score+1
coordination_mat(numbond,1)=k1
coordination_mat(numbond,2)=k2
goto 555
 endif
enddo

if(k1.ne.k2)call symm_loops(k1,k2,loop_score,connectivity_row,ncrow)
do i3=1,ncrow
numbond=numbond+1
loop_score_row(numbond)=loop_score+1
do i2=1,loop_score+1
coordination_mat(numbond,i2)=connectivity_row(i3,i2)
enddo
print*,'coordination_mat',(coordination_mat(numbond,i2),i2=1,loop_score+1)
enddo
if(k1.eq.k2)loop_score=0
goto 555
!555 print*,'tscore_before',tscore
555 if(loop_score.ne.0)tscore=tscore+1.0/loop_score
!print*,'loop_score',loop_score
!print*,'kk*,tscore',tscore,loop_score
!endif
rscore=1.0/tscore
do i2=1,ncrow
numbond1=numbond1+1
bndscore(numbond1)=rscore
enddo
!print*,'nearest neighbor is',kk,tscore,rscore
strscore=strscore+rscore
!strscore=strscore+tscore
if(loop_score.ne.0)loopsc=loopsc+prime_num(loop_score)
if(loop_score.eq.0)loopsc=loopsc+1000
!print*,strscore
!print*,'strscore,tscore,loop_score,loopsc',strscore,tscore,prime_num(loop_score),loopsc
enddo
print*,'numbond',numbond
do i1=1,numbond
print*,'loop_score_row',loop_score_row(i1)
print*,'coordination_mat',(coordination_mat(i1,i2),i2=1,loop_score_row(i1))
enddo

call coordination_val(coordination_mat,numbond,loop_score_row,bndscore,coord_score)

stsymsc(i)=strscore+coord_score
!stsymsc(i)=strscore
loopsymsc(i)=loopsc
print*,'stsymsc(i)',i,stsymsc(i),loopsymsc(i),coord_score
write(15,*)stsymsc(i)
enddo



rewind(15)
do i=1,ncqs
read(15,909)stsymsc1(i)
print*,i,stsymsc1(i)
enddo
909 format (F10.6)


jj=1
do m19=1,ncqs
!print*,'q_fac2',q_fac2(m19)
if(m19.eq.1)ssym(1)=stsymsc1(1)
j=jj
do i=1,j
if(ssym(i).eq.stsymsc1(m19))goto 373
enddo
jj=jj+1
ssym(i)=stsymsc1(m19)
!print*,qul(i)
373 enddo
nssym=jj
print*,'nssym',nssym

do i=1,jj
print*,ssym(i)
enddo

do i=1,10000
order(i)=0.0
enddo

do k3=1,nssym
do k4=1,nssym
do k5=1,k3
if(order(k5).eq.ssym(k4))goto 389
enddo
if(ssym(k4).gt.order(k3))then
order(k3)=ssym(k4)
endif
389 enddo
print*,'**order',k3,order(k3)
enddo

jj=0
do i=1,nssym
jj=jj+1
do j=1,ncqs
if(order(i).eq.stsymsc1(j))then
write(*,231)(str1(j,i1),i1=1,nae)
print*,'order(i),stsymsc(j)',order(i),stsymsc1(j)
if(nlast.ne.0)then
losc=0.0
do i1=nae-nlast,nae
do i2=1,syn(1)
if(str1(j,i1).eq.atsymset(1,i2))goto 531
enddo
goto 532
531 losc=losc+1.0/prime_num(str1(j,i1))
532 enddo
lone_score(j)=losc
endif
symq(j)=jj
!lsym(j)=loopsymsc(j)
print*,'symq(j),lsym(j)',symq(j),loopsymsc(j)
endif
enddo
enddo




print*,'sub:symmetry_cal'
do i=1,ncqs
!if(symq(i).eq.4)then
write(*,231)(str1(i,j),j=1,nae),symq(i)
!print*,lone_score(i)
!endif
enddo
print*,'sub:symmetry_cal'
!do i=1,ncqs
!write(*,231)(str2(i,j),j=1,nae)
!enddo

231 format(20I3)
!call sigma_symm(nl,ncqs,str1)
CALL SYSTEM ("rm symsig.temp")

print*,'exit symmetry_cal_sig'
200 return
end subroutine symmetry_cal_sig
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
