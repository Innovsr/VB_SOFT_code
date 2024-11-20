!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_symm_xmi_new(nl,strn,str3,ncqs,q_fac2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

common/quality/str_quality_1,str_quality_2,bondq,tqlty,bqlty,sqlty,tnqs,nssym,qulsym,symq,&
sigsym,tnqs_sig
common/infosymm/set_number,hqlty,ttqlty0,ttqlty1,ttqlty2,ttqlty3,ncqss,qul,Rid&
,mns,u1,max_set,rumset,nqset,strset,strnn,totstr,strno,incmplt,mincmplt,mincmplt_set,group_num
common/str/str5,nstr7

integer::nl,strn,strnn,ncqs,ncqss,tostr,initstr,i,i1,i2,i3,i4,i5,i6,i7,i8,i9,m119,m18,m19,m20,m21,m23,m24,count&
,qul(100),nqul,j,jj,jjj,fg,flg,ii5,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,x,y,strs
integer::i5up16,i16i7,m16m21,i5up15,i15i7,m15m21,i5up14,i14i7,m14m21,i5up13,i13i7,m13m21,i5up12,i12i7,m12m21,&
i5up11,i11i7,m11m21,i5up10,i10i7,m10m21,i5up9,i9i7,m9m21,i5up8,i8i7,m8m21,i5up7,i7i7,m7m21,i5up6,i6i7,m6m21,&
i5up5,i5i7,m5m21,i5up4,i4i7,m4m21,i5up3,i3i7,m3m21,i5up2,i2i7,m2m21,i5up1,i1i7,m1m21,i5up17,i17i7,m17m21
integer::str3(15000,20),q_fac2(15000),finalvec(15000),strset(1000),col(1000),sigsym(15000),tnqs_sig,&
ffvec2(15000,1000),bondq(15000),bondq4(15000),nqset(15000),str5(2000,20),nstr7,qfac3(15000),group_num(15000)&
,tndet,totstr,Ifail,indpnt,strno(1000),str_quality_1(15000),str_quality_2(15000),ttqlty0,ttqlty&
,tqlty,bqlty,sqlty,hqlty,tnqs,nssym,qulsym(15000),symq(15000),set_number,ttqlty1,det_inv,ttqlty2,ttqlty3
integer::rumer(15000),rumer_rad(15000),quality_fac(15000),rumset,u1,max_set,mns,mincmplt,mincmplt_set(500)
!integer,dimension(:),allocatable::qq1,qq2,qq
integer::qq1(5000),qq2(5000),qq(5000),str2(2000,20),Rid,set_num(100),sf1,sf2,incmplt
real*8::ovlp
Double Precision::D(1000)
character(10)::dd,a
character(len=100)::outfile


print*,'enter_symm_xmi_new'
Rid=0
mincmplt=0
do i=1,500
mincmplt_set(i)=0
enddo
ncqss=ncqs
strnn=strn
incmplt=1
mns=0
u1=1
max_set=75000

if(nfset.eq.3.or.nfset.eq.5)then
rumset=0
call rumer_structures(nl,str3,ncqss,rumer,rumer_rad)
call write_rumer_xmi(nl,str3,ncqss,rumer,rumer_rad,quality_fac)
endif

set_number=0
bqlty=0
tqlty=0
sqlty=0
hqlty=0
indpnt=2
!ovlpval=1.0
if(noq0.gt.strnn)then
ttqlty0=noq0
else
ttqlty0=strnn+noq0
endif
ttqlty1=strnn+noq1
ttqlty2=strnn+noq2
ttqlty3=strnn+noq3

write(*,*)'sl  structures           group_numbers'
do i=1,ncqss
write(*,231)i,(str3(i,j),j=1,nae),q_fac2(i)
!write(*,231),i,(str3(i,j),j=1,nae),q_fac2(i),str_quality_1(i),str_quality_2(i),bondq(i)
group_num(i)=q_fac2(i)
enddo
231 format(30I3)

jj=1
do m19=1,ncqss
!print*,'q_fac2',q_fac2(m19)
if(m19.eq.1)qul(1)=q_fac2(1)
j=jj
do i=1,j
if(qul(i).eq.q_fac2(m19))goto 373
enddo
jj=jj+1
qul(i)=q_fac2(m19)
!print*,qul(i)
373 enddo
nqul=jj
!print*,'nqul',nqul

do i=1,nqul
jjj=0
jj=0
do m19=1,ncqss
!print*,qul(i),q_fac2(m19)
if(qul(i).eq.q_fac2(m19))then
jjj=jjj+1
jj=m19
endif
enddo
nqset(i)=jjj
strset(i)=jj
!print*,'i,nqset(i),strset(i)',i,nqset(i),strset(i)
enddo

flg=0
totstr=0
i7=0
m21=0
i1i7=0
m1m21=0
!i5=ndet
!do i8=1,1000
!finalvec(i8)=0
!enddo
!do m19=1,10000
!do m18=1,15
!str2(m19,m18)=0
!enddo
!enddo


!i5up1=i5
i1i7=i7
m1m21=m21


do m1=1,nqul
!i5=i5up1
!i7=i1i7
!m21=m1m21
i7=0
m21=0

!print*,'loop1'
call write_symm_xmi_1(i7,m21,m1,sf1,sf2,str3,q_fac2)
if (sf1.eq.1)goto 701
if (sf2.eq.1)goto 100

!print*,'i7,m21**',i7,m21
!i5up2=i5
i2i7=i7
m2m21=m21

do m2=m1+1,nqul
!i5=i5up2
i7=i2i7
m21=m2m21

!print*,'loop2'
call write_symm_xmi_1(i7,m21,m2,sf1,sf2,str3,q_fac2)
if (sf1.eq.1)goto 702
if (sf2.eq.1)goto 100

!print*,'22 i7,m21**',i7,m21
!i5up3=i5
i3i7=i7
m3m21=m21

do m3=m2+1,nqul
!i5=i5up3
i7=i3i7
m21=m3m21

!print*,'loop3'
call write_symm_xmi_1(i7,m21,m3,sf1,sf2,str3,q_fac2)
if (sf1.eq.1)goto 703
if (sf2.eq.1)goto 100

!print*,'33 i7,m21**',i7,m21
!stop
!i5up4=i5
i4i7=i7
m4m21=m21

do m4=m3+1,nqul
!i5=i5up4
i7=i4i7
m21=m4m21

!print*,'loop4'
call write_symm_xmi_1(i7,m21,m4,sf1,sf2,str3,q_fac2)
if (sf1.eq.1)goto 704
if (sf2.eq.1)goto 100

!i5up5=i5
i5i7=i7
m5m21=m21

do m5=m4+1,nqul
!i5=i5up5
i7=i5i7
m21=m5m21

!print*,'loop5'
call write_symm_xmi_1(i7,m21,m5,sf1,sf2,str3,q_fac2)
if (sf1.eq.1)goto 705
if (sf2.eq.1)goto 100

!i5up6=i5
i6i7=i7
m6m21=m21

do m6=m5+1,nqul
!i5=i5up6
i7=i6i7
m21=m6m21

!print*,'loop6'
call write_symm_xmi_1(i7,m21,m6,sf1,sf2,str3,q_fac2)
!print*,'sf1,sf2',sf1,sf2
if (sf1.eq.1)goto 706
if (sf2.eq.1)goto 100

!i5up7=i5
i7i7=i7
m7m21=m21

do m7=m6+1,nqul
!i5=i5up7
i7=i7i7
m21=m7m21

!print*,'loop7'
call write_symm_xmi_1(i7,m21,m7,sf1,sf2,str3,q_fac2)
if (sf1.eq.1)goto 707
if (sf2.eq.1)goto 100

!i5up8=i5
i8i7=i7
m8m21=m21

do m8=m7+1,nqul
!i5=i5up8
i7=i8i7
m21=m8m21

!print*,'loop8'
call write_symm_xmi_1(i7,m21,m8,sf1,sf2,str3,q_fac2)
if (sf1.eq.1)goto 708
if (sf2.eq.1)goto 100

!i5up9=i5
i9i7=i7
m9m21=m21

do m9=m8+1,nqul
!i5=i5up9
i7=i9i7
m21=m9m21

!print*,'loop9'
call write_symm_xmi_1(i7,m21,m9,sf1,sf2,str3,q_fac2)
if (sf1.eq.1)goto 709
if (sf2.eq.1)goto 100

!i5up10=i5
i10i7=i7
m10m21=m21

do m10=m9+1,nqul
!i5=i5up10
i7=i10i7
m21=m10m21

!print*,'loop10'
call write_symm_xmi_1(i7,m21,m10,sf1,sf2,str3,q_fac2)
if (sf1.eq.1)goto 710
if (sf2.eq.1)goto 100

!i5up11=i5
i11i7=i7
m11m21=m21

do m11=m10+1,nqul
!i5=i5up11
i7=i11i7
m21=m11m21

!print*,'loop11'
call write_symm_xmi_1(i7,m21,m11,sf1,sf2,str3,q_fac2)
if (sf1.eq.1)goto 711
if (sf2.eq.1)goto 100

!i5up12=i5
i12i7=i7
m12m21=m21

do m12=m11+1,nqul
!i5=i5up12
i7=i12i7
m21=m12m21

print*,'loop12'
call write_symm_xmi_1(i7,m21,m12,sf1,sf2,str3,q_fac2)
if (sf1.eq.1)goto 712
if (sf2.eq.1)goto 100

!i5up13=i5
i13i7=i7
m13m21=m21

do m13=m12+1,nqul
!i5=i5up13
i7=i13i7
m21=m13m21

print*,'loop13'
call write_symm_xmi_1(i7,m21,m13,sf1,sf2,str3,q_fac2)
if (sf1.eq.1)goto 713
if (sf2.eq.1)goto 100

!i5up14=i5
i14i7=i7
m14m21=m21

do m14=m13+1,nqul
!i5=i5up14
i7=i14i7
m21=m14m21

!print*,'loop14',m14,sf1,sf2
!stop
call write_symm_xmi_1(i7,m21,m14,sf1,sf2,str3,q_fac2)
if (sf1.eq.1)goto 714
if (sf2.eq.1)goto 100

!i5up15=i5
i15i7=i7
m15m21=m21

do m15=m14+1,nqul
!i5=i5up15
i7=i15i7
m21=m15m21

print*,'loop15'
call write_symm_xmi_1(i7,m21,m15,sf1,sf2,str3,q_fac2)
if (sf1.eq.1)goto 715
if (sf2.eq.1)goto 100

!i5up16=i5
i16i7=i7
m16m21=m21

do m16=m15+1,nqul
!i5=i5up16
i7=i16i7
m21=m16m21

print*,'loop16'
call write_symm_xmi_1(i7,m21,m16,sf1,sf2,str3,q_fac2)
if (sf1.eq.1)goto 716
if (sf2.eq.1)goto 100

!i5up16=i5
i16i7=i7
m16m21=m21

716 enddo

715 enddo

714 enddo

713 enddo

712 enddo

711 enddo

710 enddo

709 enddo

708 enddo

707 enddo

706 enddo

705 enddo

704 enddo

703 enddo

702 enddo

701 enddo


close(21)

if(incmplt.eq.1)then
write(10,*)'----------- incomplete set -----------'
do i=1,mincmplt
 qq(i)=q_fac2(mincmplt_set(i))
 qq1(i)=str_quality_1(mincmplt_set(i))
 qq2(i)=str_quality_2(mincmplt_set(i))
 bondq4(i)=bondq(mincmplt_set(i))
    if(niao.eq.0)then
     write(10,900)qq1(i),bondq4(i),qq2(i),qq(i),'|',(str3(mincmplt_set(i),m20),m20=1,nae)
    endif
    if(niao.gt.1)then
     write(10,901)qq1(i),bondq4(i),qq2(i),qq(i),'|',1,':',niao,(str3(mincmplt_set(i),m20),m20=1,nae)
    endif
    if(niao.eq.1)then
     write(10,909)qq1(i),bondq4(i),qq2(i),qq(i),'|',1,1,(str3(mincmplt_set(i),m20),m20=1,nae)
    endif
enddo
endif

900 format(I3,x,I3,x,I3,x,I3,x,a,x,25I4)
901 format(I3,x,I3,x,I3,x,I3,x,a,x,I1,a,I3,x,25I4)
909 format(I3,x,I3,x,I3,x,I3,x,a,x,I3,I3,x,25I4)

open(unit=121,file='script10',status='unknown')
write(121,111)'rm -rf','Rumer_Sets.dat'
close(121)
CALL SYSTEM ("chmod +x script10 ")
CALL SYSTEM ("./script10 ")
CALL SYSTEM ("rm script10 ")
111 format(a,x,a)
100 return
end subroutine write_symm_xmi_new
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
