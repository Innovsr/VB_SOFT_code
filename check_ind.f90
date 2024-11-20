!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine check_ind(str3,totstr,nl,ncqs,strno,Ifailn,ffvec2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::j,i,ii,iii,str3(15000,20),d,e,elporb,str_rad_1(100),str_rad_2(100),ffvec2(15000,1000),&
Ifailn,i7,wig2,na,i1,i2,i3,totstr,ncqs,strno(15000),str_bnd_1(1000),str_bnd_2(1000),sing(1000,100),&
str_sl(1000),new_mult,nl,totbnd,totseorb,str_bnd(15000,1000),str_rad(500,20),bond_sl(10000,2),&
npstr,num_orb(20),val(20),val1(20),i4,i5,i6,l,l1,rul,bdsl(20),perm(20),Ifailn1
real*8::factorial
character(9)::sbrtn

common /chek/totbnd,totseorb,str_bnd,str_rad,bond_sl,num_orb
!character(5)::a

Ifailn=0

!do i=1,totbnd
!print*,(bond_sl(i,i1),i1=1,2)
!enddo

e=nae-nl*2
call wigner(e,wig2)
npstr=wig2

!print*,'totbnd',totbnd
iii=0

!do i=1,totstr
!write(*,100),'str3',(str3(strno(i),ii),ii=1,nae),strno(i)
!!write(*,100),'fvec',(ffvec2(strno(i),ii)*sing(i,ii),ii=1,32)
!enddo

if(nlast.ne.0)then
do ii=1,totseorb
do i=1,totstr
str_rad_1(ii)=str_rad_1(ii)+str_rad(strno(i),ii)
enddo
enddo
d=0
do ii=1,totseorb
if(str_rad_1(ii).gt.d)d=str_rad_1(ii)
enddo
if(d.lt.3) goto 121
new_mult=mult
e=nae-1-nl*2
mult=mult-1
call wigner(e,wig2)
mult=new_mult
do ii=1,totseorb
if(str_rad_1(ii).gt.wig2)then
Ifailn=1
goto 121
endif
enddo

do ii=3,d-1
do i=1,totseorb
if(ii.eq.str_rad_1(i))then
do i1=1,totstr
if(str_rad(strno(i1),i).eq.1)then
iii=iii+1
str_sl(iii)=strno(i1)
endif
enddo

do i2=1,totseorb
do i1=1,iii
str_rad_2(i2)=str_rad_2(i2)+str_rad(str_sl(i1),i2)
enddo
enddo

i3=0
do i2=1,totseorb
if(str_rad_2(i2).eq.ii)then
i3=i3+1
endif
enddo

new_mult=mult
e=nae-nl*2-i3
mult=(nlast-i3)+1
call wigner(e,wig2)
mult=new_mult
if(ii.gt.wig2)then
Ifailn=1
goto 121
endif
endif
enddo
enddo


endif



do ii=1,totbnd
!print*,'ii',ii
str_bnd_1(ii)=0
do i=1,totstr
str_bnd_1(ii)=str_bnd_1(ii)+str_bnd(strno(i),ii)
enddo
!print*,'str_bnd',str_bnd_1(ii),ii
enddo

d=0
do ii=1,totbnd
if(str_bnd_1(ii).gt.d)d=str_bnd_1(ii)
enddo

!if(totstr.lt.0)then
!do i1=1,nae-nl*2
!do i2=i1+1,nae-nl*2
!do i3=i2+1,nae-nl*2
!val(1)=num_orb(i1)
!val(2)=num_orb(i2)
!val(3)=num_orb(i3)
!
!print*,'orbs',(val(i6),i6=1,3)
!l1=0
!do i4=1,totbnd
!l=0
!do i5=1,2
!do i6=1,3
!if(bond_sl(i4,i5).eq.val(i6))then
!l=l+1
!endif
!enddo
!if(l.eq.2)then
!l1=l1+1
!val1(l1)=i4
!endif
!enddo
!enddo
!
!!print*,'orbnum',(val1(i6),i6=1,l1)
!l=0
!do i4=1,l1
!l=l+str_bnd_1(val1(i4))
!enddo
!rul=int(totstr/3)
!!print*,'lllll',l
!if(l.gt.rul*2)then
!Ifailn=1
!
!print*,'orbnum',rul,l,'>',(val1(i6),i6=1,l1)
!goto 121
!endif
!enddo
!enddo
!enddo
!
!!stop
!
!endif
!!l1=npstr/3



do ii=1,totstr
!write(*,100),'str',(str3(strno(ii),i),i=1,nae)
write(*,100),'str_bnd',strno(ii),(str_bnd(strno(ii),i),i=1,totbnd)
enddo
write(*,100),'add_val',(str_bnd_1(i),i=1,totbnd)
print*,'ddd',d,totbnd
100 format(a,50I3)

if(d.lt.3) goto 121
e=nae-2-nl*2
call wigner(e,wig2)
!print*,'wig2',wig2
do ii=1,totbnd
if(str_bnd_1(ii).gt.wig2)then
Ifailn=1
goto 121
endif
enddo



do i=1,totbnd
if(str_bnd_1(i).gt.2.and.str_bnd(strno(totstr),i).eq.1)then
iii=0
do i1=1,totstr
if(str_bnd(strno(i1),i).eq.1)then
iii=iii+1
str_sl(iii)=strno(i1)
print*,'str_sl',i,'>',iii,str_sl(iii)
endif
enddo

call all_set_gen(str_sl,iii,nl,i,Ifailn)
if(Ifailn.eq.1)goto 121
endif
enddo

!stop

ii=0
do i=1,totbnd
if(str_bnd(strno(totstr),i).eq.1)then
ii=ii+1
bdsl(ii)=i
endif
enddo

sbrtn='check_ind'
do j=2,int(ii/2)

do i1=1,ii
perm(1)=bdsl(i1)

do i2=i1+1,ii
perm(2)=bdsl(i2)
if(j.eq.2)then
print*,(perm(i),i=1,j),j
!stop
call check3(j,perm,strno,totstr,nl,Ifailn1,sbrtn)
if(Ifailn1.eq.1)goto 122
goto 101
endif

do i3=i2+1,ii
perm(3)=bdsl(i3)
if(j.eq.3)then
call check3(j,perm,strno,totstr,nl,Ifailn1,sbrtn)
if(Ifailn1.eq.1)goto 122
goto 102
endif

do i4=i3+1,ii
perm(4)=bdsl(i4)
if(j.eq.4)then
call check3(j,perm,strno,totstr,nl,Ifailn1,sbrtn)
if(Ifailn1.eq.1)goto 122
goto 103
endif

do i5=i4+1,ii
perm(5)=bdsl(i5)
if(j.eq.5)then
call check3(j,perm,strno,totstr,nl,Ifailn1,sbrtn)
if(Ifailn1.eq.1)goto 122
goto 104
endif

do i6=i5+1,ii
perm(6)=bdsl(i6)
if(j.eq.6)then
call check3(j,perm,strno,totstr,nl,Ifailn1,sbrtn)
if(Ifailn1.eq.1)goto 122
goto 105
endif

do i7=i6+1,ii
perm(7)=bdsl(i7)
if(j.eq.7)then
call check3(j,perm,strno,totstr,nl,Ifailn1,sbrtn)
if(Ifailn1.eq.1)goto 122
goto 106
endif

106 enddo
105 enddo
104 enddo
103 enddo
102 enddo
101 enddo
    enddo
 enddo

122 Ifailn=Ifailn1



121 print*,'Ifailn',Ifailn
print*,'ssssssssssssssssssssssssssssssssssssssssssssssssssssss'
!if(Ifailn.eq.1)stop
return
end subroutine check_ind
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
