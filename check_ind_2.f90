!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine check_ind_2(permut,nnstr,nl,n,Ifailn)

use commondat
implicit none

integer::i,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,e,a,nnstr,Ifailn,permut(100),num_orb(20),wig2,nl
integer::str_bnd_2(1000),totbnd,totseorb,str_bnd(15000,1000),str_rad(500,20),bond_sl(10000,2)
integer::bdsl(20),perm(20),n,ii,e1,e2,wig3,j,Ifailn1
character(9)::sbrtn

common /chek/totbnd,totseorb,str_bnd,str_rad,bond_sl,num_orb

print*,'enter check_ind_2'
print*,'permut',(permut(i),i=1,nnstr)
Ifailn=0

do i2=1,totbnd
str_bnd_2(i2)=0
do i1=1,nnstr
str_bnd_2(i2)=str_bnd_2(i2)+str_bnd(permut(i1),i2)
enddo
enddo
!
!do i2=1,iii
!!write(*,100),'str',(str3(strno(ii),i),i=1,nae)
!write(*,100),'str_bnd_2',(str_bnd(str_sl(i2),i1),i1=1,totbnd)
!enddo
!write(*,100),'add_val_2',(str_bnd_2(i1),i1=1,totbnd)
!
i3=0
do i2=1,totbnd
if(str_bnd_2(i2).eq.nnstr)then
i3=i3+1
endif
enddo

e=nae-nl*2-i3*2
call wigner(e,wig2)
!print*,'wig2ll',wig2,i3
if(nnstr.gt.wig2)then
Ifailn=1
endif
print*,'Ifailn',Ifailn
if(Ifailn.eq.1)goto 100

!ii=0
!do i=1,totbnd
!if(i.eq.n)goto 400
!if(str_bnd(permut(nnstr),i).eq.1)then
!ii=ii+1
!bdsl(ii)=i
!endif
!400 enddo
!print*,(bdsl(i),i=1,ii)
!!stop
!sbrtn='check_2'
!!if(ii)
!do j=2,int(ii/2)
!e1=nae-nl*2-j*2-2
!e2=j*2
!print*,'e1,e2',e1,e2
!call wigner2(e1,wig2)
!call wigner2(e2,wig3)
!
!print*,wig2*wig3,nnstr
!if(wig2*wig3.gt.nnstr)goto 300
!
!do i1=1,ii
!perm(1)=bdsl(i1)
!
!do i2=i1+1,ii
!perm(2)=bdsl(i2)
!if(j.eq.2)then
!!print*,(perm(i),i=1,j),j
!!stop
!call check3(j,perm,permut,nnstr,nl,Ifailn1,sbrtn)
!if(Ifailn1.eq.1)goto 122
!goto 201
!endif
!
!do i3=i2+1,ii
!perm(3)=bdsl(i3)
!if(j.eq.3)then
!call check3(j,perm,permut,nnstr,nl,Ifailn1,sbrtn)
!if(Ifailn1.eq.1)goto 122
!goto 202
!endif
!
!do i4=i3+1,ii
!perm(4)=bdsl(i4)
!if(j.eq.4)then
!call check3(j,perm,permut,nnstr,nl,Ifailn1,sbrtn)
!if(Ifailn1.eq.1)goto 122
!goto 203
!endif
!
!do i5=i4+1,ii
!perm(5)=bdsl(i5)
!if(j.eq.5)then
!call check3(j,perm,permut,nnstr,nl,Ifailn1,sbrtn)
!if(Ifailn1.eq.1)goto 122
!goto 204
!endif
!
!do i6=i5+1,ii
!perm(6)=bdsl(i6)
!if(j.eq.6)then
!call check3(j,perm,permut,nnstr,nl,Ifailn1,sbrtn)
!if(Ifailn1.eq.1)goto 122
!goto 205
!endif
!
!do i7=i6+1,ii
!perm(7)=bdsl(i7)
!if(j.eq.7)then
!call check3(j,perm,permut,nnstr,nl,Ifailn1,sbrtn)
!if(Ifailn1.eq.1)goto 122
!goto 206
!endif
!
!206 enddo
!205 enddo
!204 enddo
!203 enddo
!202 enddo
!201 enddo
!    enddo
!300 enddo
!
!122 Ifailn=Ifailn1

!print*,'Ifailn',Ifailn,wig2,i3
100 print*,'exit check_ind_2'
return
end subroutine check_ind_2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
