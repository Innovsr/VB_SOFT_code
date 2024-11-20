!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine all_set_gen(str_sl,nnstr,nl,n,Ifailn)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::i,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,ii,elporb,nl,a,nnstr,Ifailn,permut(100),&
str_sl(1000),wig2,l,k,bdsl(20),perm(20),j,n
integer::e1,e2,wig3,Ifailn1
character(9)::sbrtn
integer::totbnd,totseorb,str_bnd(15000,1000),str_rad(500,20),bond_sl(10000,2),num_orb(20)

common /chek/totbnd,totseorb,str_bnd,str_rad,bond_sl,num_orb
print*,'enter all_set_gen'
print*,'str_sl',(str_sl(i),i=1,nnstr)



ii=0
do i=1,totbnd
if(i.eq.n)goto 400
if(str_bnd(str_sl(nnstr),i).eq.1)then
ii=ii+1
bdsl(ii)=i
endif
400 enddo
print*,(bdsl(i),i=1,ii)
!stop
sbrtn='all_set'
!if(ii)
do j=2,int(ii/2)
e1=nae-nl*2-j*2-2
e2=j*2
print*,'e1,e2',e1,e2
call wigner(e1,wig2)
call wigner(e2,wig3)

print*,wig2*wig3,nnstr
if(wig2*wig3.gt.nnstr)goto 300

do i1=1,ii
perm(1)=bdsl(i1)

do i2=i1+1,ii
perm(2)=bdsl(i2)
if(j.eq.2)then
!print*,(perm(i),i=1,j),j
!stop
call check3(j,perm,str_sl,nnstr,nl,Ifailn1,sbrtn)
if(Ifailn1.eq.1)goto 500
goto 201
endif

do i3=i2+1,ii
perm(3)=bdsl(i3)
if(j.eq.3)then
call check3(j,perm,str_sl,nnstr,nl,Ifailn1,sbrtn)
if(Ifailn1.eq.1)goto 500
goto 202
endif

do i4=i3+1,ii
perm(4)=bdsl(i4)
if(j.eq.4)then
call check3(j,perm,str_sl,nnstr,nl,Ifailn1,sbrtn)
if(Ifailn1.eq.1)goto 500
goto 203
endif

do i5=i4+1,ii
perm(5)=bdsl(i5)
if(j.eq.5)then
call check3(j,perm,str_sl,nnstr,nl,Ifailn1,sbrtn)
if(Ifailn1.eq.1)goto 500
goto 204
endif

do i6=i5+1,ii
perm(6)=bdsl(i6)
if(j.eq.6)then
call check3(j,perm,str_sl,nnstr,nl,Ifailn1,sbrtn)
if(Ifailn1.eq.1)goto 500
goto 205
endif

do i7=i6+1,ii
perm(7)=bdsl(i7)
if(j.eq.7)then
call check3(j,perm,str_sl,nnstr,nl,Ifailn1,sbrtn)
if(Ifailn1.eq.1)goto 500
goto 206
endif

206 enddo
205 enddo
204 enddo
203 enddo
202 enddo
201 enddo
    enddo
300 enddo









a=((nae-nl*2)-mod(nae-nl*2,2))/2
!print*,'aa',a,nae,mod(nae-nl*2,2)
    
do i=1,a-2
elporb=nae-nl*2-i*2
!print*,'elporb',elporb
call wigner(elporb,wig2)
!print*,'wig2',wig2,elporb,nnstr
if(wig2.lt.nnstr+1)then

!print*,'sourav'
do i1=1,nnstr-1
permut(1)=str_sl(i1)
if(wig2.eq.1)then
permut(2)=str_sl(nnstr)
l=2
call check_ind_2(permut,l,nl,n,Ifailn)
if(Ifailn.eq.1)goto 600
goto 100
endif

do i2=i1+1,nnstr-1
permut(2)=str_sl(i2)
if(wig2.eq.2)then
permut(3)=str_sl(nnstr)
l=3
!print*,'permut',(permut(k),k=1,3)
call check_ind_2(permut,l,nl,n,Ifailn)
if(Ifailn.eq.1)goto 600
goto 101
endif


do i3=i2+1,nnstr-1
permut(3)=str_sl(i3)
if(wig2.eq.3)then
permut(4)=str_sl(nnstr)
l=4
call check_ind_2(permut,l,nl,n,Ifailn)
if(Ifailn.eq.1)goto 600
goto 102
endif

do i4=i3+1,nnstr-1
permut(4)=str_sl(i4)
if(wig2.eq.4)then
permut(5)=str_sl(nnstr)
l=5
call check_ind_2(permut,l,nl,n,Ifailn)
if(Ifailn.eq.1)goto 600
goto 103
endif

do i5=i4+1,nnstr-1
permut(5)=str_sl(i5)
if(wig2.eq.5)then
permut(6)=str_sl(nnstr)
l=6
call check_ind_2(permut,l,nl,n,Ifailn)
if(Ifailn.eq.1)goto 600
goto 104
endif

do i6=i5+1,nnstr-1
permut(6)=str_sl(i6)
if(wig2.eq.6)then
permut(7)=str_sl(nnstr)
l=7
call check_ind_2(permut,l,nl,n,Ifailn)
if(Ifailn.eq.1)goto 600
goto 105
endif

do i7=i6+1,nnstr-1
permut(7)=str_sl(i7)
if(wig2.eq.7)then
permut(8)=str_sl(nnstr)
l=8
call check_ind_2(permut,l,nl,n,Ifailn)
if(Ifailn.eq.1)goto 600
goto 106
endif

do i8=i7+1,nnstr-1
permut(8)=str_sl(i8)
if(wig2.eq.8)then
permut(9)=str_sl(nnstr)
l=9
call check_ind_2(permut,l,nl,n,Ifailn)
if(Ifailn.eq.1)goto 600
goto 107
endif

do i9=i8+1,nnstr-1
permut(9)=str_sl(i9)
if(wig2.eq.9)then
permut(10)=str_sl(nnstr)
l=10
call check_ind_2(permut,l,nl,n,Ifailn)
if(Ifailn.eq.1)goto 600
goto 108
endif

do i10=i9+1,nnstr-1
permut(10)=str_sl(i10)
if(wig2.eq.10)then
permut(11)=str_sl(nnstr)
l=11
call check_ind_2(permut,l,nl,n,Ifailn)
if(Ifailn.eq.1)goto 600
goto 109
endif

do i11=i10+1,nnstr-1
permut(11)=str_sl(i11)
if(wig2.eq.11)then
permut(12)=str_sl(nnstr)
l=12
call check_ind_2(permut,l,nl,n,Ifailn)
if(Ifailn.eq.1)goto 600
goto 110
endif

110 enddo
109 enddo
108 enddo
107 enddo
106 enddo
105 enddo
104 enddo
103 enddo
102 enddo
101 enddo
100 enddo


endif
enddo
!500 stop
500 Ifailn=Ifailn1
print*,'exit all_set_gen'
600 return
end subroutine all_set_gen
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
