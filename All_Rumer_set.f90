!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine All_Rumer_set(rumstr,setno,nlonep)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! N! number of Runmer set generates here 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

common/orb1/orbs1,rstr,nlpr
integer::j,i,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,n1,k,nlpr
integer::setno,nlonep,rstr(500,20),permutation(100),orbs1(20),orbs2(20),rumstr(500,20)

nlpr=nlonep
print*,'enter All_Rumer_set',setno
do i=1,setno
do i1=1,nae
rstr(i,i1)=rumstr(i,i1)
enddo
enddo
i=0
do j=nlonep*2+1,nae
i=i+1
orbs2(i)=rumstr(1,j)
enddo
n1=i

!print*,'orns_rum',(orbs2(i),i=1,n1)

i=100
do j=1,n1
if(orbs2(j).lt.i)i=orbs2(j)
enddo

orbs1(1)=i
k=1
do j=1,20
i=i+1
!print*,i
do i1=1,n1
if(i.eq.orbs2(i1))then
k=k+1
orbs1(k)=i
print*,i
endif
enddo
enddo
!print*,(orbs1(i),i=1,n1)

j=0
do i1=1,n1
permutation(1)=orbs1(i1)

do i2=1,n1
if(permutation(1).eq.orbs1(i2))goto 201
permutation(2)=orbs1(i2)
if(n1.eq.2)then
j=j+1
call rumer(permutation,n1,j,setno,nlpr)
goto 201
endif

do i3=1,n1
do i=1,2
if(permutation(i).eq.orbs1(i3))goto 202
enddo
permutation(3)=orbs1(i3)
if(n1.eq.3)then
j=j+1
call rumer(permutation,n1,j,setno,nlpr)
goto 202
endif

do i4=1,n1
do i=1,3
if(permutation(i).eq.orbs1(i4))goto 203
enddo
permutation(4)=orbs1(i4)
if(n1.eq.4)then
j=j+1
call rumer(permutation,n1,j,setno,nlpr)
goto 203
endif

do i5=1,n1
do i=1,4
if(permutation(i).eq.orbs1(i5))goto 204
enddo
permutation(5)=orbs1(i5)
if(n1.eq.5)then
j=j+1
call rumer(permutation,n1,j,setno,nlpr)
goto 204
endif

do i6=1,n1
do i=1,5
if(permutation(i).eq.orbs1(i6))goto 205
enddo
permutation(6)=orbs1(i6)
if(n1.eq.6)then
j=j+1
call rumer(permutation,n1,j,setno,nlpr)
goto 205
endif
do i7=1,n1
do i=1,6
if(permutation(i).eq.orbs1(i7))goto 206
enddo
permutation(7)=orbs1(i7)
if(n1.eq.7)then
j=j+1
call rumer(permutation,n1,j,setno,nlpr)
goto 206
endif

do i8=1,n1
do i=1,7
if(permutation(i).eq.orbs1(i8))goto 207
enddo
permutation(8)=orbs1(i8)
if(n1.eq.8)then
j=j+1
call rumer(permutation,n1,j,setno,nlpr)
goto 207
endif

do i9=1,n1
do i=1,8
if(permutation(i).eq.orbs1(i9))goto 208
enddo
permutation(9)=orbs1(i9)
if(n1.eq.9)then
j=j+1
call rumer(permutation,n1,j,setno,nlpr)
goto 208
endif

do i10=1,n1
do i=1,9
if(permutation(i).eq.orbs1(i10))goto 209
enddo
permutation(10)=orbs1(i10)
if(n1.eq.10)then
j=j+1
call rumer(permutation,n1,j,setno,nlpr)
goto 209
endif

do i11=1,n1
do i=1,10
if(permutation(i).eq.orbs1(i11))goto 210
enddo
permutation(11)=orbs1(i11)
if(n1.eq.11)then
j=j+1
call rumer(permutation,n1,j,setno,nlpr)
goto 210
endif
210 enddo
209 enddo
208 enddo
207 enddo
206 enddo
205 enddo
204 enddo
203 enddo
202 enddo
201 enddo
200 enddo
print*,'exit All_Rumer_set'

return
end subroutine All_Rumer_set
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
