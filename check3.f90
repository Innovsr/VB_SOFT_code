!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine check3(e,perm,str_sl,totstr,nl,ifailn1,sbrtn)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::i,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,e,a,nnstr,ifailn1,perm(20),wig2,nl,str_sl(1000),&
totbnd,totseorb,str_bnd(15000,1000),str_rad(500,20),bond_sl(10000,2),wig3,e2,e1,bnd(20),new_bd(100,2)
integer::totstr,bd_sl(20),num_orb(20)
character(9)::sbrtn

common /chek/totbnd,totseorb,str_bnd,str_rad,bond_sl,num_orb

print*,'***',sbrtn,'***'
do i=1,totstr
write(*,115),'check3',(str_bnd(str_sl(i),i1),i1=1,totbnd)
enddo
print*,'perm',(perm(i),i=1,e)
print*,'sl_num',(str_sl(i),i=1,totstr)
115 format(a,x,50I3)
ifailn1=0
i2=0
do i=1,e
do i1=1,2
i2=i2+1
bnd(i2)=bond_sl(perm(i),i1)
enddo
enddo
print*,'bnd',(bnd(i),i=1,i2)

a=0
do i=1,i2
!new_bd(1)=bnd(i)
do i1=i+1,i2
a=a+1
new_bd(a,1)=bnd(i)
new_bd(a,2)=bnd(i1)
enddo
enddo

do i=1,a
print*,'new_bd',(new_bd(i,i1),i1=1,2)
enddo

do i=1,a
do i3=1,totbnd
i6=0
do i4=1,2
do i5=1,2
if(new_bd(i,i4).eq.bond_sl(i3,i5))then
i6=i6+1
endif
enddo
enddo
if(i6.eq.2)then
bd_sl(i)=i3
goto 120
endif
enddo
120 enddo
print*,'bd_sl',(bd_sl(i),i=1,a)


i8=0
do i=1,totstr
i7=0
do i5=1,a
if(str_bnd(str_sl(i),bd_sl(i5)).eq.1)then
print*,'bd_sl(i5)',bd_sl(i5),i,i5
i7=i7+1
endif
enddo
if(i7.eq.e)then
print*,'i',i
i8=i8+1
endif
enddo

if(sbrtn.eq.'all_set'.or.sbrtn.eq.'check_2')then
e1=nae-nl*2-e*2-2
else
e1=nae-nl*2-e*2
endif
e2=e*2
print*,'e1,e2',e1,e2
call wigner(e1,wig2)
call wigner(e2,wig3)
if(i8.gt.wig3*wig2)ifailn1=1
print*,'ifailn1',ifailn1,nl,nae,e,wig2,wig3,i8
!stop

print*,'exit check3'
return
end subroutine check3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
