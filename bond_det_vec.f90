subroutine bond_det_vec(nl,str_num,str3,bnd_std)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::i,i1,i2,i3,l,ll,j,nl,totbnd,totseorb,str_bnd(15000,1000),str_rad(500,20)&
,bond_sl(10000,2),num_orb(20),strno(15000),totstr,str3(2000,20),str_num,bnd_str(1000),&
avg
Double Precision::bond_vec_mat(150,150),Eig(150),bond_det
real::bnd_std,var
common /chek/totbnd,totseorb,str_bnd,str_rad,bond_sl,num_orb
common /bnd_mat/bond_vec_mat

!print*,'totstr',str_num
!do i=1,str_num
!write(9,111),'str_bnd',(str3(i,j))
!enddo
111 format (a,2x,30I3)

do i=1,totbnd
!print*,'bond_sl',(bond_sl(i,i1),i1=1,2)
ll=0
do j=1,str_num
!write(*,111),'str',(str3(j,i1),i1=1,nae)
do i2=nl*2+1,nae-nlast,2
l=0
do i1=1,2
do i3=i2,i2+1
!print*,'str3(j,i3).eq.bond_sl(i,i1)',str3(j,i3),bond_sl(i,i1)
if(str3(j,i3).eq.bond_sl(i,i1))then
l=l+1
!print*,l
goto 300
endif
enddo
300 enddo
!print*,'l',l
if(l.eq.2)then
ll=ll+1
!print*,'ll*',ll
endif
enddo
enddo
!print*,'ll',ll
bnd_str(i)=ll
enddo

avg=0
do i=1,totbnd
avg=avg+bnd_str(i)
enddo
avg=avg/totbnd

var=0
do i=1,totbnd
var=var+(bnd_str(i)-avg)**2
enddo

bnd_std=sqrt(var/totbnd)

!write(9,*)(bnd_str(i),i=1,totbnd)
!write(*,*),'bnd_str',bnd_std

!stop

!do i=1,str_num
!write(*,111),'str',(str3(i,j),j=1,nae)
!enddo
!do i=1,totbnd
!write(*,111),'bond_sl',(bond_sl(i,i1),i1=1,2)
!enddo
!
!stop

!do i=1,150
!do i1=1,150
!bond_vec_mat(i,i1)=0.0
!enddo
!enddo

!do i=1,totstr
!do i2=1,totbnd
!str_bnd(strno(i),i2)=str_bnd(strno(i),i2)+1
!enddo
!enddo

!do i=1,totstr
!do i1=i,totstr
!do i2=1,totbnd
!bond_vec_mat(i,i1)=bond_vec_mat(i,i1)+str_bnd(strno(i),i2)*str_bnd(strno(i1),i2)
!
!enddo
!enddo
!enddo
!
!
!do i=1,totstr
!do i1=i,totstr
!bond_vec_mat(i1,i)=bond_vec_mat(i,i1)
!
!enddo
!enddo
!
!!do i=1,totstr
!!write(9,111),'bond_vec_mat',(int(bond_vec_mat(i,i1)),i1=1,totstr)
!!enddo 
!
!call MatLDR('bnv',strno,totstr,Eig)
!
!!write(9,*)(Eig(i),i=1,totstr)
!
!bond_det=1.0
!do i=1,totstr
!bond_det=bond_det*Eig(i)
!enddo
!write(9,*),'bond_det',bond_det

!stop
return
end subroutine bond_det_vec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
