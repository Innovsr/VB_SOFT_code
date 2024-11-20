!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine charge_cntr(str2,nl,strnum)
use commondat
implicit none

integer::i,i1,i2,i3,i4,iii,j,j1,nl,str2(15000,20),strnum,nn(10)
real*8::cent_coord(5000,3),coordx(100),coordy(100),coordz(100),chrg_cntrx,chrg_cntry,chrg_cntrz,charge_p_cntr(50),&
chrg_cntr(1000,3)

common/coordinate/coordx,coordy,coordz

!i4=0

do i=1,atom
print*,'coord',coordx(i),coordy(i),coordz(i)
enddo
do i=1,strnum
j1=0
print*,'str2',(str2(i,i1),i1=1,nae)
do i1=1+nl*2,(nae-nlast),2
iii=0
nn(1)=0
nn(2)=0
do i2=i1,i1+1
do i3=1,atom
do i4=1,atn(active_atoms(i3))
if(str2(i,i2).eq.atoset(active_atoms(i3),i4))then
iii=iii+1
if(mod(i2,2).eq.0)nn(2)=active_atoms(i3)
if(mod(i2,2).eq.1)nn(1)=active_atoms(i3)
goto 101
endif
enddo
enddo
101 enddo
!if(nn(2).eq.nn(1))goto 202
if(iii.ne.2) goto 202
print*,'nn(1),nn(2)',i,nn(1),nn(2)

j1=j1+1
cent_coord(j1,1)=(coordx(nn(1))+coordx(nn(2)))/2.0
cent_coord(j1,2)=(coordy(nn(1))+coordy(nn(2)))/2.0
cent_coord(j1,3)=(coordz(nn(1))+coordz(nn(2)))/2.0
charge_p_cntr(j1)=2.0

202 enddo
!if (nlast.ne.0)then
!do i1=nae-nlast+1,nae
!do i3=1,atom
!do i4=1,atn(active_atoms(i3))
!if(str2(i,i1).eq.atoset(active_atoms(i3),i4))then
!nn(1)=active_atoms(i3)
!print*,'radical_atom',nn(1)
!j1=j1+1
!cent_coord(j1,1)=coordx(nn(1))
!cent_coord(j1,2)=coordy(nn(1))
!cent_coord(j1,3)=coordz(nn(1))
!charge_p_cntr(j1)=1.0
!endif
!enddo
!enddo
!enddo
!endif

do i1=1,j1
print*,'cent_coord',(cent_coord(i1,i2),i2=1,3)
enddo
chrg_cntrx=0
chrg_cntry=0
chrg_cntrz=0
do i1=1,j1
chrg_cntrx=chrg_cntrx+charge_p_cntr(i1)*cent_coord(i1,1)
chrg_cntry=chrg_cntry+charge_p_cntr(i1)*cent_coord(i1,2)
chrg_cntrz=chrg_cntrz+charge_p_cntr(i1)*cent_coord(i1,3)
enddo
chrg_cntr(i,1)=chrg_cntrx/(nae-nl*2)
chrg_cntr(i,2)=chrg_cntry/(nae-nl*2)
chrg_cntr(i,3)=chrg_cntrz/(nae-nl*2)
print*,'chrg_cntr',(chrg_cntr(i,j),j=1,3)
enddo


!stop

return
end subroutine charge_cntr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
