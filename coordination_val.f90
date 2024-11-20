!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine coordination_val(coordination_mat,numbond,loop_score_row,bndscore,coord_score)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::i,i1,i2,i3
integer::coordination_mat(20,20),numbond,loop_score_row(20),new_row,loop_score_new(20)
real*8::bndscore(20),coord_score,c_score,score
common/lp/loop_score_new

print*,'enter coordination_val ********************************'
do i=1,numbond
print*,'coordination_mat',(coordination_mat(i,i2),i2=1,loop_score_row(i))
enddo

do i=1,20
loop_score_new(i)=loop_score_row(i)
enddo

coord_score=0.0
do i=1,numbond
!print*,'coordination_mat,i',i,'>',(coordination_mat(i,i2),i2=1,loop_score_row(i))
if (coordination_mat(i,1).ne.0)then
c_score=0.0
do i1=i+1,numbond
!print*,'coordination_mat,i1',i1,'>',(coordination_mat(i1,i2),i2=1,loop_score_row(i1))
if(i.ne.i1)then
if (coordination_mat(i1,1).ne.0)then
score=0.0
do i2=1,loop_score_row(i)
do i3=1,loop_score_row(i1)
print*,'i2,i3',i2,i3
if (coordination_mat(i,i2).eq.coordination_mat(i1,i3))then
!score=score+(bndscore(i)/(loop_score_row(i)-1))*((loop_score_row(i)-1)-(new_row(i2,i)-1))+ &
!(bndscore(i1)/(loop_score_row(i1)-1))*((loop_score_row(i1)-1)-(new_row(i3,i1)-1))
score=score+(bndscore(i)/new_row(i2,i)+bndscore(i1)/new_row(i3,i1))
print*,'bndscore(i),bndscore(i1),(i,i2),(i1,i3),new_row(i2,i),new_row(i3,i1)',&
bndscore(i),bndscore(i1),'(',i,i2,')','(',i1,i3,')',new_row(i2,i),new_row(i3,i1)
!print*,'i,i1,i2,i3,score',i,i1,i2,i3,score
endif
enddo
enddo
if(score.ne.0.0)c_score=c_score+1.0/score
!c_score=c_score+score
!print*,'c_score',c_score
endif
endif
enddo
!if(c_score.ne.0.0)coord_score=coord_score+1.0/c_score
if(c_score.ne.0.0)coord_score=coord_score+c_score
print*,'coord_score',coord_score
print*,'*************'
endif
enddo


return
end subroutine coordination_val
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
