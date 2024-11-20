!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine eq_dstr_set(tns,lnp,npstr,str2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::i,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,&
l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16,l17,l18,l19,l20,j
integer::tns,lnp,npstr,fail,strsln(1000),str2(15000,20),num_str,numpstr
common/nst/num_str,numpstr

numpstr=npstr
num_str=0

print*,'enter eq_dstr_set'
!write(9,*),'tns,lnp,npstr',tns,lnp,npstr
!call check_str_bond(lnp,tns,str2,npstr)

do i=1,tns
write(*,200)i,(str2(i,i1),i1=1,nae)
enddo
200 format (30I5)

!stop
do i1=1,tns
do i=1,200
strsln(i)=0
enddo
j=0
j=j+1
strsln(j)=i1

l1=j
do i2=i1+1,tns
do i=l1+1,j
strsln(i)=0
enddo
j=l1
j=j+1
strsln(j)=i2
call eq_dst_check(j,strsln,fail,lnp,str2)
if (fail.eq.1) goto 102

l2=j

do i3=i2+1,tns
do i=l2+1,j
strsln(i)=0
enddo
j=l2
j=j+1
strsln(j)=i3
call eq_dst_check(j,strsln,fail,lnp,str2)
if (fail.eq.1) goto 103

l3=j


do i4=i3+1,tns
do i=l3+1,j
strsln(i)=0
enddo
j=l3
j=j+1
strsln(j)=i4
call eq_dst_check(j,strsln,fail,lnp,str2)
if (fail.eq.1) goto 104

l4=j


do i5=i4+1,tns
do i=l4+1,j
strsln(i)=0
enddo
j=l4
j=j+1
strsln(j)=i5 
call eq_dst_check(j,strsln,fail,lnp,str2)
if (fail.eq.1) goto 105

l5=j


do i6=i5+1,tns
do i=l5+1,j
strsln(i)=0
enddo
j=l5
j=j+1
strsln(j)=i6
call eq_dst_check(j,strsln,fail,lnp,str2)
if (fail.eq.1) goto 106

l6=j


do i7=i6+1,tns
do i=l6+1,j
strsln(i)=0
enddo
j=l6
j=j+1
strsln(j)=i7
call eq_dst_check(j,strsln,fail,lnp,str2)
if (fail.eq.1) goto 107

l7=j

do i8=i7+1,tns
do i=l7+1,j
strsln(i)=0
enddo
j=l7
j=j+1
strsln(j)=i8
call eq_dst_check(j,strsln,fail,lnp,str2)
if (fail.eq.1) goto 108

l8=j

do i9=i8+1,tns
do i=l8+1,j
strsln(i)=0
enddo
j=l8
j=j+1
strsln(j)=i9
call eq_dst_check(j,strsln,fail,lnp,str2)
if (fail.eq.1) goto 109

l9=j

do i10=i9+1,tns
do i=l9+1,j
strsln(i)=0
enddo
j=l9
j=j+1
strsln(j)=i10
call eq_dst_check(j,strsln,fail,lnp,str2)
if (fail.eq.1) goto 110

l10=j

do i11=i10+1,tns
do i=l10+1,j
strsln(i)=0
enddo
j=l10
j=j+1
strsln(j)=i11
call eq_dst_check(j,strsln,fail,lnp,str2)
if (fail.eq.1) goto 111

l11=j

do i12=i11+1,tns
do i=l11+1,j
strsln(i)=0
enddo
j=l11
j=j+1
strsln(j)=i12
call eq_dst_check(j,strsln,fail,lnp,str2)
if (fail.eq.1) goto 112

l12=j

do i13=i12+1,tns
do i=l12+1,j
strsln(i)=0
enddo
j=l12
j=j+1
strsln(j)=i13
call eq_dst_check(j,strsln,fail,lnp,str2)
if (fail.eq.1) goto 113

l13=j

do i14=i13+1,tns
do i=l13+1,j
strsln(i)=0
enddo
j=l13
j=j+1
strsln(j)=i14
call eq_dst_check(j,strsln,fail,lnp,str2)
if (fail.eq.1) goto 114

l14=j

do i15=i14+1,tns
do i=l14+1,j
strsln(i)=0
enddo
j=l14
j=j+1
strsln(j)=i15
call eq_dst_check(j,strsln,fail,lnp,str2)
if (fail.eq.1) goto 115

l15=j

do i16=i15+1,tns
do i=l15+1,j
strsln(i)=0
enddo
j=l15
j=j+1
strsln(j)=i16
call eq_dst_check(j,strsln,fail,lnp,str2)
if (fail.eq.1) goto 116

l16=j

116 enddo
115 enddo
114 enddo
113 enddo
112 enddo
111 enddo
110 enddo
109 enddo
108 enddo
107 enddo
106 enddo
105 enddo
104 enddo
103 enddo
102 enddo
enddo


print*,'exit eq_dstr_set'
return
end subroutine eq_dstr_set
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
