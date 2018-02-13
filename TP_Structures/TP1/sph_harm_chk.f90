program check
implicit none
integer :: i, j, k, l, m
real, parameter :: pi = acos(-1.0_8)
real(8) :: s_harm, theta, phi, x, y, z
character (13) :: flnm = "sph_harm.0000"
print '("Entrer les valeurs de l et m : ",$)' ; read *, l, m
write(flnm(10:11),'(i2.2)') l ; write(flnm(12:13),'(i2.2)') m
print '("Fichier de résultats : ",a13)' , flnm

open(1,file=flnm)
do i = 1, 180
  theta = i*pi/180
  do j = 1, 360
    phi = j*pi/180
    write(1,'(3g13.5)') theta, phi, s_harm(l,m,theta,phi)
  enddo
enddo
end program check


real(8) function s_harm(l,m,theta,phi)
implicit none
integer :: l, m, fact
real(8) :: theta, phi, assoc_legendre
real(8), parameter :: pi = acos(-1._8)
if ( m==0 ) then
  s_harm = sqrt((2*l+1)/(4*pi))*assoc_legendre(l,0,cos(theta))
else
  s_harm = sqrt(2*(2*l+1)*fact(l-m)/(4*pi*fact(l+m)))* &
              assoc_legendre(l,m,cos(theta))*cos(m*phi)
endif
end function s_harm

include 'legendre.f90'

integer function fact(l)
implicit none
integer :: k, l, i
i = 1
if ( l > 1 ) then
 do k = 2, l
  i = i*k
 enddo
endif
fact = i
end function fact


