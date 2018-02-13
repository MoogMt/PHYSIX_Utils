program hydrogenion
implicit none
integer, parameter :: np = 6000 ! number of radial points
integer, parameter :: nx = 200 ! number of x points
integer, parameter :: ny = 200 ! number of y points
integer, parameter :: nz = 200 ! number of z points
real(8), parameter :: a0 = 1._8 ! Borh's radius
real(8), dimension(:,:), allocatable :: rad
real(8), dimension(:,:), allocatable :: rad2
real(8), dimension(:), allocatable :: c
real(8), dimension(:,:,:), allocatable :: psi
integer :: i, j, l, m, n, fact, ix, iy, iz, ll, mm
real(8) :: r, dr, norm, theta, phi, x, y, z, ratom, dx
real(8) :: laguerre, s_harm, psiatom
real(8) :: rho, vol, vkin, vpot, pot, etot
character(9) :: flnm="hydro.0_0"
namelist /input/ n,l,m,dx

n = 1
l = 0
m = 0
dx = 1._8/20

open(10,file='atom.input') ; read(10,nml=input) ; close(10)
open(17,file='radial_functions.out') 
open(18,file='energies.out') 
open(19,file='density_3D.out') 
if(n.lt.1.or.n.gt.3) then
  write(*,*) "ERROR !! n cannot be smaller than 1 or larger than 3"
  stop
endif
if(l.ge.n) then
  write(*,*) "ERROR !! l cannot be greater than or equal to n"
  stop
endif
if(abs(m).gt.abs(l)) then
  write(*,*) "ERROR !! m is bound by l"
  stop
endif
if(dx.lt.0) then
  write(*,*) "ERROR !! dx cannot be negative"
  stop
endif
if(dx.gt.0.5) then
  write(*,*) "WARNING, dx = ", dx, " is quite large !"
endif
if(dx*nx.gt.20) then
  write(*,*) "WARNING, x = ", dx*nx, " is quite large !"
endif

ll = l
mm = m

allocate(rad(0:n-1,0:np),c(0:n-1))
rad = 0._8
do l = 0, n-1
  c(l) =  dble(fact(n-l-1))/(2*n*fact(n+l)**3)*(2/(n*a0))**(l+1.5_8)
enddo
dr = 30._8/np
do i = 0, np
  r = i*dr
  do l = 0, n-1
    rad(l,i) =  c(l)*r**l*exp(-r/(n*a0))*laguerre(n-l-1,2*l+1,2*r/(n*a0))
  enddo
enddo
do l = 0, n-1
  norm = 0._8
  do i = 0, np
    norm = norm + (i*dr*rad(l,i))**2*dr
  enddo
  rad(l,:) = rad(l,:)/sqrt(norm)
enddo

do i = 0,np
  r = i*dr
  write(17,*) r,(rad(l,i),l=0,n-1)
enddo

allocate(psi(-nx:nx,-ny:ny,-nz:nz))

psi = 0._8

rho = 0._8
vkin = 0._8
vpot = 0._8
vol = 0._8


l = ll
m = mm

do ix = -nx, nx-1
  x = (0.5_8 + ix) * dx
  do iy = -ny, ny-1
    y = (0.5_8 + iy) * dx
    do iz = -nz, nz-1
      z = (0.5_8 + iz) * dx
      psiatom = 0._8
      ratom = sqrt(x**2 + y**2 + z**2)
      theta = acos(z/ratom)
      if(abs(x).gt.1.e-3) then
        phi = atan2(y,x)
      else
        phi = 0.
      endif
      i = nint(ratom/dr)
      if(i.le.np) then
        psiatom = rad(l,i)*s_harm(l,m,theta,phi)
!       psiatom = rad(0,i)*s_harm(0,0,0._8,0._8)
      endif
      psi(ix,iy,iz)=psiatom
      rho = rho + (psi(ix,iy,iz)**2) * (dx**3)
      vol = vol + (dx**3)
      vpot = vpot + (dx**3) * psi(ix,iy,iz)**2 * pot(ratom,dx)
    enddo
    write(19,'(2f8.4,2f12.6)') x,y,psi(ix,iy,0),psi(ix,iy,0)**2
  enddo
  write(19,*)
enddo

do ix = -nx+1, nx-2
  x = (0.5_8 + ix) * dx
  do iy = -ny+1, ny-2
    y = (0.5_8 + iy) * dx
    do iz = -nz+1, nz-2
      z = (0.5_8 + iz) * dx
        vkin = vkin -1._8*dx*psi(ix,iy,iz)*                         &
                (psi(ix-1,iy,iz)+psi(ix+1,iy,iz)-2*psi(ix,iy,iz) &
                +psi(ix,iy-1,iz)+psi(ix,iy+1,iz)-2*psi(ix,iy,iz) &
                +psi(ix,iy,iz-1)+psi(ix,iy,iz+1)-2*psi(ix,iy,iz))
    enddo
  enddo
enddo

etot = vpot+vkin

write(18,'(f6.3,3f12.4)') rho, vpot/rho, vkin/rho, etot/rho
write(*,'(''Linear size of the integration volume (bohr) = '',f12.4)') dx*nx*2 
write(*,'(''Wavefunction normalization      = '',f12.4)')  rho 
write(*,'(''Av. value of pot Energy (Ryd)   = '',f12.4)')  vpot/rho
write(*,'(''Av. value of kin Energy (Ryd)   = '',f12.4)')  vkin/rho
write(*,'(''Av. value of total Energy (Ryd) = '',f12.4)')  etot/rho

stop

end program hydrogenion

include 'laguerre.f90'
include 'sph_harm.f90'

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

real(8) function pot(u,dr)
implicit none
real(8) :: u, dr
if(abs(u).le.dr) pot = - 2._8/(2._8 * dr)
if(abs(u).gt.dr) pot = - 2._8/u
end function pot

