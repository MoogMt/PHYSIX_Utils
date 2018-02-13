program hydrogenion
implicit none
integer, parameter :: np = 6000 ! number of radial points
integer, parameter :: nx = 200 ! number of x points
integer, parameter :: ny = 100 ! number of y points
integer, parameter :: nz = 100 ! number of z points
real(8), parameter :: a0 = 1._8 ! Borh's radius
real(8), dimension(:,:), allocatable :: rad
real(8), dimension(:,:), allocatable :: rad2
real(8), dimension(:), allocatable :: c
real(8), dimension(:,:,:,:), allocatable :: psi
integer :: i, j, l, m, n, fact, ix, iy, iz, imol, nmol, idens
real(8) :: r, dr, norm, theta, phi, rmol, x, y, z, rleft, rright, dx
real(8) :: laguerre, s_harm, psileft, psiright
real(8) :: rho(2), vol, vkin(2), vpot(2), pot, etot(2)
real(8) :: rmin, rmax, deltar, rdens
namelist /input/ rmin, rmax, deltar, dx, rdens
rmin = 0.5_8
rmax = 10.0_8
deltar = 1._8/2
dx = 1._8/20
rdens = 2.5_8

open(10,file='mol_2px.input') ; read(10,nml=input) ; close(10)
open(8,file='energies_mol_2px.out')
open(7,file='density_3D_mol_2px.out')

if(rmin.lt.deltar/2.or.rmin.gt.rmax) then
  write(*,*) "ERROR !! bad choice for rmin"
  stop
endif

if(rmax.gt.nx*dx) then
  write(*,*) "ERROR !! bad choice for rmax"
  stop
endif

if(rmax.gt.nx*dx*0.75) then
  write(*,*) "WARNING !! rmax = ",rmax," is probably too large "
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

if(deltar.lt.0) then
  write(*,*) "ERROR !! deltar cannot be negative"
  stop
endif
if(deltar.gt.0.5) then
  write(*,*) "WARNING, deltar = ", deltar, " is quite large !"
endif

n = 2
allocate(rad(0:n-1,0:np),c(0:n-1))
rad = 0._8
do l = 0, n-1
  c(l) =  dble(fact(n-l-1))/(2*n*fact(n+l)**3)*(2/(n*a0))**(l+1.5_8)
! write(6,*) c(l)
! c(l) =  (dble(fact(n-l-1))/(2*n*fact(n+l)**3))**(0.5_8)*(2/(n*a0))**(l+1.5_8)
! write(6,*) c(l)
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

deallocate(c)

n = 2
allocate(rad2(0:n-1,0:np),c(0:n-1))
rad2 = 0._8

do l = 0, n-1
  c(l) =  dble(fact(n-l-1))/(2*n*fact(n+l)**3)*(2/(n*a0))**(l+1.5_8)
! c(l) =  (dble(fact(n-l-1))/(2*n*fact(n+l)**3))**(0.5_8)*(2/(n*a0))**(l+1.5_8)
! write(6,*) c(l)
enddo


dr = 30._8/np
do i = 0, np
  r = i*dr
  do l = 0, n-1
    rad2(l,i) =  c(l)*r**l*exp(-r/(n*a0))*laguerre(n-l-1,2*l+1,2*r/(n*a0))
  enddo
enddo
do l = 0, n-1
  norm = 0._8
  do i = 0, np
    norm = norm + (i*dr*rad2(l,i))**2*dr
  enddo
  rad2(l,:) = rad2(l,:)/sqrt(norm)
enddo

allocate(psi(2,-nx:nx,-ny:ny,-nz:nz))
psi = 0._8

nmol=(rmax-rmin)/deltar
idens=nint((rdens-rmin)/deltar)

do imol = 0,nmol

rmol = rmin + imol*deltar
!write(6,*) imol, rmol

rho = 0._8
vkin = 0._8
vpot = 0._8
vol = 0._8
dx = 1._8/20

do ix = -nx, nx-1
  x = (0.5_8 + ix) * dx
  do iy = -ny, ny-1
    y = (0.5_8 + iy) * dx
    do iz = -nz, nz-1
      z = (0.5_8 + iz) * dx
      psileft = 0._8
      psiright = 0._8

      rleft = sqrt((x+rmol/2)**2 + y**2 + z**2)
      theta = acos(z/rleft)
      if(abs(x+rmol/2).gt.1.e-3) then
        phi = atan2(y,(x+rmol/2))
      else
        phi = 0.
      endif
      i = nint(rleft/dr)
      if(i.le.np) then
        psileft = rad(1,i)*s_harm(1,1,theta,phi)
      endif

      rright = sqrt((x-rmol/2)**2 + y**2 + z**2)
      theta = acos(z/rright)
      if(abs(x-rmol/2).gt.1.e-3) then
        phi = atan2(y,(x-rmol/2))
      else
        phi = 0.
      endif
      i = nint(rright/dr)
      if(i.le.np) then
        psiright = rad2(1,i)*s_harm(1,1,theta,phi)
      endif
!     psi(1,ix,iy,iz)=psileft
!     psi(2,ix,iy,iz)=psiright
      psi(1,ix,iy,iz)=(psileft+psiright)/sqrt(2._8)
      psi(2,ix,iy,iz)=(psileft-psiright)/sqrt(2._8)
      rho(1) = rho(1) + (psi(1,ix,iy,iz)**2) * (dx**3)
      rho(2) = rho(2) + (psi(2,ix,iy,iz)**2) * (dx**3)
      vol = vol + (dx**3)
      vpot(1) = vpot(1) + (dx**3) * psi(1,ix,iy,iz)**2 * (pot(rleft,dx)+pot(rright,dx))
      vpot(2) = vpot(2) + (dx**3) * psi(2,ix,iy,iz)**2 * (pot(rleft,dx)+pot(rright,dx))
    enddo
 if(imol.eq.idens) write(7,'(2f8.4,2f12.6)') x,y,psi(1,ix,iy,0)**2,psi(2,ix,iy,0)**2
  enddo
 if(imol.eq.idens)  write(7,*)
enddo

do ix = -nx+1, nx-2
  x = (0.5_8 + ix) * dx
  do iy = -ny+1, ny-2
    y = (0.5_8 + iy) * dx
    do iz = -nz+1, nz-2
      z = (0.5_8 + iz) * dx
      i = nint(rleft/dr)
        vkin(1) = vkin(1) -1._8*dx*psi(1,ix,iy,iz)*                         &
                (psi(1,ix-1,iy,iz)+psi(1,ix+1,iy,iz)-2*psi(1,ix,iy,iz) &
                +psi(1,ix,iy-1,iz)+psi(1,ix,iy+1,iz)-2*psi(1,ix,iy,iz) &
                +psi(1,ix,iy,iz-1)+psi(1,ix,iy,iz+1)-2*psi(1,ix,iy,iz))
        vkin(2) = vkin(2) -1._8*dx*psi(2,ix,iy,iz)*                         &
                (psi(2,ix-1,iy,iz)+psi(2,ix+1,iy,iz)-2*psi(2,ix,iy,iz) &
                +psi(2,ix,iy-1,iz)+psi(2,ix,iy+1,iz)-2*psi(2,ix,iy,iz) &
                +psi(2,ix,iy,iz-1)+psi(2,ix,iy,iz+1)-2*psi(2,ix,iy,iz))
    enddo
  enddo
enddo

etot(1) = vpot(1)+vkin(1)
etot(2) = vpot(2)+vkin(2)

write(8,'(f6.3,14f10.4)') rmol, rho, vpot, vkin, etot/rho, etot/rho+2/rmol
write(6,'(f6.3,14f10.4)') rmol, rho, vpot, vkin, etot/rho, etot/rho+2/rmol

enddo
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

