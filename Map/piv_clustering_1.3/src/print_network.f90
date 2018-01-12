subroutine print_network(nc,link)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!   Planar representation of a graph in a svg file                     !
!   using the Fruchterman-Reingold algorithm:                          !
!   SOFTWARE—PRACTICE AND EXPERIENCE, VOL. 21, 1129-1164 (1991)        !
!                                                                      !
!   In short, nodes repel each other with electrostatic force and      !
!   links between nodes are springs. A damped dynamics is performed,   !
!   reducing progressively the temperature.                            !
!   No walls are imposed, and the cooling function is a Fermi-Dirac.   !
!   Typically 300-5000 iterations give a relaxed system.               !
!                                                                      !
!   The input file contains pairs of nodes i j connected by a link.    !
!   Due to initial randomization, each run produces a different plot.  !
!   Nodes are colored from red to yellow to cyan with increasing       !
!   number of links.                                                   !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none
integer::nc,link(nc,nc)
!
integer :: niter
integer :: i,j,k
integer :: ic,jc,kc,nl,il,it,maxlink,minlink
integer,allocatable::nlink(:),ilink(:,:),nlinksort(:),bond(:,:)
real*8 :: kappa,kk,d,dt,temp0,temp,avec,avelink,ccoeff,cloc
real*8, allocatable :: posc(:,:),disp(:,:),dpos(:,:,:),ddpos(:,:)
real*8 :: minx,maxx,miny,maxy
integer :: np,x1,x2,y1,y2,radius,color(256,3),tmp,itmp,jmax
integer, allocatable :: isort(:)
!
write(*,*) "*** analyzing the network of clusters ***"
write(*,*) "(links are defined by Apollonius theorem applied to centers)"
!
allocate(posc(2,nc))
allocate(disp(2,nc))
allocate(nlink(nc))
allocate(dpos(2,nc,nc))
allocate(ddpos(nc,nc))
nl=sum(link(:,:))/2
allocate(bond(2,nl))
allocate(ilink(nc,nl))
nlink(1:nc)=0
il=0
do ic=1,nc-1
  do jc=ic+1,nc
    if(link(ic,jc).eq.1) then
      il=il+1
      bond(1,il)=ic
      bond(2,il)=jc
      nlink(ic)=nlink(ic)+1
      nlink(jc)=nlink(jc)+1
      ilink(ic,nlink(ic))=jc
      ilink(jc,nlink(jc))=jc
    endif
  enddo
enddo
avelink=sum(nlink(:))/dble(nc)
minlink=minval(nlink,1)
maxlink=maxval(nlink,1)
write(*,'(a,i10,f10.2,2i10)') ' links: tot aver min max =',nl,avelink,minlink,maxlink
ccoeff=0.d0
do ic=1,nc
  if(nlink(ic).eq.0)then
    write(*,*) 'ERROR: unexpectedly, this cluster has no neighbors:',ic,'  Exiting.'
    stop
  endif
  cloc=0.d0
  do j=1,nlink(ic)-1
    do k=j+1,nlink(ic)
      jc=ilink(ic,j)
      kc=ilink(ic,k)
      if(link(jc,kc).eq.1)cloc=cloc+1.d0
    enddo
  enddo
  if(cloc.gt.0.d0) then
    cloc=cloc/(0.5d0*dble(nlink(ic)*(nlink(ic)-1)))
    ccoeff=ccoeff+cloc
  endif
enddo
ccoeff=ccoeff/dble(nc)
write(*,'(a,f8.5)') ' global clustering coefficient of the network =',ccoeff
!
! ****** assign initial random positions to nodes
kk=1.d0/dble(nc*10)
kappa=sqrt(kk)
call init_random
do ic=1,nc
  call random_number(posc(1,ic))
  call random_number(posc(2,ic))
  posc(:,ic)=posc(:,ic)*kappa*10.d0+(1.d0-kappa*10.d0)*0.5
enddo
!
! ****** evolve the dynamics of the nodes
if(nc.le.500)then
  niter=10000
else
  niter=2000
endif
ddpos(:,:)=0.d0
temp0=kappa*0.1d0 
write(*,'(a,f10.7,a,f10.7,a,i6)') ' equil. dist. =',kappa,'  initial T =',temp0,'  niter =',niter
do it=1,niter
!  temp=temp0/(1.d0+dble((it-1)*100)/dble(niter))
  temp=temp0*1.d0/(1.d0+dexp((dble(it)/dble(niter)-0.6d0)*8.d0))
  ! distance matrix
  do ic=1,nc-1
    do jc=ic+1,nc
      dpos(1:2,ic,jc)=posc(1:2,ic)-posc(1:2,jc)
      dpos(1:2,jc,ic)=-dpos(1:2,ic,jc)
      ddpos(ic,jc)=sum(dpos(1:2,ic,jc)**2)
      ddpos(jc,ic)=ddpos(ic,jc)
    enddo
  enddo
  ! repulsive forces
  do ic=1,nc
    disp(1:2,ic)=0.d0
    do jc=1,nc
      if (ic.eq.jc) cycle
      if(ddpos(ic,jc).eq.0.d0) then
        ! random kick to avoid explosion
        call random_number(d)
        disp(1,ic)=disp(1,ic)+(d-0.5d0)*kappa*0.001d0
        call random_number(d)
        disp(2,ic)=disp(2,ic)+(d-0.5d0)*kappa*0.001d0
      else
        disp(1:2,ic)=disp(1:2,ic)+dpos(1:2,ic,jc)*kk/ddpos(ic,jc)
      endif
    enddo
  enddo
  ! attractive forces
  do il=1,nl
    ic=bond(1,il)
    jc=bond(2,il)
    d=ddpos(ic,jc)
    if (d.eq.0.d0) cycle
    d=sqrt(d)
    disp(1:2,ic)=disp(1:2,ic)-dpos(1:2,ic,jc)*d/kappa
    disp(1:2,jc)=disp(1:2,jc)+dpos(1:2,ic,jc)*d/kappa
  enddo
  ! evolve damped dynamics
  avec=0.d0
  do ic=1,nc
    d=sum(disp(1:2,ic)**2)
    if (d.eq.0.d0) cycle
    d=sqrt(d)
    dt=d
    if (dt>temp) dt=temp
    posc(1:2,ic)=posc(1:2,ic)+disp(1:2,ic)*dt/d
    avec=avec+dt
!    do i=1,2
!      if (posc(i,ic)<0.d0) posc(i,ic)=0.d0
!      if (posc(i,ic)>1.d0) posc(i,ic)=1.d0
!    enddo
  enddo 
  avec=avec/dble(nc)
  if (niter.ge.5) then
    if (mod(it,niter/5).eq.0) then
      write(*,'(a,i6,2f10.7)') '  step,T,K=',it,temp,avec
    endif
  else
    write(*,'(a,i6,2f10.7)') '  step,T,K=',it,temp,avec
  endif
enddo
!
! ****** prints output to a svg file
! 2400*2400 points in 20cm*20cm = 300 dpi
np=2400
radius=np/(nc*2)
if(radius>80) radius=80
if(radius<16) radius=16
call svg_open(30, "network.svg", 20,20, np,np)
minx=minval(posc(1,:),1)
maxx=maxval(posc(1,:),1)
miny=minval(posc(2,:),1)
maxy=maxval(posc(2,:),1)
do ic=1,nc
  posc(1,ic)=radius*2+(np-radius*4)*(posc(1,ic)-minx)/(maxx-minx)
  posc(2,ic)=radius*2+(np-radius*4)*(posc(2,ic)-miny)/(maxy-miny)
enddo
do il=1,nl
  ic=bond(1,il)
  jc=bond(2,il)
  x1=nint(posc(1,ic))
  y1=nint(posc(2,ic))
  x2=nint(posc(1,jc))
  y2=nint(posc(2,jc))
  call svg_line(30,x1,y1,x2,y2,4,0,0,0)
enddo
! color scale red-yellow-cyan
do i=1,128
  color(i,1)=255
  color(i,2)=i*2-1
  color(i,3)=0
  color(i+128,1)=256-i*2
  color(i+128,2)=255
  color(i+128,3)=i*2-1
enddo
! sort nodes based on number of bonds
allocate(nlinksort(nc),isort(nc))
nlinksort=nlink
do ic=1,nc
  isort(ic)=ic
enddo
jmax=nc-1
do ic=1,nc-1
  tmp=1000000
  do jc=1,jmax
    if (nlinksort(jc).lt.nlinksort(jc+1)) cycle
    tmp=nlinksort(jc)
    nlinksort(jc)=nlinksort(jc+1)
    nlinksort(jc+1)=tmp
    itmp=isort(jc)
    isort(jc)=isort(jc+1) 
    isort(jc+1)=itmp
  enddo
  if(tmp.eq.1000000)exit
enddo
do ic=1,nc
  jc=isort(ic)
  x1=nint(posc(1,jc))
  y1=nint(posc(2,jc))
  i=nint(dble(nlink(jc)-minlink+1)*256.d0/dble(maxlink-minlink+1))
  call svg_circle(30,x1,y1,radius,4,0,0,0, &
 & color(i,1),color(i,2),color(i,3))
enddo
call svg_close(30)
write(*,*) 'planar representation of network written in file network.svg'
!
! ****** clean memory
deallocate(nlinksort)
deallocate(isort)
deallocate(posc)
deallocate(disp)
deallocate(bond)
deallocate(ilink)
deallocate(nlink)
deallocate(dpos)
deallocate(ddpos)
!
end subroutine print_network
!
!===============================================================================
!           Routines to write simple svg vector graphics on file.
!-------------------------------------------------------------------------------
! usage:         open svg, add some shapes and close.
! compilation:   just compile together with your source.
!-------------------------------------------------------------------------------
!                                 Fabio Pietrucci, Aug 2007
!                                 fabio.pietrucci@gmail.com
!===============================================================================
subroutine svg_open(iunit,filename,dx_cm,dy_cm,dx_viewbox,dy_viewbox)
!
  implicit none
  character :: filename*11
  integer :: iunit,dx_cm,dy_cm,dx_viewbox,dy_viewbox
!
  open(iunit,file=filename,status='unknown')
  write(iunit,'(a)') '<?xml version="1.0" standalone="no"?>' 
  write(iunit,'(a)') '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"'
  write(iunit,'(a)') '  "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">'
  write(iunit,'(a,i4,a,i4,a,4i6,a)') &
   '<svg width="',dx_cm,'cm" height="',dy_cm,'cm" viewBox="', &
   0,0,dx_viewbox,dy_viewbox,'"'
  write(iunit,'(a)') '     xmlns="http://www.w3.org/2000/svg" version="1.1">'
  write(iunit,'(a)') '  <desc> generated by svg_writer.f90 </desc>'
!
end subroutine svg_open
!===============================================================================
subroutine svg_close(iunit)
! 
  implicit none
  integer :: iunit
!
  write(iunit,'(a)') '</svg>'
  close(iunit)
!
end subroutine svg_close
!===============================================================================
subroutine svg_line(iunit,x1,y1,x2,y2,stroke_width,r,g,b)
!
  implicit none
  integer :: iunit,x1,y1,x2,y2,stroke_width,r,g,b
!
  write(iunit,'(4(a,i4),a)') &
   '<line x1="',x1,'" y1="',y1,'" x2="',x2,'" y2="',y2,'" '
  write(iunit,'(a,i4,3(a,i3),a)') &
   ' stroke-width="',stroke_width,'" stroke="rgb(',r,',',g,',',b,')"/>'
!
end subroutine svg_line
!===============================================================================
subroutine svg_rect(iunit,x1,y1,dx,dy,stroke_width,r,g,b,r2,g2,b2)
!
  implicit none
  integer :: iunit,x1,y1,dx,dy,stroke_width,r,g,b,r2,g2,b2
!
  write(iunit,'(4(a,i4),a)') &
   '<rect x="',x1,'" y="',y1,'" width="',dx,'" height="',dy,'" '
  write(iunit,'(a,i4,3(a,i3),a)') &
   ' stroke-width="',stroke_width,'" stroke="rgb(',r,',',g,',',b,')"'
  write(iunit,'(3(a,i3),a)') &
   ' fill="rgb(',r2,',',g2,',',b2,')"/>'
!
end subroutine svg_rect
!===============================================================================
subroutine svg_circle(iunit,x,y,radius,stroke_width,r,g,b,r2,g2,b2)
!
  implicit none
  integer :: iunit,x,y,radius,stroke_width,r,g,b,r2,g2,b2
!
  write(iunit,'(3(a,i4),a)') &
   '<circle cx="',x,'" cy="',y,'" r="',radius,'"'
  write(iunit,'(a,i4,3(a,i3),a)') &
   ' stroke-width="',stroke_width,'" stroke="rgb(',r,',',g,',',b,')"'
  write(iunit,'(3(a,i3),a)') &
   ' fill="rgb(',r2,',',g2,',',b2,')"/>'
!
end subroutine svg_circle
!===============================================================================
subroutine svg_text(iunit,x,y,font_size,r,g,b,text)
!
  implicit none
  integer :: iunit,x,y,font_size,r,g,b,length
  character :: text*1000
!
  write(iunit,'(3(a,i4),a,i3,3(a,i3),a)') &
   '<text x="',x,'" y="',y,'" font-family="Helvetica" font-size="', &
   font_size,'" fill="rgb(',r,',',g,',',b,')">' 
  write(iunit,'(a)') text
  write(iunit,'(a)') '</text>'
!
end subroutine svg_text
!===============================================================================
subroutine init_random
  integer :: i, n, clock
  integer, dimension(:), allocatable :: seed
  real*8 :: r
  !
  call RANDOM_SEED(size = n)
  allocate(seed(n))
  call SYSTEM_CLOCK(count=clock)
  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  call random_seed(put = seed)
  deallocate(seed)
  call random_number(r)
  !!!!!!!
  do i=1,100
    call random_number(r)
  enddo
  !!!!!!!
end subroutine init_random
!===============================================================================
