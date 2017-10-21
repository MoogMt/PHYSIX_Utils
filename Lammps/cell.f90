program cell
  !============
  ! Variables
  !======================================
  ! File Status
  integer :: status
  ! Max and Min of the box
  real :: xmax=0, xmin=0
  real :: ymax=0, ymin=0
  real :: zmax=0, zmin=0
  ! Position variables
  real :: x,y,z
  !=====================================

  ! Opening files
  open(1,file="temp")
  open(2,file="cell")

  do
     ! Reading file
     read(1,*,iostat=status) x, y, z
     ! If EOF exit
     if( status < 0 ) exit
     if( x > xmax ) then
        xmax = x
     endif
     if( y > ymax ) then
        ymax = y
     endif
     if( z > zmax ) then
        zmin = z
     endif 
     if( x < xmin ) then
        xmin = x
     endif
     if( y < ymin ) then
        ymin = y
     endif
     if( z < zmin ) then
        zmin = z
     endif 
  enddo

  ! Write results
  write(2,*) xmin, xmax, "xlo, xhi"
  write(2,*) ymin, ymax, "ylo, yhi"
  write(2,*) zmin, zmax, "zlo, zhi"
  
  close(1)
  close(2)
  
endprogram cell
