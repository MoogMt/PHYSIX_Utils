PROGRAM XYZTOGRO

  character(len=20) :: atom, str 
  integer :: i, j, t, Natoms, Nblock, Nstep
  real :: x, y, z

  Nblock = 4
  Natoms = 432
  Nstep = 1100

  open(10,file="film_run20.xyz")
  open(20,file="output.dat")


  do t=1, Nstep
    read(10,*)
    read(10,*)
    do i=1, Natoms/Nblock
      do j=1, Nblock
        read(10,*) atom, x, y, z
        if(atom=="N") then
          write(20,*) "NH1", x/10, y/10, z/10
        end if
        if(atom=="H") then
          write(20,*) "HN"//trim(str(j-1))//"", x/10, y/10, z/10
        end if
      end do
    end do
  end do

  close(10)
  close(20)

END PROGRAM XYZTOGRO

!*********************************************************************!

character(len=20) FUNCTION str(p)
    
  integer, intent(in) :: p
  write (str, *) p
  str = adjustl(str)

END FUNCTION str
