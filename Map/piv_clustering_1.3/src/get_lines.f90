!----------------------------------------------------------------------!
! Get the number of lines in a file
!----------------------------------------------------------------------!
function get_N_lines(unit_number)
  implicit none

  integer::unit_number,get_N_lines
  integer::ios = 0, N_lines = 0

  rewind(unit_number)

  do while(ios == 0)
    read(unit_number,*,iostat=ios)
    if(ios == 0) N_lines = N_lines +1
  enddo

  if(N_lines == 0) then
    write(*,'(a,i4)') 'ERROR : non-existing or empty file in unit :', unit_number
    stop
  endif

  get_N_lines = N_lines

  rewind(unit_number)
end function get_N_lines
!----------------------------------------------------------------------!
