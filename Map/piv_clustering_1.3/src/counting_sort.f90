! simple and effective sorting method, useful when the possible values 
! of the array are integer and do not span a huge range (2 bytes is good, 4 is a lot...)

subroutine counting_sort(n,array, tmin, tmax)
  implicit none
    integer::sizeofshortint
#ifdef FOURBYTES
    parameter(sizeofshortint=4)
#else
    parameter(sizeofshortint=2)
#endif
  integer, intent(in) :: n
  integer(kind=sizeofshortint), dimension(n), intent(inout) :: array
  integer(kind=sizeofshortint), intent(in) :: tmin, tmax

  integer, dimension(tmin:tmax) :: cnt
  integer :: i, z

  cnt=0
  do i=1,size(array)
    cnt(array(i))=cnt(array(i))+1
  enddo

  z = 1
  
  do i = tmin, tmax
     do while ( cnt(i) > 0 )
        array(z) = i
        z = z + 1
        cnt(i) = cnt(i) - 1
     enddo
  enddo
  
end subroutine counting_sort


