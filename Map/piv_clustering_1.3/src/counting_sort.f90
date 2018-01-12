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

!  forall(i=tmin:tmax)
!     cnt(i) = count(array == i)
!  end forall

  z = 1
  do i = tmin, tmax
     do while ( cnt(i) > 0 )
        array(z) = i
        z = z + 1
        cnt(i) = cnt(i) - 1
     enddo
  enddo
end subroutine counting_sort


!program c
!  implicit none
!  integer::sizeofshortint
!  parameter(sizeofshortint=4)
!  integer(kind=sizeofshortint), dimension(64620) :: array
!  integer :: i,t1,t2
!  real :: r
!
!  do i=1,size(array)
!    call random_number(r)
!    array(i)=nint(r*10000)
!  enddo
!
!  call system_clock(t1)
!
!  do i=1,10
!  call counting_sort(size(array),array,minval(array),maxval(array))
!  enddo
!
!  call system_clock(t2)
!  write(*,*) "timing =",t2-t1
!
!  write(99,'(i6)') array
!
!end program c
