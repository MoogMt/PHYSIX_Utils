!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine distance(r1,r2,dr,h,hi)
!implicit none
!double precision :: r1(3),r2(3),dr(3),h(3,3)
!double precision :: s1(3),s2(3),ds(3),hi(3,3)
!integer :: k
!! pass to scaled coordinates
!do k=1,3
!  s1(k)=sum(hi(k,:)*r1(:)) 
!  s2(k)=sum(hi(k,:)*r2(:)) 
!!  s1(k)=s1(k)-nint(s1(k)) ! put back into box
!!  s2(k)=s2(k)-nint(s2(k)) ! put back into box
!enddo
!! distance from minimum image
!do k=1,3
!  ds(k)=s1(k)-s2(k)
!  ds(k)=ds(k)-nint(ds(k))
!enddo
!! scale back to angstrom
!do k=1,3
!  dr(k)=sum(h(k,:)*ds(:)) 
!enddo
!!
!end subroutine distance
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine invert(a,ainv,ok_flag,det)

  implicit none
  
  double precision, dimension(3,3), intent(in)  :: a
  double precision, dimension(3,3), intent(out) :: ainv
  logical, intent(out) :: ok_flag
  
  double precision, parameter :: eps = 1.0d-10
  double precision :: det
  double precision, dimension(3,3) :: cofactor
  
  det =   a(1,1)*a(2,2)*a(3,3)  &
       - a(1,1)*a(2,3)*a(3,2)  &
       - a(1,2)*a(2,1)*a(3,3)  &
       + a(1,2)*a(2,3)*a(3,1)  &
       + a(1,3)*a(2,1)*a(3,2)  &
       - a(1,3)*a(2,2)*a(3,1)
  
  if (abs(det) .le. eps) then
     ainv = 0.0d0
     ok_flag = .false.
     return
  end if
  
  cofactor(1,1) = +(a(2,2)*a(3,3)-a(2,3)*a(3,2))
  cofactor(1,2) = -(a(2,1)*a(3,3)-a(2,3)*a(3,1))
  cofactor(1,3) = +(a(2,1)*a(3,2)-a(2,2)*a(3,1))
  cofactor(2,1) = -(a(1,2)*a(3,3)-a(1,3)*a(3,2))
  cofactor(2,2) = +(a(1,1)*a(3,3)-a(1,3)*a(3,1))
  cofactor(2,3) = -(a(1,1)*a(3,2)-a(1,2)*a(3,1))
  cofactor(3,1) = +(a(1,2)*a(2,3)-a(1,3)*a(2,2))
  cofactor(3,2) = -(a(1,1)*a(2,3)-a(1,3)*a(2,1))
  cofactor(3,3) = +(a(1,1)*a(2,2)-a(1,2)*a(2,1))
  
  ainv = transpose(cofactor) / det
  
  ok_flag = .true.
  
  return
  
end subroutine invert
