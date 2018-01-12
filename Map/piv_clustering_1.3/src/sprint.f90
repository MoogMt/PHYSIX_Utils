subroutine sprint(n,contacts,v)
!
! it returns (non-sorted) "sprint" coordinates
! see F.Pietrucci & W.Andreoni, Phys.Rev.Lett. 107, 085504, 2011.
!
  implicit none
  integer::n,it
  double precision::contacts(n,n),v(n)
  double precision::lastv(n),diff
  double precision,parameter::eps=1.d-10
!
  ! initial guess of principal eigenvector
  v(1:n)=1.d0
  ! power method
  lastv(1:n)=0.d0
  do it=1,1000
    v=matmul(contacts,v)
    v=v/sqrt(sum(v(1:n)**2))
    diff=sum((v(1:n)-lastv(1:n))**2)
    !debug!    write(*,'(a,e16.6)') "diff",diff
    if(diff.le.eps)exit
    lastv=v
  enddo
  if(diff.gt.eps)then
    write(*,*) "ERROR: power method for SPRINT coordinates not converged after 1000 iterations: diff=",diff
  endif
  ! convergence
  v=sqrt(dble(n))*matmul(contacts,v) ! non-sorted sprint coordinates
  return
!
end subroutine sprint
