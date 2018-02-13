program essai
implicit none
real(8) :: x, dx, laguerre
integer :: n, p, i

dx = 0.1
print '("Enter n: ",$)' ; read *, n
do i = -20, 100
 write(10,*) i*dx, laguerre(n,0,i*dx), laguerre(n,1,i*dx), laguerre(n,2,i*dx)
enddo
end program essai

include 'laguerre.f90'
