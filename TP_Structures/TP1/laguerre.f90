real(8) function laguerre(n,p,x)
implicit none
integer :: n, p
real(8) :: x
integer :: k
real(8) :: lp, l, lm
if ( n == 0 ) then
   laguerre = 1._8
else if ( n == 1 ) then
   laguerre = 1._8+p-x
else
   lm = 1._8 ; l = 1._8+p-x
   do k = 2, n
     lp = ((2*k+1+p-x)*l-(k+p)*lm)/(k+1)
     lm = l ; l = lp
   enddo
   laguerre = l
endif
end function laguerre
