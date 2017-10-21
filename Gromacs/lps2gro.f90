program read

  !========
  ! ATOM I
  !========================================
  integer :: atomi,moli,typei,nxi,nyi,nzi
  real :: xi,yi,zi,chargei
  !========================================

  !========
  ! ATOM J
  !=============================================
  integer :: atomj, molj, typej, nxj, nyj, nzj
  real :: xj, yj, zj, chargej
  !=============================================

  !=======
  ! DUMMY
  !==============
  integer :: i,j
  !==============

  !--------
  ! OUTPUT
  !------------------------------------
  open(2,file='/home/moog/input.gro')
  !------------------------------------
  
  !--------
  ! INPUT
  !-------------------------------------
  open(1,file='/home/moog/data.start')
  !-------------------------------------

  !----------------
  ! FIrst lines
  !-----------------
  write(2,*) "CO2 Meta"
  write(2,*) 2592
  !--------------------
  
  !---------------
  ! Reading input
  !------------------------------------------------------------------------------------
  do
     read(1,*,iostat=i) atomi, moli, typei, chargei, xi, yi , zi, nxi, nyi, nzi
     if( i < 0 ) exit
     if( typei .eq. 1 ) then
         if ( xi < 0 ) x=x+30.336
         if ( yi < 0 ) y=y+30.336
         if ( zi < 0 ) z=z+30.336
        write(2,'(i5,2a5,i5,3f8.3,3f8.4)') moli, "COO", "C", atomi, xi, yi ,zi
        open(3,file='/home/moog/data.start2')
        do
           read(3,*,iostat=j) atomj, molj, typej, chargej, xj, yj , zj, nxj, nyj, nzj
           if ( j < 0 ) exit
           if ( xj < 0 ) x=x+30.336
           if ( yj < 0 ) y=y+30.336
           if ( zj < 0 ) z=z+30.336
           if ( typej .eq. 2 ) then
             if ( molj .eq. moli ) then
                 write(2,'(i5,2a5,i5,3f8.3,3f8.4)') molj, "COO", "O", atomj, xj, yj ,zj
              endif
           endif
        enddo
        close(3)
     endif
  enddo
  !------------------------------------------------------------------------------------

  !----------------
  ! Closing files
  !----------------
  close(1)
  close(2)
  !----------------
  
end program read
