PROGRAM RANDOM

  integer :: N_blocks, N_N, N_O, N_Hn, N_Ho, N_atoms

  character (len=1), dimension(:), allocatable :: id, id_new
  real, dimension(:,:), allocatable :: pos, Hn, Ho, O, N

  N_blocks = 16
  N_N = 86
  N_Hn = 3*N_N
  N_O = 42
  N_Ho = 2*N_O
  N_atoms = N_N + N_Hn + N_O + N_Ho

  allocate(id(N_N+N_O))
  id(:) = ""
  allocate(id_new(N_N+N_O))
  id_new(:) = ""
  allocate(pos(N_N+N_O,3))
  pos(:,:) = 0.
  allocate(N(N_N,3))
  N(:,:) = 0.
  allocate(O(N_O,3))
  O(:,:) = 0.
  allocate(Hn(N_Hn,3))
  Hn(:,:) = 0.
  allocate(Ho(N_Ho,3))
  Ho(:,:) = 0. 

  CALL READ(N_atoms, N_N, N_O, pos, id)

  CALL RANDOM_NO(N_atoms, N_N, N_O, pos, id, id_new, N, O)

  CALL RANDOM_Ho(N_O, N_Ho, Ho, O)

  CALL RANDOM_Hn(N_N, N_Hn, Hn, N)

  CALL PRINT_CELL(N_blocks, N_N, N_Hn, N_O, N_Ho, O, Ho, N, Hn)

  CALL PRINT_GRO(N_blocks, N_N, N_Hn, N_O, N_Ho, O, Ho, N, Hn)

END PROGRAM RANDOM

!*****************************************************************************!

SUBROUTINE READ(N_atoms, N_N, N_O, pos, id)

  integer, intent(in) :: N_atoms, N_N, N_O
  
  character (len=1), dimension(N_N+N_O), intent(out) :: id
  real, dimension(N_N+N_O,3), intent(out) :: pos 

  character (len=10) :: atom, a
  integer :: i, j, b
  real :: x, y, z 

  open(10,file="start.gro")

  read(10,*)
  read(10,*)

  j = 0

  do i=1, N_atoms

    read(10,*) a, atom, b, x, y, z
 
    if(atom == "NH1") then
      j = j + 1
      id(j) = "N"
      pos(j,1) = x
      pos(j,2) = y
      pos(j,3) = z
    end if
     
    if(atom == "OH1") then
      j = j + 1
      id(j) = "O"
      pos(j,1) = x
      pos(j,2) = y
      pos(j,3) = z
    end if

  end do

  close(10)

END SUBROUTINE READ

!*****************************************************************************!

SUBROUTINE RANDOM_NO(N_atoms, N_N, N_O, pos, id, id_new, N, O)

  integer, intent(in) :: N_atoms, N_N, N_O
  integer :: i, clock, j, j_N, j_O
  integer, dimension(12) :: seed
  real :: r

  character (len=1), dimension(N_N+N_O), intent(in) :: id
  real, dimension(N_N+N_O,3), intent(in) :: pos 
  character (len=1), dimension(N_N+N_O), intent(out) :: id_new
  real, dimension(N_N,3), intent(out) :: N
  real, dimension(N_O,3), intent(out) :: O

  CALL SYSTEM_CLOCK(COUNT=clock)
  seed = clock
  CALL RANDOM_SEED(PUT=seed)

  j = 0
  j_N = 0
  j_O = 0 
 
  do while(j_N+j_O <= N_N+N_O)

    CALL RANDOM_NUMBER(r)
    rand = int(r*(4-1)) + 1

    if(rand == 1) then
      if(j_N < N_N) then
        j = j + 1
        id_new(j) = "N"
        j_N = j_N + 1
      end if
    end if
    if(rand == 2) then
      if(j_O < N_O) then
        j = j + 1
        id_new(j) = "O"
        j_O = j_O + 1
      end if
    end if
    if(rand == 3) then
      if(j_N < N_N) then
        j = j + 1
        id_new(j) = "N"
        j_N = j_N + 1
      end if
    end if

    if(j_N+j_O == N_N+N_O) exit

  end do

  j_N = 0
  j_O = 0

  do i=1, N_N+N_O
    if(id_new(i) == "N") then
      j_N = j_N + 1
      N(j_N,1) = pos(i,1)
      N(j_N,2) = pos(i,2)
      N(j_N,3) = pos(i,3)
    end if
    if(id_new(i) == "O") then
      j_O = j_O + 1
      O(j_O,1) = pos(i,1)
      O(j_O,2) = pos(i,2)
      O(j_O,3) = pos(i,3)
    end if
  end do

END SUBROUTINE RANDOM_NO

!*****************************************************************************!

SUBROUTINE RANDOM_Ho(N_O, N_Ho, Ho, O)

  integer, intent(in) :: N_O, N_Ho
  real, dimension(N_O,3), intent(in) :: O
  real, dimension(N_Ho,3), intent(out) :: Ho

  integer :: i, j, clock
  integer, dimension(12) :: seed
  real :: r

  call SYSTEM_CLOCK(COUNT=clock)
  seed = clock
  call RANDOM_SEED(put=seed)

  do i=1, N_O

    CALL RANDOM_NUMBER(r)
    j = int(r*(7-1)) + 1
      
      if(j == 1) then
  
        Ho((i-1)*2+1,1) = O(i,1) + 0.059
        Ho((i-1)*2+1,2) = O(i,2)
        Ho((i-1)*2+1,3) = O(i,3) + 0.076
        

        Ho((i-1)*2+2,1) = O(i,1) + 0.059
        Ho((i-1)*2+2,2) = O(i,2)
        Ho((i-1)*2+2,3) = O(i,3) - 0.076
  
      end if
 
      ! -x
      if(j == 2) then

        Ho((i-1)*2+1,1) = O(i,1) - 0.059
        Ho((i-1)*2+1,2) = O(i,2)
        Ho((i-1)*2+1,3) = O(i,3) + 0.076

        Ho((i-1)*2+2,1) = O(i,1) - 0.059
        Ho((i-1)*2+2,2) = O(i,2)
        Ho((i-1)*2+2,3) = O(i,3) - 0.076

      end if

      ! +y
      if(j == 3) then

        Ho((i-1)*2+1,1) = O(i,1) 
        Ho((i-1)*2+1,2) = O(i,2) + 0.059
        Ho((i-1)*2+1,3) = O(i,3) + 0.076

        Ho((i-1)*2+2,1) = O(i,1) 
        Ho((i-1)*2+2,2) = O(i,2) + 0.059
        Ho((i-1)*2+2,3) = O(i,3) - 0.076

      end if
  
      ! -y
      if(j == 4) then
 
        Ho((i-1)*2+1,1) = O(i,1)
        Ho((i-1)*2+1,2) = O(i,2) - 0.059
        Ho((i-1)*2+1,3) = O(i,3) + 0.076

        Ho((i-1)*2+2,1) = O(i,1) 
        Ho((i-1)*2+2,2) = O(i,2) - 0.059
        Ho((i-1)*2+2,3) = O(i,3) - 0.076

      end if

      ! +z
      if(j == 5) then

        Ho((i-1)*2+1,1) = O(i,1) - 0.076
        Ho((i-1)*2+1,2) = O(i,2)
        Ho((i-1)*2+1,3) = O(i,3) + 0.059

        Ho((i-1)*2+2,1) = O(i,1) + 0.076
        Ho((i-1)*2+2,2) = O(i,2) 
        Ho((i-1)*2+2,3) = O(i,3) + 0.059

      end if

      !-z
      if(j == 6) then

        Ho((i-1)*2+1,1) = O(i,1) - 0.076
        Ho((i-1)*2+1,2) = O(i,2)
        Ho((i-1)*2+1,3) = O(i,3) - 0.059

        Ho((i-1)*2+2,1) = O(i,1) + 0.076
        Ho((i-1)*2+2,2) = O(i,2) 
        Ho((i-1)*2+2,3) = O(i,3) - 0.059

    end if

  end do

END SUBROUTINE RANDOM_Ho

!*****************************************************************************!

SUBROUTINE RANDOM_Hn(N_N, N_Hn, Hn, N)

  integer, intent(in) :: N_N, N_Hn
  real, dimension(N_N,3), intent(in) :: N
  real, dimension(N_Hn,3), intent(out) :: Hn

  integer :: i, j, clock
  real :: r

  integer, dimension(12) :: seed

  call system_clock(COUNT=clock)
  seed = clock
  call random_seed(put=seed)

  do i=1, N_N

    CALL RANDOM_NUMBER(r)
    j = int(r*(9-1)) + 1
 
      if(j == 1) then
        Hn((i-1)*3+1,1) = 0.061 + N(i,1)
        Hn((i-1)*3+1,2) = 0.061 + N(i,2)
        Hn((i-1)*3+1,3) = -0.061 + N(i,3)        
         
        Hn((i-1)*3+2,1) = 0.061 + N(i,1)
        Hn((i-1)*3+2,2) = -0.061 + N(i,2)
        Hn((i-1)*3+2,3) = 0.061 + N(i,3)
  
        Hn((i-1)*3+3,1) = -0.061 + N(i,1)
        Hn((i-1)*3+3,2) = -0.061 + N(i,2)
        Hn((i-1)*3+3,3) = -0.061 + N(i,3)
      end if

      if(j == 2) then
        Hn((i-1)*3+1,1) = 0.061 + N(i,1)
        Hn((i-1)*3+1,2) = 0.061 + N(i,2)
        Hn((i-1)*3+1,3) = -0.061 + N(i,3)        
         
        Hn((i-1)*3+2,1) = 0.061 + N(i,1)
        Hn((i-1)*3+2,2) = -0.061 + N(i,2)
        Hn((i-1)*3+2,3) = 0.061 + N(i,3)
  
        Hn((i-1)*3+3,1) = -0.061 + N(i,1)
        Hn((i-1)*3+3,2) = 0.061 + N(i,2)
        Hn((i-1)*3+3,3) = 0.061 + N(i,3)
      end if

      if(j == 3) then
        Hn((i-1)*3+1,1) = 0.061 + N(i,1)
        Hn((i-1)*3+1,2) = 0.061 + N(i,2)
        Hn((i-1)*3+1,3) = -0.061 + N(i,3)

        Hn((i-1)*3+2,1) = -0.061 + N(i,1)
        Hn((i-1)*3+2,2) = -0.061 + N(i,2)
        Hn((i-1)*3+2,3) = -0.061 + N(i,3)
    
        Hn((i-1)*3+3,1) = -0.061 + N(i,1)
        Hn((i-1)*3+3,2) = 0.061 + N(i,2)
        Hn((i-1)*3+3,3) = 0.061 + N(i,3)
      end if

      if(j == 4) then
        Hn((i-1)*3+1,1) = 0.061 + N(i,1)
        Hn((i-1)*3+1,2) = -0.061 + N(i,2)
        Hn((i-1)*3+1,3) = 0.061 + N(i,3)
  
        Hn((i-1)*3+2,1) = -0.061 + N(i,1)
        Hn((i-1)*3+2,2) = -0.061 + N(i,2)
        Hn((i-1)*3+2,3) = -0.061 + N(i,3)
    
        Hn((i-1)*3+3,1) = -0.061 + N(i,1)
        Hn((i-1)*3+3,2) = 0.061 + N(i,2)
        Hn((i-1)*3+3,3) = 0.061 + N(i,3)
      end if

      if(j == 5) then
        Hn((i-1)*3+1,1) = 0.061 + N(i,1)
        Hn((i-1)*3+1,2) = -0.061 + N(i,2)
        Hn((i-1)*3+1,3) = -0.061 + N(i,3)
  
        Hn((i-1)*3+2,1) = -0.061 + N(i,1)
        Hn((i-1)*3+2,2) = 0.061 + N(i,2)
        Hn((i-1)*3+2,3) = -0.061 + N(i,3)
   
        Hn((i-1)*3+3,1) = -0.061 + N(i,1)
        Hn((i-1)*3+3,2) = -0.061 + N(i,2)
        Hn((i-1)*3+3,3) = 0.061 + N(i,3)
      end if

      if(j == 6) then
        Hn((i-1)*3+1,1) = 0.061 + N(i,1)
        Hn((i-1)*3+1,2) = 0.061 + N(i,2)
        Hn((i-1)*3+1,3) = 0.061 + N(i,3)
 
        Hn((i-1)*3+2,1) = -0.061 + N(i,1)
        Hn((i-1)*3+2,2) = 0.061 + N(i,2)
        Hn((i-1)*3+2,3) = -0.061 + N(i,3)
   
        Hn((i-1)*3+3,1) = -0.061 + N(i,1)
        Hn((i-1)*3+3,2) = -0.061 + N(i,2)
        Hn((i-1)*3+3,3) = 0.061 + N(i,3)
      end if

      if(j == 7) then
        Hn((i-1)*3+1,1) = 0.061 + N(i,1)
        Hn((i-1)*3+1,2) = 0.061 + N(i,2)
        Hn((i-1)*3+1,3) = 0.061 + N(i,3)      
  
        Hn((i-1)*3+2,1) = 0.061 + N(i,1)
        Hn((i-1)*3+2,2) = -0.061 + N(i,2)
        Hn((i-1)*3+2,3) = -0.061 + N(i,3)

        Hn((i-1)*3+3,1) = -0.061 + N(i,1)
        Hn((i-1)*3+3,2) = -0.061 + N(i,2)
        Hn((i-1)*3+3,3) = 0.061 + N(i,3)
      end if

      if(j == 8) then
        Hn((i-1)*3+1,1) = 0.061 + N(i,1)
        Hn((i-1)*3+1,2) = 0.061 + N(i,2)
        Hn((i-1)*3+1,3) = 0.061 + N(i,3)        
         
        Hn((i-1)*3+2,1) = 0.061 + N(i,1)
        Hn((i-1)*3+2,2) = -0.061 + N(i,2)
        Hn((i-1)*3+2,3) = -0.061 + N(i,3)
  
        Hn((i-1)*3+3,1) = -0.061 + N(i,1)
        Hn((i-1)*3+3,2) = 0.061 + N(i,2)
        Hn((i-1)*3+3,3) = -0.061 + N(i,3)
      end if

  end do

END SUBROUTINE RANDOM_Hn

!*****************************************************************************!

SUBROUTINE PRINT_CELL(N_blocks, N_N, N_Hn, N_O, N_Ho, O, Ho, N, Hn)

  integer, intent(in) :: N_blocks, N_O, N_N, N_Hn, N_Ho
  real, dimension(N_O,3), intent(in) :: O
  real, dimension(N_Ho,3), intent(in) :: Ho
  real, dimension(N_N,3), intent(in) :: N
  real, dimension(N_Hn,3), intent(in) :: Hn

  character (len=10) :: a, atom_t0
  integer :: i, j, b, j_N, j_Hn, j_O, j_Ho

  j_N = 0
  j_Hn = 0
  j_O = 0
  j_Ho = 0

  open(10,file="start.gro")
  open(20,file="random_pos.cell")  

  read(10,*)
  read(10,*)

  write(20,*) "%BLOCK LATTICE_ABC"
  write(20,*) "13.3088 13.3088 13.3088"
  write(20,*) "90 90 90"
  write(20,*) "%ENDBLOCK LATTICE_ABC"

  write(20,*) " "

  write(20,*) "%BLOCK POSITIONS_ABS"


  Do i=1, N_blocks-2

    do j=1, 24
    
      read(10,*) a, atom_t0, b

      if(atom_t0 =="NH1") then
        j_N = j_N + 1
        write(20, *) "N", N(j_N,1)*10, N(j_N,2)*10, N(j_N,3)*10
      end if

     if(atom_t0 == "HN1" .or. & 
        atom_t0 == "HN2" .or. &
        atom_t0 == "HN3") then
        j_Hn = j_Hn + 1
        write(20,*) "H", Hn(j_Hn,1)*10, Hn(j_Hn,2)*10, Hn(j_Hn,3)*10
      end if

    end do

    do j=25, 30

      read(10,*) a, atom_t0, b

      if(atom_t0 == "OH1") then
        j_O = j_O + 1
        write(20,*) "O", O(j_O,1)*10, O(j_O,2)*10, O(j_O,3)*10
      end if

      if(atom_t0 == "HO1" .or. &
         atom_t0 == "HO2") then
         j_Ho = j_Ho + 1
         write(20,*) "H", Ho(j_Ho,1)*10, Ho(j_Ho,2)*10, Ho(j_Ho,3)*10
      end if

    end do
    
  End Do

  Do i=N_blocks-1, N_blocks-1

    do j=1, 8  

      read(10,*) a, atom_t0, b

      if(atom_t0 =="NH1") then
        j_N = j_N + 1
        write(20, *) "N", N(j_N,1)*10, N(j_N,2)*10, N(j_N,3)*10
      end if

      if(atom_t0 == "HN1" .or. & 
        atom_t0 == "HN2" .or. &
        atom_t0 == "HN3") then
        j_Hn = j_Hn + 1
        write(20,*) "H", Hn(j_Hn,1)*10, Hn(j_Hn,2)*10, Hn(j_Hn,3)*10
      end if

    end do
   
    do j=9, 26

      read(10,*) a, atom_t0, b

      if(atom_t0 == "OH1") then
        j_O = j_O + 1
        write(20,*) "O", O(j_O,1)*10, O(j_O,2)*10, O(j_O,3)*10
      end if

      if(atom_t0 == "HO1" .or. &
         atom_t0 == "HO2") then
         j_Ho = j_Ho + 1
         write(20,*) "H", Ho(j_Ho,1)*10, Ho(j_Ho,2)*10, Ho(j_Ho,3)*10
      end if

    end do

  End Do

  Do i=N_blocks, N_blocks

    do j=1, 24
  
      read(10,*) a, atom_t0, b

        if(atom_t0 == "OH1") then
          j_O = j_O + 1
          write(20,*) "O", O(j_O,1)*10, O(j_O,2)*10, O(j_O,3)*10
        end if

        if(atom_t0 == "HO1" .or. &
          atom_t0 == "HO2") then
          j_Ho = j_Ho + 1
          write(20,*) "H", Ho(j_Ho,1)*10, Ho(j_Ho,2)*10, Ho(j_Ho,3)*10
        end if

    end do
 
  End Do

  write(20,*) "%ENDBLOCK POSITIONS_ABS"

  close(10)
  close(20)

END SUBROUTINE PRINT_CELL

!*****************************************************************************!

SUBROUTINE PRINT_GRO(N_blocks, N_N, N_Hn, N_O, N_Ho, O, Ho, N, Hn)

  integer, intent(in) :: N_blocks, N_O, N_N, N_Hn, N_Ho
  real, dimension(N_O,3), intent(in) :: O
  real, dimension(N_Ho,3), intent(in) :: Ho
  real, dimension(N_N,3), intent(in) :: N
  real, dimension(N_Hn,3), intent(in) :: Hn

  character (len=10) :: a, atom_t0
  integer :: i, j, b, j_N, j_Hn, j_O, j_Ho

  j_N = 0
  j_Hn = 0
  j_O = 0
  j_Ho = 0

  open(10,file="start.gro")
  open(30,file="random_pos.dat")  

  read(10,*)
  read(10,*)

  write(30,*) "13.3088 13.3088 13.3088"

  Do i=1, N_blocks-2

    do j=1, 24
    
      read(10,*) a, atom_t0, b

      if(atom_t0 =="NH1") then
        j_N = j_N + 1
        write(30, *) "NH1", N(j_N,1)*10, N(j_N,2)*10, N(j_N,3)*10
      end if

     if(atom_t0 == "HN1") then
        j_Hn = j_Hn + 1
        write(30,*) "HN1", Hn(j_Hn,1)*10, Hn(j_Hn,2)*10, Hn(j_Hn,3)*10
      end if

      if(atom_t0 == "HN2") then
        j_Hn = j_Hn + 1
        write(30,*) "HN2", Hn(j_Hn,1)*10, Hn(j_Hn,2)*10, Hn(j_Hn,3)*10
      end if     

      if(atom_t0 == "HN3") then
        j_Hn = j_Hn + 1
        write(30,*) "HN3", Hn(j_Hn,1)*10, Hn(j_Hn,2)*10, Hn(j_Hn,3)*10
      end if    

    end do

    do j=25, 30

      read(10,*) a, atom_t0, b

      if(atom_t0 == "OH1") then
        j_O = j_O + 1
        write(30,*) "OH1", O(j_O,1)*10, O(j_O,2)*10, O(j_O,3)*10
      end if

      if(atom_t0 == "HO1") then
        j_Ho = j_Ho + 1
        write(30,*) "HO1", Ho(j_Ho,1)*10, Ho(j_Ho,2)*10, Ho(j_Ho,3)*10
      end if

      if(atom_t0 == "HO2") then
        j_Ho = j_Ho + 1
        write(30,*) "HO2", Ho(j_Ho,1)*10, Ho(j_Ho,2)*10, Ho(j_Ho,3)*10
      end if

    end do
    
  End Do

  Do i=N_blocks-1, N_blocks-1

    do j=1, 8
    
      read(10,*) a, atom_t0, b

      if(atom_t0 =="NH1") then
        j_N = j_N + 1
        write(30, *) "NH1", N(j_N,1)*10, N(j_N,2)*10, N(j_N,3)*10
      end if

     if(atom_t0 == "HN1") then
        j_Hn = j_Hn + 1
        write(30,*) "HN1", Hn(j_Hn,1)*10, Hn(j_Hn,2)*10, Hn(j_Hn,3)*10
      end if

      if(atom_t0 == "HN2") then
        j_Hn = j_Hn + 1
        write(30,*) "HN2", Hn(j_Hn,1)*10, Hn(j_Hn,2)*10, Hn(j_Hn,3)*10
      end if     

      if(atom_t0 == "HN3") then
        j_Hn = j_Hn + 1
        write(30,*) "HN3", Hn(j_Hn,1)*10, Hn(j_Hn,2)*10, Hn(j_Hn,3)*10
      end if    

    end do

    do j=9, 26

      read(10,*) a, atom_t0, b

      if(atom_t0 == "OH1") then
        j_O = j_O + 1
        write(30,*) "OH1", O(j_O,1)*10, O(j_O,2)*10, O(j_O,3)*10
      end if

      if(atom_t0 == "HO1") then
        j_Ho = j_Ho + 1
        write(30,*) "HO1", Ho(j_Ho,1)*10, Ho(j_Ho,2)*10, Ho(j_Ho,3)*10
      end if

      if(atom_t0 == "HO2") then
        j_Ho = j_Ho + 1
        write(30,*) "HO2", Ho(j_Ho,1)*10, Ho(j_Ho,2)*10, Ho(j_Ho,3)*10
      end if

    end do

  End Do

  Do i=N_blocks, N_blocks

    do j=1, 24

      read(10,*) a, atom_t0, b

      if(atom_t0 == "OH1") then
        j_O = j_O + 1
        write(30,*) "OH1", O(j_O,1)*10, O(j_O,2)*10, O(j_O,3)*10
      end if

      if(atom_t0 == "HO1") then
        j_Ho = j_Ho + 1
        write(30,*) "HO1", Ho(j_Ho,1)*10, Ho(j_Ho,2)*10, Ho(j_Ho,3)*10
      end if

      if(atom_t0 == "HO2") then
        j_Ho = j_Ho + 1
        write(30,*) "HO2", Ho(j_Ho,1)*10, Ho(j_Ho,2)*10, Ho(j_Ho,3)*10
      end if

    end do

  End Do

  close(10)
  close(30)

END SUBROUTINE print_gro
