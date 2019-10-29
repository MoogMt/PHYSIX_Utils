module xyz

contains
! ------------------------------------------------------------------------------
! Read and write XYZ files
! ------------------------------------------------------------------------------
subroutine readxyztraj( file_path, str_len, n_steps, n_atoms, traj )

  implicit none

  integer,intent(in)::str_len
  character(len=str_len),intent(in):: file_path
  integer,intent(in):: n_steps
  integer,intent(in):: n_atoms

  integer:: atom, step
  character(len=4) :: dummy_name

  double precision, dimension(:,:,:),allocatable,intent(out):: traj

  allocate( traj(n_steps, n_atoms, 3) )

  open(24,file=file_path,status="old")
  do step=1,n_steps
    ! Read the two dummy lines
    read(24,*)
    read(24,*)
    ! Read one step worth of trajectory
    do atom=1,n_atoms
      read(24,*) dummy_name, traj(step,atom,1),traj(step,atom,2),traj(step,atom,3)
    enddo
  enddo
  close(24)

end subroutine readxyztraj

end module xyz
