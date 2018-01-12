!----------------------------------------------------------------------!
!                                                                      !
!----------------------------------------------------------------------!
subroutine read_blocks(mpisize,svf_basename,step,step_max,reduced_sorted_distance_vector,n_steps,vec_size,max_distance,min_distance)
  implicit none
  integer:: mpisize,step,step_max
  character(len=*)::svf_basename
  integer:: unitr,i,j,tmpint,vec_size,n_steps,real_step
#ifdef STREAM
  integer:: pos
#else
  integer:: record
#endif
  logical::step_read
  character(len=200)::fn
  double precision:: max_distance,min_distance,shortint2real


  integer::sizeofshortint
#ifdef FOURBYTES
  parameter(sizeofshortint=4)
  integer(kind=4),dimension(n_steps,vec_size),intent(inout)::reduced_sorted_distance_vector
#else
  parameter(sizeofshortint=2)
  integer(kind=2),dimension(n_steps,vec_size),intent(inout)::reduced_sorted_distance_vector
#endif


  do i=1,mpisize
    unitr=1000+i
    write(fn,*), i
    fn=trim(adjustl(svf_basename))//trim(adjustl(fn))

#ifdef STREAM
    open(unit=unitr,file=trim(adjustl(fn)), form='unformatted',access='stream')
  enddo
  do i=1,n_steps
    real_step=step+(i-1)
    if(real_step.le.step_max) then
      pos=(floor((real_step-0.001)/mpisize))*vec_size*sizeofshortint+1
      unitr=1001+mod(step+(i-1)-1,mpisize)
      do j=1,vec_size
        read(unitr,pos=pos) reduced_sorted_distance_vector(i,j)
        pos=pos+sizeofshortint
      enddo
    else
      exit
    endif
  enddo
#else
    open(unit=unitr,file=trim(adjustl(fn)), form='unformatted',access='direct',recl=vec_size*sizeofshortint)
  enddo
  do i=1,n_steps
    real_step=step+(i-1)
    if(real_step.le.step_max) then
      record=floor((real_step-0.001)/mpisize)+1
      unitr=1001+mod(step+(i-1)-1,mpisize)
      read(unitr,rec=record), (reduced_sorted_distance_vector(i,j),j=1,vec_size)
      ! DEBUG print'(a,15f8.3)','XXX ',(shortint2real(reduced_sorted_distance_vector(i,j),max_distance,min_distance),j=1,vec_size)
    else
      exit
    endif
  enddo
#endif

  do i=1,mpisize
    close(1000+i)
  enddo

end subroutine read_blocks
!------------------------------------------------------------------------------------------!
