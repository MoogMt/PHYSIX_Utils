module clustering_mod

contains
  
  !----------------------------------------------------------------------------------!
  ! K-Medoids Algorithm
  !----------------------------------------------------------------------------------!
  subroutine kmedoids_algorithm( frame2frame, n_steps, n_clusters, cluster_size, cluster_centers, cluster_members, slow )
    
    !------------------------------------------------------------------------
    ! Algorithm from H.S.Park and C.H.Jun, Expert Syst. Appl. 36, 3336, 2009.
    !------------------------------------------------------------------------
    
    implicit none

    !------------------
    ! Choice of method
    !-----------------------------------
    logical, intent(in) :: slow  ! Slow
    !-----------------------------------
    
    ! Variables
    !---------------------------------------------------------------------------------------
    integer, intent(in)::n_steps                                          ! Number of steps
    double precision , dimension(n_steps,n_steps), intent(in)::frame2frame  ! F2F matrix
    !---------------------------------------------------------------------------------------
    
    !-----------
    ! Clusters
    !---------------------------------------------------------------------------------------
    integer,intent(in)::n_clusters                                   ! Number of clusters
    integer:: max_cluster_size                                       ! Maximum cluster size
    integer, allocatable, dimension(:), intent(out)::cluster_size       ! Cluster sizes
    integer, allocatable, dimension(:), intent(out)::cluster_centers    ! Cluster centers
    integer, allocatable, dimension(:,:), intent(out)::cluster_members  ! Clusters members
    integer, allocatable, dimension(:,:) :: cluster_members_tmp
    integer, allocatable, dimension(:) :: cluster_ranking
    integer, allocatable, dimension(:) :: cluster_centers_tmp
    !---------------------------------------------------------------------------------------

    !----
    ! ?
    !----------------------------------------------------------------
    integer :: i,j,ii,k,newcenter,i2cluster(n_steps),cluster2i(n_clusters,n_steps)
    double precision :: dmin,cost,oldcost,dtmp,r,prob(n_steps)
    integer :: good(n_steps)
    logical :: found
    integer :: ass(n_steps)
    double precision , allocatable, dimension (:,:) :: vector_temp
    !------------------------------------------------------------------------------
    
    !---------------------------------
    ! Init cluster related variables
    !------------------------------------------------------------------------------------------
    max_cluster_size=n_steps
    allocate(cluster_size(n_clusters),cluster_centers(n_clusters),cluster_members(n_clusters,max_cluster_size))
  !--------------------------------------------------------------------------------------------
  
    !----------------------
    ! Init cluster message
    !-------------------------------------------------------------------------------------------
    write(*,*) "using k-medoids algorithm ..."
    write(*,*) "  (H.S.Park and C.H.Jun, Expert Syst. Appl. 36, 3336, 2009)"
    write(*,*) "... with k-means++ initialization of centers"
    write(*,*) "  (D.Arthur and S.Vassilvitskii, Proc.18th annual ACM-SIAM symposium on Discrete algorithms, 1027, 2007)"
    !--------------------------------------------------------------------------------------------

    !---------------------
    ! Init random number
    !-------------------------
    call init_random_seed
    call random_number(r)
    !-------------------------
    
    !-------
    ! Init
    !---------------------------------------------
    cluster_centers(1)=int(r*dble(n_steps))+1
    good(1:n_steps)=1
    good(cluster_centers(1))=0
    !---------------------------------------------

    !--------------------
    ! Fast Method
    !------------------------------------------------------------------------------------------
    if ( .not. slow ) then
       !-------------------
       ! Loop over clusters
       !-------------------------------------------------------------------
       do i=2,n_clusters
          ! Read paper
          !-----------------------------------------------------------------
          do k=1,n_steps
             ! if (good(k)==0) cycle
             dmin=1.d10
             !-------------------------------------------------------------
             do j=1,i-1
                dtmp=frame2frame(k,cluster_centers(j))
                if(dtmp.lt.dmin) dtmp=dmin
             enddo
             !--------------------------------------------------------------
             prob(k)=dmin*dmin
          enddo
          !------------------------------------------------------------------
          dtmp=sum(prob(:))
          prob(:)=prob(:)/dtmp
          do k=2,n_steps
             prob(k)=prob(k)+prob(k-1)
          enddo
          found=.false.
          do while(.not.found)
             call random_number(r)
             if (r.lt.prob(1)) then
                if (good(1)==0) then
                   cycle
                endif
                cluster_centers(i)=1
                good(1)=0
                found=.true.
             else
                do k=2,n_steps
                   if (good(k)==0) cycle
                   if (r.ge.prob(k-1).and.r.lt.prob(k)) then
                      cluster_centers(i)=k
                      good(k)=0
                      found=.true.
                      exit
                   endif
                enddo
             endif
          enddo
       enddo
       !------------------------------------------------------------------
    else
       do i=2,n_clusters
          ! Lire en memoire i-1 vecteurs avec les dist
          
          ! Read paper
          do k=1,n_steps
             ! if (good(k)==0) cycle
             dmin=1.d10
             do j=1,i-1
                dtmp=frame2frame(k,cluster_centers(j))
                if(dtmp.lt.dmin) dtmp=dmin
             enddo
             prob(k)=dmin*dmin
          enddo
          dtmp=sum(prob(:))
          prob(:)=prob(:)/dtmp
          do k=2,n_steps
             prob(k)=prob(k)+prob(k-1)
          enddo
          found=.false.
          do while(.not.found)
             call random_number(r)
             if (r.lt.prob(1)) then
                if (good(1)==0) then
                   cycle
                endif
                cluster_centers(i)=1
                good(1)=0
                found=.true.
             else
                do k=2,n_steps
                   if (good(k)==0) cycle
                   if (r.ge.prob(k-1).and.r.lt.prob(k)) then
                      cluster_centers(i)=k
                      good(k)=0
                      found=.true.
                      exit
                   endif
                enddo
             endif
          enddo
       enddo
    endif
    !-------------------------------------------------------------------------------
    
    !-----------------------
    ! Out message for debug
    !-----------------------------------------------------------------------
    write(*,'(a,1000i6)') "kmeans++:",cluster_centers(1:n_clusters) !debug
    !-----------------------------------------------------------------------
    
    !------------
    ! Initialize
    !---------------
    oldcost=1.d10
    !---------------
    
    !------
    ! Loop
    !-------------------------------------------------------------------------------------------
    do
       !---------------------
       ! Voronoi assignment
       !---------------------
       do i=1,n_clusters
          cluster_size(i)=0
       enddo
       !----------------------
       
       !-----------------------------------------------------------
       do i=1,n_steps
          dmin=1.d10
          do j=1,n_clusters
             if (frame2frame(i,cluster_centers(j))<dmin) then
                dmin=frame2frame(i,cluster_centers(j))
                i2cluster(i)=j
             endif
          enddo
          j=i2cluster(i)
          cluster_size(j)=cluster_size(j)+1      
          cluster2i(j,cluster_size(j))=i
       enddo
       !-------------------------------------------------------------
       
       !---------
       cost=0.d0
       !---------
       
       !----------------------------------------------------------
       do i=1,n_steps
          cost=cost+frame2frame(i,cluster_centers(i2cluster(i)))
       enddo
       !----------------------------------------------------------
       
       !------------------------------------
       write(*,'(a,f18.8)') "  cost =",cost
       !-------------------------------------
       
       !--------------------------
       if (cost.eq.oldcost) exit
       !--------------------------
       
       !------------
       oldcost=cost
       !------------
       
       !---------------
       ! Update medoids
       !----------------------------------------------------------------------------
       do i=1,n_clusters
          dmin=1.d10
          do j=1,cluster_size(i)
             dtmp=sum(frame2frame(cluster2i(i,j),cluster2i(i,1:cluster_size(i))))
             if (dtmp<dmin) then
                dmin=dtmp
                newcenter=cluster2i(i,j)
             endif
          enddo
          cluster_centers(i)=newcenter    
       enddo
       !-----------------------------------------------------------------------------
       
    enddo
    !---------------------------------------------------------------------------------------------

    !--------------
    ! Convergence
    !------------------------------------------------------------------------------------
    write(*,*) "  converged"
    cluster_members(:,:)=cluster2i(:,:)
    !------------------------------------------------------------------------------------

    !------------------------------------------
    ! Sorting clusters from largest to smallest
    !------------------------------------------------------------------------------------
    allocate(cluster_ranking(n_clusters))
    do i=1,n_clusters
       cluster_ranking(i)=i
    enddo
    !------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------
    do i=1,n_clusters-1
       do j=i+1,n_clusters
          if(cluster_size(j).gt.cluster_size(i)) then
             ii=cluster_size(j)
             cluster_size(j)=cluster_size(i)
             cluster_size(i)=ii
             ii=cluster_ranking(j)
             cluster_ranking(j)=cluster_ranking(i)
             cluster_ranking(i)=ii
          endif
       enddo
    enddo
    !------------------------------------------------------------------------------------
    
    allocate(cluster_centers_tmp(n_clusters))
    allocate(cluster_members_tmp(n_clusters,max_cluster_size))
    cluster_centers_tmp=cluster_centers
    cluster_members_tmp=cluster_members

    !------------------------------------------------------------------------------------
    do i=1,n_clusters
       cluster_centers(i)=cluster_centers_tmp(cluster_ranking(i))
       cluster_members(i,1:max_cluster_size)=cluster_members_tmp(cluster_ranking(i),1:max_cluster_size)
    enddo
    !------------------------------------------------------------------------------------
    
    deallocate(cluster_ranking,cluster_centers_tmp,cluster_members_tmp)

    !------------------------------------------------------------------------------------
    ! check that all frames are assigned and no frame is unassigned
    !------------------------------------------------------------------------------------
    ass=0
    do i=1,n_clusters
       do j=1,cluster_size(i)
          ass(cluster_members(i,j))=ass(cluster_members(i,j))+1
       enddo
    enddo
    !------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------
    do i=1,n_steps
       if (ass(i).ne.1) then
          write(*,*) "ERROR in clustering: frame",i," is assigned ",ass(i)," times !!"
          stop
       endif
    enddo
    !--------------------------------------------------------------------------------------

    !------------------
    write(*,*) "DONE"
    !------------------

    !-----------------------
    ! Clusters Log writting
    !--------------------------------------------------------------------------
    open(99,file="log-clusters",status="unknown")
    write(99,'(a,1000i5)') "centers: ",cluster_centers
    do i=1,n_clusters
       write(99,'(i5,a,1000i5)') i,":  ",cluster_members(i,1:cluster_size(i))
    enddo
    close(99)
    !--------------------------------------------------------------------------
    
  end subroutine kmedoids_algorithm
  !----------------------------------------------------------------------------------------------!

  !==============================================================================================!
  
  !----------------------------------------------------------------------------------------------!
  ! Daura Algorithm
  !----------------------------------------------------------------------------------------------!
  subroutine daura_algorithm( frame2frame_orig, n_steps, n_clusters,actual_cluster_size, cluster_centers,actual_cluster_members, cutoff, slow )
    
    !------------------------------------------------------------------------
    ! Algorithm from Daura et al, Angew. Chem. Int. Ed. 38, 236-240, 1999
    !------------------------------------------------------------------------
    implicit none
    
    integer,intent(in)::n_steps
    double precision,intent(in)::cutoff
    double precision,dimension(n_steps,n_steps),intent(inout)::frame2frame_orig
    integer,intent(out)::n_clusters
    integer,allocatable,dimension(:),intent(out)::actual_cluster_size
    integer,allocatable,dimension(:),intent(out)::cluster_centers
    integer,allocatable,dimension(:,:),intent(out)::actual_cluster_members
    integer:: max_n_clusters,max_cluster_size,n_neighbours,m_max,n_frames_left,max_neighbours,pcs,acs
    integer,allocatable,dimension(:)::possible_cluster_size,cluster_ranking,cluster_centers_tmp
    integer,allocatable,dimension(:,:)::possible_cluster_members
    double precision,allocatable,dimension(:,:)::frame2frame,frame2frame_save,possible_n2cluster
    double precision,allocatable,dimension(:,:)::closest_center
    double precision:: d
    integer::i,j,l,n,m,ii,jj,jj1,jj2,tmp_int
    
    max_n_clusters=n_steps
    max_cluster_size=n_steps
    
    allocate(cluster_centers(max_n_clusters))
    
    allocate(possible_n2cluster(n_steps,max_n_clusters))
    allocate(possible_cluster_members(max_n_clusters,max_cluster_size))
    allocate(possible_cluster_size(max_n_clusters))
    allocate(closest_center(n_steps,2)) ! for each frame (index_closest_center, value)
    allocate(frame2frame(n_steps,n_steps),frame2frame_save(n_steps,n_steps))
    
    frame2frame=frame2frame_orig
    
    possible_n2cluster=-0.9d0
    closest_center=100000000.d0
    
    write(*,*) "using daura algorithm + final Voronoi assignment"
    write(*,*) "  (Daura et al, Angew. Chem. Int. Ed. 38, 236-240, 1999)"

    ! frame2frame will change depending on the clusterization process (neighbours will be removed once they attributed to a center)
    ! frame2frame_save doesn't change; once clusterization is finished and cluster centers have been determined, frame2frame_save will serve
    ! - to attribute each frame to the closest cluster center (which might be different from the one it was first attributed to)
    do m=1,n_steps
       do n=m+1,n_steps
          if(frame2frame(m,n).gt.cutoff) then
             frame2frame(m,n)=-0.d0
             frame2frame(n,m)=-0.d0
          endif
          frame2frame_save(m,n)=frame2frame(m,n)
          frame2frame_save(n,m)=frame2frame(n,m)
       enddo
       frame2frame(m,m)=1.d0
       frame2frame_save(m,m)=-0.9d0
    enddo

    n_frames_left=n_steps
    n_clusters=0
    print '(a,$)', 'Clustering...'
    do while (n_frames_left.gt.0)
       max_neighbours=-1
       m_max=0
       do m=1,n_steps
          if(frame2frame(m,m).gt.0.d0) then
             n_neighbours=count(frame2frame(m,:).gt.0.d0)-1 ! -1 because frame2frame(m,m) is not really a neighbour
             if(n_neighbours.gt.max_neighbours) then
                m_max=m
                max_neighbours=n_neighbours
             endif
          endif
       enddo
       !From here the cluster center IS chosen; the next loop only determines ALL its possible neighbours
       if(max_neighbours.ge.0) then
          n_clusters=n_clusters+1
          n_frames_left=n_frames_left-1 ! also decremented lower in the loop
          cluster_centers(n_clusters)=m_max
          closest_center(m_max,1)=dble(m_max)+0.1d0 ! so that when int() it really goes to the correct index
          closest_center(m_max,2)=0.d0
          possible_cluster_members(n_clusters,1)=m_max
          frame2frame(m_max,m_max)=-0.9d0
          pcs=1
          do n=1,n_steps
             if(frame2frame_save(m_max,n).gt.0.d0) then  !frame2frame_save(m_max,m_max) was initialized to 0
                if(frame2frame(m_max,n).gt.0.d0) then
                   n_frames_left=n_frames_left-1
                endif
                d=frame2frame_save(m_max,n)
                pcs=pcs+1
                possible_cluster_members(n_clusters,pcs)=n
                possible_n2cluster(n,n_clusters)=d
                if(d.lt.closest_center(n,2)) then   !checks if the new center is the closest center
                   closest_center(n,1)=dble(n_clusters)+0.1d0
                   closest_center(n,2)=d
                endif
                frame2frame(n,:)=-0.9d0   !remove the neigbours used to determine this center form the neighbour list
                frame2frame(:,n)=-0.9d0
             endif
          enddo
          possible_cluster_size(n_clusters)=pcs
          frame2frame(m_max,:)=-0.9d0
          frame2frame(:,m_max)=-0.9d0
          !print*,'NEW STEP around',m_max
          !do m=1,n_steps
          !  print '(i3,a,$)',m,' : '
          !  do n=1,n_steps
          !    print '(i3,$)',int(frame2frame(m,n))
          !  enddo
          !  print*,
          !enddo
       endif
    enddo
    !Restore frame2frame before returning
    do m=1,n_steps
       frame2frame(m,m)=1.d0
       do n=m+1,n_steps
          frame2frame(m,n)=frame2frame_save(m,n)
          frame2frame(n,m)=frame2frame_save(n,m)
       enddo
    enddo
    deallocate(possible_n2cluster,frame2frame_save)
    print*, 'DONE'
    
    ! Here we compute the ACTUAL cluster size; meaning we attribute the frame to the closest center (daura)
    allocate(actual_cluster_size(n_clusters))
    allocate(actual_cluster_members(n_clusters,possible_cluster_size(1)))
    actual_cluster_size=1
    max_n_clusters=0
    do i=1,n_clusters
       pcs=possible_cluster_size(i)
       actual_cluster_members(i,1)=cluster_centers(i)
       do j=2,pcs
          jj=possible_cluster_members(i,j)
          if(int(closest_center(jj,1)).eq.i) then
             actual_cluster_size(i)=actual_cluster_size(i)+1
             actual_cluster_members(i,actual_cluster_size(i))=jj
          endif
       enddo
       if(actual_cluster_size(i).gt.max_n_clusters) max_n_clusters=actual_cluster_size(i)
    enddo
    deallocate(possible_cluster_members,possible_cluster_size)
    
    !here we sort the actual_cluster_members inside each cluster by increasing distance from the center
    do i=1,n_clusters
       acs=actual_cluster_size(i)
       do ii=2,acs-1
          do jj=ii+1,acs
             jj1=actual_cluster_members(i,ii) 
             jj2=actual_cluster_members(i,jj) 
             if (closest_center(jj1,2) > closest_center(jj2,2)) then ! the second dimention of closest_center indicates the distance from the center.
                actual_cluster_members(i,jj)=jj1
                actual_cluster_members(i,ii)=jj2
             endif
          enddo
       enddo
    enddo
    deallocate(closest_center)
    
    ! Here we rearange the clusters so as to have the biggest one as the first cluster
    ! First step is to compute the ranking and exchange the size; the ranking is used after to make the matching of the cluster_members and cluster_centers themselves
    allocate(cluster_ranking(n_clusters))
    do i=1,n_clusters
       cluster_ranking(i)=i
    enddo
    
    do i=1,n_clusters-1
       do j=i+1,n_clusters
          if(actual_cluster_size(j).gt.actual_cluster_size(i)) then
             tmp_int=actual_cluster_size(j)
             actual_cluster_size(j)=actual_cluster_size(i)
             actual_cluster_size(i)=tmp_int
             tmp_int=cluster_ranking(j)
             cluster_ranking(j)=cluster_ranking(i)
             cluster_ranking(i)=tmp_int
          endif
       enddo
    enddo

    !Then exchange both cluster_centers and actual_cluster_members
    allocate(possible_cluster_members(n_clusters,max_n_clusters))
    allocate(cluster_centers_tmp(n_clusters))
    cluster_centers_tmp=cluster_centers
    do i=1,n_clusters
       possible_cluster_members(i,1:max_n_clusters)=actual_cluster_members(cluster_ranking(i),1:max_n_clusters)
       cluster_centers(i)=cluster_centers_tmp(cluster_ranking(i))
    enddo
    deallocate(actual_cluster_members,cluster_ranking,cluster_centers_tmp)
    allocate(actual_cluster_members(n_clusters,max_n_clusters))
    
    actual_cluster_members=possible_cluster_members
    deallocate(possible_cluster_members)
    
    open(99,file="log-clusters",status="unknown")
    write(99,'(a,1000i5)') "centers: ",cluster_centers(1:n_clusters)
    do i=1,n_clusters
       write(99,'(i5,a,1000i5)') i,":  ",actual_cluster_members(i,1:actual_cluster_size(i))
    enddo
    close(99)
    
    return
    
end subroutine daura_algorithm
!------------------------------------------------------------------------------------------------!

!-----------------------------------
! Computing clustering coefficients
!------------------------------------------------------------------------------------------------!
subroutine compute_clustering_coefficients(frame2frame,n_steps, n_clusters,cluster_members,cluster_size, clustering_coefficients,cutoff_clcoeff)
    implicit none
    integer,intent(in)::n_clusters,n_steps
    double precision,dimension(n_steps,n_steps),intent(in)::frame2frame
    double precision,intent(in)::cutoff_clcoeff
    integer,dimension(n_clusters),intent(in)::cluster_size
    integer,allocatable,dimension(:,:),intent(in)::cluster_members
    ! Although its intent is IN, the second dimension of cluster_members is not readily available (it is the possible maximum of cluster_size (before the actual
    ! size of the clusters were actaully computed)) and therefore I CHOOSE to make it allocatable although it will not be modified nor should it be
    double precision,dimension(n_clusters),intent(out)::clustering_coefficients

    integer::i,j,k,acs,jj,kk
  
    clustering_coefficients=0.d0
    if(cutoff_clcoeff==0.d0)return
    do i=1,n_clusters
      acs=cluster_size(i)
      if(acs.gt.1) then
        do j=2,acs-1
          jj=cluster_members(i,j)
          do k=j+1,acs
            kk=cluster_members(i,k)
            if( (frame2frame(jj,kk).gt.0.d0 .and. frame2frame(jj,kk).lt.cutoff_clcoeff) .or. (frame2frame(kk,jj).gt.0.d0 .and. frame2frame(kk,jj).lt.cutoff_clcoeff) ) then
              clustering_coefficients(i)=clustering_coefficients(i)+1.d0
            endif
          enddo
        enddo
        clustering_coefficients(i)=clustering_coefficients(i)/(dble(acs*(acs-1)/2))
      else
         clustering_coefficients(i)=1.d0
       endif
    enddo
end subroutine compute_clustering_coefficients
!------------------------------------------------------------------------------------------------!


!-----------------------------------
! Initiate seed for random number
!-----------------------------------
subroutine init_random_seed

  !------------------
  integer i, clock
  integer seed (100)
  !------------------

  !-----
  n=100
  !-----

  !------------------------------
  call random_seed(size = n)
  call system_clock(count=clock)
  !------------------------------

  !---------------------------------
  do i=1,100
     seed (i) = clock + 37 * (i-1)
  enddo
  !---------------------------------

  !----------------------------
  call random_seed(put = seed)
  !----------------------------
  
end subroutine init_random_seed
!-----------------------------------

end module clustering_mod

