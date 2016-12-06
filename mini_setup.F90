! proper header
! description

subroutine mini_setup(ierr)
    use mini_constants
    use mini_common_mod
    use mini_global_data_mod
    ! do I need this as well?
    use mini_parallel_data_mod
    use mini_grid_mod
    use mini_buffers_mod
    ! end "do i need this as well"
    use mpi

    implicit none

    integer, intent(inout) :: ierr

    integer  irows(0:mini_parallel%procs_num)
    ! counters
    integer :: index_arr,neib,ii,j,k,indx(LAPDIR),numb !TODO - not 12, but 6!
    integer :: num_p,neibs_num,icount,nelem,mxelem,ishell
    ! shorthand notation
    integer :: nnodes,masterid,inode
    ! potentially misleading
    integer :: comm,nord2
    ! msg for neighbor exchange bookkeeping
    integer :: msgtype 
    ! for checking 
    integer :: alccheck, mpinfo
    ! for creating neib, probably not needed 
    integer, allocatable :: tmp2(:,:)
    integer :: mpi_stat(MPI_STATUS_SIZE)

    ! legacy support
    call create_group_layout(mini_parallel,1)

    nnodes   = mini_parallel%group_size
    masterid = mini_parallel%group_master
    comm     = mini_parallel%group_comm
    nord2    = mini_grid%norder
    msgtype  = 0

    !TODO: end with BCASTING hardcoded values as a proactive measure

    call MPI_BCAST(mini_grid%step          , 3,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
    call MPI_BCAST(mini_grid%rmax          , 1,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
    call MPI_BCAST(mini_grid%h_2           , 1,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
    call MPI_BCAST(mini_grid%hcub          , 1,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
    call MPI_BCAST(mini_grid%hcub2         , 1,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
    call MPI_BCAST(mini_grid%d_shape_param , 3,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
    call MPI_BCAST(mini_grid%i_shape_param , 1,MPI_INTEGER,masterid,comm,mpinfo)
    call MPI_BCAST(mini_grid%domain_shape  , 1,MPI_INTEGER,masterid,comm,mpinfo)
    call MPI_BCAST(mini_grid%norder        , 1,MPI_INTEGER,masterid,comm,mpinfo)
    call MPI_BCAST(mini_grid%shift, 3,MPI_DOUBLE_PRECISION, mini_parallel%masterid,mini_parallel%comm,mpinfo)

    call MPI_BCAST(mini_grid%n1,1,MPI_INTEGER,masterid,comm,mpinfo)
    call MPI_BCAST(mini_grid%n2,1,MPI_INTEGER,masterid,comm,mpinfo)
    call MPI_BCAST(mini_grid%n3,1,MPI_INTEGER,masterid,comm,mpinfo)
    !
    !  Create communicator arrays.
    !
    allocate(mini_parallel%irows(0:nnodes))

    allocate(mini_parallel%irecvp(0:nnodes))
    allocate(mini_parallel%jsendp(0:nnodes))

    ! this is a simplified form of the parsec "super" subroutine
    ! sets up the domain decomposition etc.
    ! master only
    ierr = 0
    if (mini_parallel%iamgmaster) then
    call mini_grid_partition(ierr)
    endif
    ! collective call
    call mini_exit_err(ierr)

    ! simple BCAST to share the results of partition
    ! some data duplication here but at least it works ;)
    call MPI_BCAST(mini_parallel%ndim,   1, MPI_INTEGER, masterid,comm,mpinfo)
    call MPI_BCAST(mini_parallel%nwedge, 1, MPI_INTEGER, masterid,comm,mpinfo)
    call MPI_BCAST(mini_parallel%ldn,    1, MPI_INTEGER, masterid,comm,mpinfo)
    call MPI_BCAST(mini_parallel%irows, nnodes+1, MPI_INTEGER,masterid, comm,mpinfo)
    call MPI_BCAST(mini_grid%nxyz, 1, MPI_INTEGER, masterid, comm, mpinfo)

    irows(0:nnodes) = mini_parallel%irows

    
    neibs_num=LAPDIR


    allocate(mini_parallel%neibs(neibs_num*mini_grid%norder,mini_parallel%ldn), stat=alccheck)
    ierr = 0 
    if (alccheck /= 0) then 
        ierr = 1000
        ! ERROR CODE 1000 :: MEMORY ALLOCATION FAILURE FOR NEIBS
    else
        mini_parallel%neibs = 0
    endif
    !call mini_exit_err(ierr) !avod sync here
    allocate(mini_parallel%pint(mini_parallel%ldn), stat=alccheck)
    ierr = 0 
    if (alccheck /= 0) then 
        ierr = 1001
        ! ERROR CODE 1001 :: MEMORY ALLOCATION FAILURE FOR PINT OR SENROWS
    endif

    allocate(mini_parallel%senrows(mini_parallel%nwedge), stat=alccheck)
    ierr = 0 
    if (alccheck /= 0) then 
        ierr = 1001
        ! ERROR CODE 1001 :: MEMORY ALLOCATION FAILURE FOR PINT OR SENROWS
    endif

    if (.not. mini_parallel%iamgmaster) then
        mini_grid%nwedge = mini_parallel%nwedge
        call mini_grid%set_wedge(nnodes,ierr)
        call mini_grid%set_ist(nnodes)
    endif
    call mini_exit_err(ierr) !test every allocation failure from before

    call MPI_BCAST(mini_grid%kx, mini_grid%nwedge, MPI_INTEGER, masterid,comm,mpinfo)
    call MPI_BCAST(mini_grid%ky, mini_grid%nwedge, MPI_INTEGER, masterid,comm,mpinfo)
    call MPI_BCAST(mini_grid%kz, mini_grid%nwedge, MPI_INTEGER, masterid,comm,mpinfo)

    !  Creating Neibs!!
    !  For each point in the sphere, make an array which points to the
    !  position of all neighbors to be used for derivation. There are
    !  3*norder of those for each point.
    !
    allocate(tmp2(neibs_num*mini_grid%norder,mini_parallel%ldn),stat=alccheck)
    ierr = 0 
    if (alccheck /= 0) then 
        ierr = 1002
        ! ERROR CODE 1002 :: MEMORY ALLOCATION FAILURE FOR tmp
    endif
    call mini_exit_err(ierr)

    index_arr = 0
   do inode = 0,nnodes-1
    !  In parallel environment, send lists of neighbor points to the
    !  PEs that are going to handle them.
    msgtype = msgtype + 1
    num_p = irows(inode+1) - irows(inode)
    tmp2 = 0
    if (mini_parallel%iamgmaster) then
      ii = 0
      do numb = irows(inode), irows(inode+1)-1
        index_arr = index_arr+1
        ii = ii + 1
        neib = 0
        do ishell = 1,nord2

        !  Doing first the main axes neighbors

        indx(1) = mini_grid%indexg(mini_grid%kx(numb)-ishell,mini_grid%ky(numb)       ,mini_grid%kz(numb))
        indx(2) = mini_grid%indexg(mini_grid%kx(numb)+ishell,mini_grid%ky(numb)       ,mini_grid%kz(numb))
        indx(3) = mini_grid%indexg(mini_grid%kx(numb)       ,mini_grid%ky(numb)-ishell,mini_grid%kz(numb))
        indx(4) = mini_grid%indexg(mini_grid%kx(numb)       ,mini_grid%ky(numb)+ishell,mini_grid%kz(numb))
        indx(5) = mini_grid%indexg(mini_grid%kx(numb)       ,mini_grid%ky(numb)       ,mini_grid%kz(numb)-ishell)
        indx(6) = mini_grid%indexg(mini_grid%kx(numb)       ,mini_grid%ky(numb)       ,mini_grid%kz(numb)+ishell)
        do j = 1, LAPDIR
        tmp2(neib+j,ii) = mini_grid%rindex(indx(j))
        enddo

        neib = neib + neibs_num
        enddo
      enddo

        if (mini_parallel%group_iam == inode) then
            mini_parallel%neibs = tmp2(:,:)
        else
            call MPI_SEND(tmp2, neibs_num*mini_grid%norder*mini_parallel%ldn, MPI_INTEGER, inode,msgtype, comm, mpinfo)
#ifdef DEBUG
             write(9,*) 'neibs sent to node ',inode
#endif
        endif
    else
        if (mini_parallel%group_iam == inode) then
            call MPI_RECV(tmp2,neibs_num*mini_grid%norder*mini_parallel%ldn,MPI_INTEGER,masterid,msgtype, comm,mpi_stat,mpinfo)
            mini_parallel%neibs = tmp2(:,:)
#ifdef DEBUG
            write(9,*) 'GOT DA NEIBS'
#endif
        endif
    endif
   enddo
    call MPI_Barrier(comm,mpinfo)
    deallocate(tmp2)

    !  Master PE broadcasts indexw
    call MPI_BCAST(mini_grid%nxmax,1,MPI_INTEGER,masterid,comm,mpinfo)
    call MPI_BCAST(mini_grid%nxmin,1,MPI_INTEGER,masterid,comm,mpinfo)
    call MPI_BCAST(mini_grid%nymax,1,MPI_INTEGER,masterid,comm,mpinfo)
    call MPI_BCAST(mini_grid%nymin,1,MPI_INTEGER,masterid,comm,mpinfo)
    call MPI_BCAST(mini_grid%nzmax,1,MPI_INTEGER,masterid,comm,mpinfo)
    call MPI_BCAST(mini_grid%nzmin,1,MPI_INTEGER,masterid,comm,mpinfo)

    nelem = (mini_grid%nxmax-mini_grid%nxmin+1)*(mini_grid%nymax-mini_grid%nymin+1)* (mini_grid%nzmax-mini_grid%nzmin+1)
    ierr = 0 
    if (.not. mini_parallel%iamgmaster) then
      allocate (mini_grid%indexw (mini_grid%nxmin:mini_grid%nxmax, mini_grid%nymin:mini_grid%nymax,mini_grid%nzmin:mini_grid%nzmax) ,stat=alccheck)
        if (alccheck /= 0) then 
            ierr = 1003
            ! ERROR CODE 1003 :: Non master process could not allocate indexw
            ! indexw should not really be used in the benchmark BTW since no symmetry
        endif


    endif
    ! hack for checking if any process failed
    call mini_exit_err(ierr)

    call MPI_BCAST(mini_grid%ndim,1,MPI_INTEGER,masterid,comm,mpinfo)
    call MPI_BCAST(mini_grid%nwedge,1,MPI_INTEGER,masterid,comm,mpinfo)
    if (.not. mini_parallel%iamgmaster) then
        mini_grid%indexw = 0
    endif
    call MPI_BCAST(mini_grid%indexw,nelem,MPI_INTEGER,masterid,comm,mpinfo)

    if (.not. mini_parallel%iamgmaster) then
        call mini_grid%set_ndim(mini_grid%ndim,ierr)
    endif
    call mini_exit_err(ierr)

    call MPI_BCAST(mini_grid%rindex,mini_grid%ndim+1,MPI_INTEGER,masterid,comm,mpinfo)

    call MPI_BCAST(mini_grid%fx,mini_grid%ndim,MPI_INTEGER,masterid,comm,mpinfo)
    call MPI_BCAST(mini_grid%fy,mini_grid%ndim,MPI_INTEGER,masterid,comm,mpinfo)
    call MPI_BCAST(mini_grid%fz,mini_grid%ndim,MPI_INTEGER,masterid,comm,mpinfo)

    !  define local dimension mydim
    mini_parallel%mydim = irows(mini_parallel%group_iam+1) - irows(mini_parallel%group_iam)
    if (mini_parallel%verbosity == DEBUGEACH) then
    write(9,*) 'my group_iam', mini_parallel%group_iam
    write(9,*) 'The local dimension is :', mini_parallel%mydim
    endif

    !AJB: might as well remove this someday
    call MPI_Barrier(mini_parallel%comm,mpinfo)

  !  Find boundary information for exchange
    call mini_parallel%comm_neib(mini_grid%norder)

    if (mini_parallel%verbosity == DEBUGEACH) & write(9,*) 'Information to be sent to neighbors:'

    icount = 0
    mxelem = 0
    do inode = 0, nnodes - 1
        nelem = mini_parallel%jsendp(inode+1)-mini_parallel%jsendp(inode)
        if (mini_parallel%verbosity == DEBUGEACH) & write(9,*) 'send PE ',inode,' #rows:',nelem
        mxelem = mxelem + nelem
        if (nelem/=0) icount = icount +1
    enddo

    mini_parallel%maxcomm   = mxelem
    mini_parallel%countcomm = icount
    mini_parallel%maxcomm2  = mini_parallel%irecvp(nnodes)-mini_parallel%irecvp(0)

    if (mini_parallel%verbosity == DEBUGEACH) then
    write(9,*) 'I communicate with ',icount,' processors, a total of ',mini_parallel%maxcomm,' elements'
    write(9,*) 'Also, I receive from these procs: ', mini_parallel%maxcomm2, ' elements'

        if (mini_parallel%maxcomm /= mini_parallel%maxcomm2) then
            write(9,*) ''
            write(9,*) 'WARNING, ASSYMETRIC COMM PATTERN ENCOUNTERED'
            write(9,*) ''
        endif
    endif

    !  create the graph communicator and update local comm data
    !   when that's done
    call mini_parallel%topo_aware()

    do ii=1,mini_parallel%blksize !hardcoded 
#ifdef DEBUG
    write(9,*) 'setting up buffer number',ii
#endif

    call mini_buffers(ii)%setup(mini_parallel,ierr)

#ifdef DEBUG
    write(9,*) 'buffer status - protect = ', mini_buffers(ii)%protect
#endif
    end do
#ifdef DEBUG
                     ! write(*,*) 'DEBUG:: after buffer setup'
                     ! write(*,*) 'DEBUG:: I am', mini_parallel%iam
                     ! write(*,*) 'DEBUG:: rmax is now', mini_grid%rmax
                     ! write(*,*) 'DEBUG:: stepin is now', mini_grid%stepin
#endif

    if (ierr /=0) then
        ierr = ierr + 1000
        ! ERROR CODE 11Y0 :: BUFFER ALLOCATION ERROR
    endif
    call mini_exit_err(ierr)

    !what is the true length of workvec...
    !this limits the system size because each PE
    !has to hold a vector with the enitre system length
    !it is going to be changed
    !allocate(workvec(mini_parallel%nwedge+1),stat=alccheck)
    allocate(workvec(mini_parallel%workvec_edge+1),stat=alccheck)
    if (alccheck /= 0) then 
        ierr = 1004
        ! ERROR CODE 1004 :: MEMORY ALLOCATION FAILURE FOR GLOBAL WORKVEC
    endif

    allocate(mockpot(mini_parallel%ldn),stat=alccheck)
    if (alccheck /= 0) then 
        ierr = 1005
        ! ERROR CODE 1005 :: MEMORY ALLOCATION FAILURE FOR LOCAL MOCKPOT
    endif

    !TODO - this should not be hardcoded
    allocate(wvfuncs(mini_parallel%ldn,mini_parallel%num_wvfunc),stat=alccheck)
    if (alccheck /= 0) then 
        ierr = 1006
        ! ERROR CODE 1006 :: MEMORY ALLOCATION FAILURE FOR WVFUNCS
    endif

    call mini_exit_err(ierr)

    ! fill the wavefunctions with random numbers
    call my_random_array(wvfuncs,mini_parallel%ldn,mini_parallel%mydim,mini_parallel%num_wvfunc)
    ! fill the potential with random numbers
    call my_random_array(mockpot,mini_parallel%ldn,mini_parallel%mydim,1)

  !now we setup the laplacian by filling up the coefficents with random numbers
     call mini_grid%setup_mock_laplacian()

end subroutine mini_setup
