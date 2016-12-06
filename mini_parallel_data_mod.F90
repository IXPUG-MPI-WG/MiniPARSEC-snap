! proper header
! description

module mini_parallel_data_mod
    use mini_constants

    implicit none

    type mini_parallel_data !{{{
        ! This is mostly copy paste from PARSEC
        ! Irrelevant parts are supposed to be removed,
        ! for example, groups.

        ! verbosity level (0,1,2)
        integer :: verbosity
        ! total number of processors
        integer :: procs_num
        ! ID of master
        integer :: masterid
        ! true for master processor, false for all other ones
        logical :: iammaster
        ! rank of all the processors in the groups
        integer :: iam
        ! communicator for inter-processor communication
        integer :: comm

        ! handle of the world group (containing all processors)
        integer :: world_handle
        ! total number of groups
        integer :: groups_num
        ! number of processors per group ( group_size * groups_num = procs_num)
        integer :: group_size
        ! order of group to which each processor belongs to
        integer :: mygroup
        ! true for the master processors of this group (each group has one
        ! and only one master)
        logical :: iamgmaster
        ! rank of the master processor of this group
        integer :: group_master
        ! rank of processors within a group
        integer :: group_iam
        ! handle of a group
        integer :: group_handle
        ! communicator within a group
        integer :: group_comm
        ! optimized comm for matvec
        integer :: group_comm_topo
        ! plan of the groups: all processors belonging to j-th group
        ! have parallel%iam = parallel%gmap(parallel%group_iam,j)
        integer, dimension(:,:), pointer :: gmap
        ! rank of processors in the group of masters
        integer :: gmaster_iam
        ! handle of the group of masters
        integer :: gmaster_handle
        ! communicator in the group of masters
        integer :: gmaster_comm
        ! size of workvec without +1 element for external elements
        integer :: workvec_edge

        ! actual size of hamiltonian
        integer ndim, nwedge
        ! local size of hamiltonian (equal to number of grid points
        ! held by each processor)
        integer mydim
        ! dimension of distributed grid arrays; in order to account
        ! for points outside the domain, choose ldn > max(numb)
        integer ldn

        ! Position of the first block of rows (=points in irreducible
        ! wedge) that each processor holds. The number of rows in
        ! processor ipe is irows(ipe+1)-irows(ipe). Global
        integer, dimension (:), pointer :: irows
        ! arrays for inter-processor communication
        ! they are all local (i.e., different values on different procs)
        ! see more about the meaning of these variables in comm_neigh.F
        integer, dimension (:), pointer :: irecvp, jsendp
        integer, dimension (:), pointer :: senrows
        integer, dimension (:), pointer :: pint
        integer :: inter1,inter2
        ! neighbor points to each point belonging to a computing PE;
        ! this index array indicates the position, in the 1-d array, 
        ! of the ith neighbor of point j (i from 1 to norder*3, because
        ! there are three directions, j from 1 to ndim). used for
        ! calculating derivatives; notice that neibs points to the
        ! index in irreducible wedge, not full grid!
        integer, dimension (:,:), pointer :: neibs
        ! local size for communication buffers (edges,size of send,size of receive)
        integer :: countcomm, maxcomm, maxcomm2
        ! local information about communication:
        integer, allocatable, dimension (:) :: &
        sources,destinations,sendcounts,recvcounts,rdispls,sdispls
        ! logical place to put blksize since solver structure does not exist
        integer :: blksize
        ! hardcoded width of wavefunction array
        integer :: num_wvfunc
        ! hardcoded number of outer loop benchmark iteration
        integer :: num_outer

    contains
        procedure :: init => mini_parallel_data_init
        procedure :: destroy => mini_parallel_data_destroy
        procedure :: comm_init => mini_parallel_data_setup
        procedure :: comm_neib => mini_comm_neib
        procedure :: topo_aware => mini_parallel_topo_aware
        procedure :: itac_init => init_trace_data

    end type mini_parallel_data !}}}

contains

    subroutine mini_parallel_data_init(mini_parallel)!{{{
        implicit none

        class(mini_parallel_data) :: mini_parallel

        nullify (mini_parallel%irows)
        nullify (mini_parallel%irecvp)
        nullify (mini_parallel%jsendp)
        nullify (mini_parallel%senrows)
        nullify (mini_parallel%pint)
        nullify (mini_parallel%neibs)
        nullify (mini_parallel%gmap)

        mini_parallel%num_wvfunc = 0


    end subroutine mini_parallel_data_init!}}}

    subroutine mini_parallel_data_destroy(mini_parallel)!{{{
        implicit none

        class(mini_parallel_data) :: mini_parallel

        if (associated (mini_parallel%irows))   deallocate (mini_parallel%irows)
        if (associated (mini_parallel%irecvp))  deallocate (mini_parallel%irecvp)
        if (associated (mini_parallel%jsendp))  deallocate (mini_parallel%jsendp)
        if (associated (mini_parallel%senrows)) deallocate (mini_parallel%senrows)
        if (associated (mini_parallel%pint))    deallocate (mini_parallel%pint)
        if (associated (mini_parallel%neibs))   deallocate (mini_parallel%neibs)
        if (associated (mini_parallel%gmap))    deallocate (mini_parallel%gmap)

    end subroutine mini_parallel_data_destroy!}}}

    subroutine mini_parallel_data_setup(mini_parallel,ierr)!{{{
        use mpi
        implicit none

        class(mini_parallel_data) :: mini_parallel
        integer, intent(inout) :: ierr

        integer mpinfo
        integer thread_provided
        integer thread_requested
#ifdef OMP
        thread_requested = MPI_THREAD_FUNNELED
#else
        thread_requested = MPI_THREAD_SINGLE
#endif
        ! Initialise MPI, get size and my id, create communicator
        call MPI_INIT_THREAD(thread_requested,thread_provided,mpinfo)

        if (thread_provided < thread_requested) then
            write(*,*) 'Error: MPI init does not support requested thread level '
            call MPI_ABORT(MPI_COMM_WORLD,-1,mpinfo)
        endif
        ierr = 0

        ! Save comm world
        mini_parallel%comm = MPI_COMM_WORLD
        ! Determine the size of this group. This group should have total
        ! number of PEs as its size
        call MPI_COMM_SIZE(mini_parallel%comm, mini_parallel%procs_num, mpinfo)
        ! Among this group define the PE with rank 0 as the master PE
        mini_parallel%masterid = 0
        ! Now determine the ranks
        call MPI_COMM_RANK(mini_parallel%comm, mini_parallel%iam, mpinfo)
        if (mini_parallel%masterid == mini_parallel%iam) then
            mini_parallel%iammaster = .true.
        else
            mini_parallel%iammaster = .false.
        endif
        ! until the group stuff is removed
        mini_parallel%group_size = mini_parallel%procs_num

    end subroutine mini_parallel_data_setup!}}}

    subroutine mini_parallel_topo_aware(mini_parallel) !{{{
        use mpi
        implicit none
        class(mini_parallel_data) :: mini_parallel
        integer :: graph_count, nnodes, inode, ipack, irank0,irank1
        integer, dimension(mini_parallel%group_size) :: nelem_in,nelem_out
        !type(MPI_Info) :: info !supposedly has some hints on the graph.
        integer :: info
        logical :: reorder   !wether to reorder ranks in new comm
        integer :: mpierr    !normal fortran mpi error int
        !
        !type(MPI_Comm) :: comm_dist_graph !the new communicator. 
        integer :: comm_dist_graph, comm_dist_graph_reorder !the new communicator. 
        !
        !   Topology data
        !
        integer :: graph_degree !size of arrays needed to describe local communication
        !   Actual local communication data
        integer :: max_in_degree,max_out_degree
        integer :: sources(0:mini_parallel%countcomm-1)
        integer :: source_weights(0:mini_parallel%countcomm-1)
        integer :: destinations(0:mini_parallel%countcomm-1)
        integer :: destination_weights(0:mini_parallel%countcomm-1)
        ! integer arrays with the j-th element containing the number
        ! of elements to send to neighbor j
        integer :: sendcounts(0:mini_parallel%countcomm-1)
        integer :: recvcounts(0:mini_parallel%countcomm-1)
        ! integer arrays with the j-th element containing the displacement to
        ! send/recvbuff for the j-th neighbor data
        integer :: sdispls(0:mini_parallel%countcomm-1)
        integer :: rdispls(0:mini_parallel%countcomm-1)
        !
        ! verification
        !
        integer :: source_check
        integer :: destination_check
        logical :: weight_check
        logical :: lopen
        integer :: writeto

        ! TODO : test if you can actually call topo right now or is it too early

        ! is each pe writing to separate file?
        if (mini_parallel%verbosity == DEBUGEACH) then
        lopen = .TRUE.
        writeto = 9
        endif


        graph_degree = mini_parallel%countcomm
        nnodes = mini_parallel%group_size

        do inode = 0, nnodes-1
            nelem_in(inode+1) = mini_parallel%jsendp(inode+1)-mini_parallel%jsendp(inode)
            nelem_out(inode+1) = mini_parallel%irecvp(inode+1)-mini_parallel%irecvp(inode)
        enddo

        graph_count=0
        do inode = 0,nnodes-1
            if (nelem_in(inode+1) /= 0) then
                sources(graph_count)=inode
                source_weights(graph_count)=nelem_in(inode+1)
                graph_count= graph_count + 1
            endif
        enddo

        graph_count=0
        do inode = 0,nnodes-1
            if (nelem_out(inode+1) /= 0) then
                destinations(graph_count)=inode
                destination_weights(graph_count)=nelem_out(inode+1)
                graph_count= graph_count + 1
            endif
        enddo

        reorder = .FALSE.
        ! call MPI_INFO_CREATE(info, mpierr)
        info=MPI_INFO_NULL
        call MPI_DIST_GRAPH_CREATE_ADJACENT(mini_parallel%group_comm,graph_degree,sources,source_weights,&
            graph_degree,destinations,destination_weights, info,reorder, comm_dist_graph, mpierr )
        call MPI_DIST_GRAPH_NEIGHBORS_COUNT(comm_dist_graph,source_check,destination_check,weight_check, mpierr)

        if (lopen) then
            write(writeto,*) ' NOTE: NEIGHBORHOOD GRAPH COMMUNICATOR WAS CREATED!'
            write(writeto,*) 'topo_aware: I am now aware that I communicate with',source_check, 'PEs'
            if (source_check-destination_check /= 0) then
                write(writeto,*) 'how come you are sending and receiving from different number of PEs?'
                write(writeto,*) 'This is VERY WRONG.'
                ! oh if we could just throw an exception here
            endif
            write(writeto,*) ''
            !call MPI_TOPO_TEST(mini_parallel%group_comm,inode,mpinfo)
            !write(9,*) 'topo_test returns:',inode
        endif

        !actually change stuff for communication:
        !we call the graph_neighbors to make sure mpi got the message about the communication pattern
        mini_parallel%group_comm_topo = comm_dist_graph
        max_in_degree=mini_parallel%countcomm
        max_out_degree=mini_parallel%countcomm
        call MPI_DIST_GRAPH_NEIGHBORS(mini_parallel%group_comm_topo, max_in_degree, sources, source_weights,&
            max_out_degree, destinations, destination_weights, mpierr)

        do inode = 0,max_in_degree-1
            sendcounts(inode) = mini_parallel%jsendp(sources(inode)+1)-mini_parallel%jsendp(sources(inode))
        enddo
        sdispls(0)=0
        do ipack=1,max_in_degree-1
            sdispls(ipack)=sdispls(ipack-1)+sendcounts(ipack-1)
        enddo

        do inode = 0,max_out_degree-1
            recvcounts(inode) = mini_parallel%irecvp(destinations(inode)+1)-mini_parallel%irecvp(destinations(inode))
        enddo

        rdispls(0)=0
        do ipack=1,max_out_degree-1
            rdispls(ipack)=rdispls(ipack-1)+recvcounts(ipack-1)
        enddo

        allocate(mini_parallel%sendcounts(0:mini_parallel%countcomm-1))
        allocate(mini_parallel%sources(0:mini_parallel%countcomm-1))
        allocate(mini_parallel%sdispls(0:mini_parallel%countcomm-1))
        !allocate(mini_parallel%source_weights(0:mini_parallel%countcomm-1))
        mini_parallel%sendcounts=sendcounts
        mini_parallel%sources=sources
        mini_parallel%sdispls=sdispls
        !mini_parallel%source_weights=source_weights

        allocate(mini_parallel%recvcounts(0:mini_parallel%countcomm-1))
        allocate(mini_parallel%destinations(0:mini_parallel%countcomm-1))
        allocate(mini_parallel%rdispls(0:mini_parallel%countcomm-1))
        !allocate(mini_parallel%destination_weights(0:mini_parallel%countcomm-1))

        mini_parallel%recvcounts=recvcounts
        mini_parallel%destinations=destinations
        mini_parallel%rdispls=rdispls
        !mini_parallel%destination_weights=destination_weights

        ! tests:
        call MPI_DIST_GRAPH_CREATE_ADJACENT(mini_parallel%group_comm,graph_degree,sources,source_weights,&
            graph_degree,destinations,destination_weights, &
            info,.TRUE., comm_dist_graph_reorder, mpierr )
        call MPI_COMM_RANK(comm_dist_graph, irank0, mpierr)
        call MPI_COMM_RANK(comm_dist_graph_reorder, irank1, mpierr)
        call MPI_COMM_FREE(comm_dist_graph_reorder,mpierr)
        if (lopen) then
            write(writeto,*) ' NEIGHBORHOOD COMMUNICATION DATA UPDATED FOR EACH PE '
            if (irank0 /= irank1) then
                write(writeto,*) 'topo_aware: I am PE rank',irank0,'but MPI would like me to be rank',irank1,'!!'
            else
                write(writeto,*) 'topo_aware: I am PE rank',irank0,'both in the group comm and the topo comm'
            endif
        endif
    end subroutine mini_parallel_topo_aware !}}}

    !legacy subroutine, to be removed
    subroutine create_group_layout(mini_parallel,groups_num_in) !{{{
        use mpi
        implicit none
        type (mini_parallel_data), intent (inout) :: mini_parallel
        integer, intent(in) :: groups_num_in

        mini_parallel%group_iam = mini_parallel%iam
        mini_parallel%group_comm = mini_parallel%comm
        mini_parallel%mygroup = 1

        mini_parallel%group_master = 0
        if (mini_parallel%group_iam == mini_parallel%group_master) then
            mini_parallel%iamgmaster = .true.
        else
            mini_parallel%iamgmaster = .false.
        endif

    end subroutine create_group_layout !}}}

    subroutine psum(vec,nmax,nnodes,comm) !{{{

        ! the PARSEC allreduce wrapper for dp 
        use mini_constants
        use mpi
        implicit none
        ! dimension of array vec
        integer, intent(in) :: nmax
        ! number of procs available
        integer, intent(in) :: nnodes
        ! communicator
        integer, intent(in) :: comm
        ! The original vector whose value is partial for each PE
        real(dp), intent(inout) :: vec(nmax)

        integer mpinfo

        if (nnodes == 1) return
        if (nmax < 1) return
        call MPI_allREDUCE(MPI_IN_PLACE,vec,nmax,MPI_DOUBLE_PRECISION,MPI_SUM,comm,mpinfo)
    end subroutine psum !}}}

    subroutine int_allreduce(ivec,nmax,nnodes,comm) !{{{
        use mini_constants
        use mpi
        implicit none
        !
        !  Input/Output variables:
        !
        !  dimension of array vec
        integer, intent(in) :: nmax
        !  number of procs available
        integer, intent(in) :: nnodes
        !  communicator
        integer, intent(in) :: comm
        !  The original vector whose value is partial for each PE
        integer, intent(inout) :: ivec(nmax)
        !  Work variables:
        !  exit code for mpi calls
        integer mpinfo
        !  ---------------------------------------------------------------
        if (nnodes == 1) return
        call MPI_allREDUCE(MPI_IN_PLACE,ivec,nmax,MPI_INTEGER,MPI_SUM,comm,mpinfo)
    end subroutine int_allreduce !}}

    subroutine mini_comm_neib(mini_parallel,norder)
        !  Rewrites parallel%neibs with a local one with the properties:
        !
        !              |--> parallel%neibs(neigh,row)-ioffset for all 
        !              |     neighbors "neigh" local on this processor
        !              |
        !  parallel%neibs(neigh,row)= --> position of neigh, (> mydim) in the boundary 
        !              |    values appended at the end of the current vector.
        !              |    This is the "neigh" neighbor of the local "row" 
        !              |   
        !              |--> ndim+1 for all unused neighbors 
        use mini_constants
        use mpi
        implicit none
        !
        !  Input/Output variables:
        !
        !  parallel computation related data
        class(mini_parallel_data) :: mini_parallel
        !  number of neighbors used on one side in numerical derivative
        integer, intent(in) :: norder
        !
        !  Work variables:
        !
        integer nnodes,iam,node,mpinfo
        integer ielem,jelem,istart,jstart,irow,myrow,itemp,mydim &
            ,ioffset,ncur,ish,their_row

        logical, allocatable :: sentto(:,:)
        integer, allocatable :: received(:)

        integer, dimension(:), allocatable :: &
            irows,isendcount,irecvcount

        integer, allocatable :: recrows(:)

        integer :: neibs_num

        integer inter1, inter2

        integer i, alcstat , ierr

        !
        !  External functions:
        !
!        integer, external :: lochome

        !---------------------------------------------------------------
        iam = mini_parallel%group_iam
        nnodes = mini_parallel%group_size
        mydim   = mini_parallel%mydim

        allocate(irows(0:nnodes))
        allocate(sentto(mini_parallel%mydim,0:nnodes))
        allocate(isendcount(0:nnodes))
        allocate(irecvcount(0:nnodes))
        allocate(received(mini_parallel%nwedge))
        allocate(recrows(mini_parallel%nwedge))

        irows(:) = mini_parallel%irows
        ioffset = irows(iam) - 1

        neibs_num = 6 

        received(:) = 0
        recrows(:) = 0
        mini_parallel%senrows(:) = 0
        mini_parallel%irecvp(:) = 0
        mini_parallel%jsendp(:) = 0
        isendcount(:) = 0
        irecvcount(:) = 0

        sentto(:,:) = .false.
        inter1 = 0

        ncur = 0
        !  for all neighbors of this row: 
        !  2*(3+lap_dir_num)*norder (first shell first)

        do ish = 1, neibs_num*norder
        do myrow = 1, mydim
        irow = mini_parallel%neibs(ish,myrow)
        !  find PE # of neighbor irow
        node = mini_lochome(irow,nnodes,irows,iam)

        if (node /= iam) then
            if (.not.sentto(myrow,node)) then
                mini_parallel%jsendp(node) = mini_parallel%jsendp(node) + 1
                ! AJB: ish=6 is last element in the first shell
                ! We don't use this kind of thing around here no more
                !if (ish < 7) parallel%jp1(node) = parallel%jp1(node)+1
                sentto(myrow,node) = .true.
            endif
            if (received(irow) ==  0) then
                mini_parallel%irecvp(node) = mini_parallel%irecvp(node) + 1
                ! AJB: ish=6 is last element in the first shell
                ! We don't use this kind of thing around here no more
                !if (ish < 7) parallel%ip1(node) = parallel%ip1(node)+1
                ncur = ncur + 1
                mini_parallel%senrows(ncur) = irow
                ! AJB: it seems that we add 1 so that node 0 wont screw up the loop?
                received(irow) = node+1
            endif
        endif
        enddo                  ! myrow = 1, mydim
        enddo                     ! ish = 1,2*(3+parallel%lap_dir_num)*norder


        !
        !  Find where info from (irecvp) and to (jsendp) each proc starts
        !
        ielem = mini_parallel%irecvp(0)
        jelem = mini_parallel%jsendp(0)
        mini_parallel%irecvp(0) = 1
        mini_parallel%jsendp(0) = 1
        do node = 1, nnodes
        istart = mini_parallel%irecvp(node-1) + ielem
        jstart = mini_parallel%jsendp(node-1) + jelem
        ielem = mini_parallel%irecvp(node)
        jelem = mini_parallel%jsendp(node)
        mini_parallel%irecvp(node) = istart
        mini_parallel%jsendp(node) = jstart
        enddo
        !
        !  Use received and senrows arrays to paste together the columns 
        !  I will need. Reset received to zero wherever it has something
        !  else. Rows are according to the global permutation.
        !
        do i = 1, ncur
        irow = mini_parallel%senrows(i)
        node = received(irow)-1
        received(irow) = 0
        recrows( mini_parallel%irecvp(node) ) = mini_parallel%senrows(i)
        mini_parallel%irecvp(node) = mini_parallel%irecvp(node) + 1
        enddo
        !
        !  Restore the ip index to starting points in the appended vector
        !  for each node. This is local now (eg irecvp(0) = 1 + mydim)
        !
        do i = nnodes-1,1,-1
        mini_parallel%irecvp(i) = mini_parallel%irecvp(i-1) + mydim
        enddo
        mini_parallel%irecvp(0) = 1 + mydim

        !  Go through mini_parallel%neibs array to reset the neigh global
        !  indexing to a local. Notice the order of the loops must be the
        !  same as in previous loops over ish, myrow.
        !  IMPORTANT: Since this is not 1-to-1 mapping, the same row may
        !  be accessed as a neighbor of many rows in mini_parallel%neibs. To
        !  remember where this row resides in the new order, the received
        !  array now records the position of the first time this row
        !  (neighbor) is encountered.
        !
        do ish = 1,neibs_num*norder
        !  for this neighbor of all rows
        do myrow = 1, mydim
        irow = mini_parallel%neibs(ish,myrow)
        node = mini_lochome(irow,nnodes,irows,iam)
        if (node /= iam) then
            if (received(irow) ==  0) then
                ! First time encountered. Change irecvp, received
                received(irow)  = mini_parallel%irecvp(node)
                ! if(mini_parallel%irecvp(node) > mini_parallel%irecvp(nnodes)+mydim) then
                !     write(9,*) 'WARNING?! in comm_neigh:'
                !     write(9,*) 'ipn>> ',ish,myrow,irow,mini_parallel%irecvp(node),mini_parallel%irecvp(nnodes)
                ! endif

                mini_parallel%irecvp(node) = mini_parallel%irecvp(node) + 1
            else
                ! This row exists already in the list
                if(received(irow) > mini_parallel%irecvp(nnodes)+mydim .or. received(irow) < 1) then 
                    ! ! AJB: this seems to be a problem too!
                    ! write(9,*) 'WARNING?! in comm_neigh:'
                    ! write(9,*) 'exs>> ',ish ,myrow,irow,received(irow)
                endif
            endif !received==0
            mini_parallel%neibs(ish,myrow) = received(irow)

        else
            ! AJB: prune the tree towards nwedge+1
            if (irow == mini_parallel%nwedge+1) cycle
            mini_parallel%neibs(ish,myrow) = irow - ioffset 
        endif
        enddo                  ! myrow = 1, mydim
        enddo

        !
        !  Here would be the place to find and permute the inter1 interior points in pint
        ! AJB: pint is not a permutation on this branch 
        !
        do i = 1, mydim
        mini_parallel%pint(i) = i
        enddo

        do myrow = 1,mydim
            do ish = 1,neibs_num*norder
            irow = mini_parallel%neibs(ish,myrow)
            if (irow > mydim) goto 10 
            enddo
            inter1 = inter1+1      ! interior node 
            itemp  = mini_parallel%pint(myrow)   ! permute to the begining !AJB ...maybe not...
            mini_parallel%pint(myrow) = mini_parallel%pint(inter1)
            mini_parallel%pint(inter1)= itemp
10  continue
        enddo
        inter2=0

        !
        !  Restore again the irecvp index to starting points, but in the 
        !  non-appended vector to be used in the alltoall
        !
        do i = nnodes-1,1,-1
        mini_parallel%irecvp(i) = mini_parallel%irecvp(i-1)-mydim
        enddo
        mini_parallel%irecvp(0) = 1 
        !
        !  Communicate the information from other processors' irecvp into 
        !  the appropriate jsendp locations. This is an alltoallv step 
        !  where the j-th section of recrows (from irecvp(j) to irecvp(j+1)-1) on 
        !  processor i, goes to the j-th processor senrows between its 
        !  jsendp(i) and jsendp(i+1)-1.
        !
        do i=0,nnodes-1
        isendcount(i) = mini_parallel%irecvp(i+1)-mini_parallel%irecvp(i)
        irecvcount(i) = mini_parallel%jsendp(i+1)-mini_parallel%jsendp(i)
        enddo
        do i=0,nnodes
        mini_parallel%irecvp(i) = mini_parallel%irecvp(i)-1
        mini_parallel%jsendp(i) = mini_parallel%jsendp(i)-1
        enddo

        call MPI_Alltoallv(recrows,isendcount, mini_parallel%irecvp, MPI_INTEGER, &
            mini_parallel%senrows,irecvcount, mini_parallel%jsendp, MPI_INTEGER, &
            mini_parallel%group_comm,mpinfo)

        !
        !  Finally, restore irecvp and jsendp indices 
        !
        do i = 0, nnodes
        mini_parallel%irecvp(i) = mini_parallel%irecvp(i) + 1 + mydim
        mini_parallel%jsendp(i) = mini_parallel%jsendp(i) + 1
        enddo

        !  printout
#ifdef DEBUG
        write(9,*) 'INTERIOR POINTS Mvec, PreCond (', inter1,inter2, ') '
        if (inter1 == 0) then
         write(9,*) 'WARNING, NO MATVEC INTERIOR POINTS'
         write(9,*) 'CURRENT SYSTEM SIZE/PE RATIO NOT RECOMMENDED'
        endif
#endif
        mini_parallel%inter1=inter1
        mini_parallel%inter2=inter2
        !  Change the senrows to local numbering (subtract the offset)
        do i = 1, mini_parallel%jsendp(nnodes)-1
        mini_parallel%senrows(i) = mini_parallel%senrows(i) - ioffset
        enddo

        mini_parallel%workvec_edge = mydim + mini_parallel%irecvp(nnodes)-mini_parallel%irecvp(0)
        ! Experimental - change the neib structure so that instead of nwedge+1
        ! the external points are at mydim+ irecvp(nnodes)-irecvp(0) + 1
        do myrow = 1,mydim
            do ish = 1,neibs_num*norder
            irow = mini_parallel%neibs(ish,myrow)
            if (irow > mini_parallel%nwedge) then
                mini_parallel%neibs(ish,myrow) = mini_parallel%workvec_edge + 1
            endif
            enddo
        enddo




        deallocate(received)
        deallocate(recrows)
        deallocate(irows)
        deallocate(sentto)
        deallocate(isendcount)
        deallocate(irecvcount)
    end subroutine mini_comm_neib

integer function mini_lochome(irow,nnodes,irows,iam)

  implicit none
  !
  !  Input/Output variables:
  !
  integer, intent(in) :: irow, iam, nnodes
  integer, intent(in) :: irows(0:nnodes)
  !
  !  Work variables:
  !
  integer i
  !---------------------------------------------------------------
  mini_lochome = iam
  do i = 0, nnodes-1
     if (irows(i) <= irow .and. irow < irows(i+1)) then
        mini_lochome = i
        exit
     endif
  enddo

end function mini_lochome

!=============================Intel TraceAnalyzer API ====================================!{{{
            subroutine init_trace_data (mini_parallel,vtierr) 
                implicit none
#ifdef ITAC
                include 'VT.inc'
                include 'mini_vt.inc'
                class(mini_parallel_data) :: mini_parallel
                integer vtierr
                vtierr = 0
            call VTCLASSDEF( 'SETUP',     main_class ,vtierr )
            call VTCLASSDEF( 'STENCIL',     lap_class ,vtierr )
            call VTCLASSDEF( 'NOT-STENCIL',     non_lap_class ,vtierr )
            call VTCLASSDEF( 'BUFFERS',     buffering_class ,vtierr )
            call VTFUNCDEF( 'Buffering'     , buffering_class, vt_buffer , vtierr              )
            call VTFUNCDEF( 'Prebuffering'     , buffering_class, vt_prebuffer , vtierr              )
            call VTFUNCDEF( 'Diag'     , non_lap_class, vt_diag , vtierr              )
            call VTFUNCDEF( 'Mock-laplacian (inner)'     , lap_class, vt_mocklap_inner , vtierr              )
            call VTFUNCDEF( 'Mock-laplacian (outer)'     , lap_class, vt_mocklap_outer , vtierr              )
            call VTFUNCDEF( 'Outside Kernel'     , non_lap_class, vt_outside , vtierr              )
#else
                class(mini_parallel_data) :: mini_parallel
                integer, intent(inout):: vtierr
            vtierr = 0
#endif
            end subroutine init_trace_data 
!=============================Intel TraceAnalyzer API ====================================!}}}


end module mini_parallel_data_mod
