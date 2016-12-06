!--------------------------------------------------------------------------!
! proper header
!
!> Doxygen description
!! This module contains the buffer type and functions
!! Because the buffers are the one involved in the
!! ghost exchange, they are also responsible for the communication
!! using the topologically aware communicator
!
!> @todo:    consider MPI memory (de)allocation routines for faster operations?
!!           This would be done in whatever alloc is pointing at.
!--------------------------------------------------------------------------!
module mini_buffers_mod
    use mini_constants
    implicit none

    type mini_buffer !{{{
        !>send buffer 
        real(dp), dimension(:),pointer :: send
        !>receive buffer 
        real(dp), dimension(:),pointer :: receive
        !>crude write protect
        logical :: protect
        !
        integer :: vec_num
        integer :: request_number
#ifdef NOCOLLECTIVE
        integer, allocatable, dimension (:) :: request
#else
        integer :: request
#endif

    contains
        procedure :: init => init_buffer
        procedure :: destroy => destroy_buffer
        procedure :: setup => setup_buffer
        procedure :: alloc => allocate_buffer 
        procedure :: wprotect => buff_write_protect
        procedure :: wrelease => buff_write_release
        procedure :: packing => buff_pack
        procedure :: unpacking => buff_unpack
#ifdef BLOCKING
        procedure :: buffcomm => buff_comm_block
        procedure :: testing => buff_nothing
        procedure :: waiting => buff_nothing
#else
        procedure :: buffcomm => buff_comm_non_block
        procedure :: testing => buff_test_non_block
        procedure :: waiting => buff_wait_non_block
#endif

    end type mini_buffer !}}}


contains
    subroutine init_buffer(buffer) !{{{
        implicit none
        class (mini_buffer) :: buffer

        nullify (buffer%send) 
        nullify (buffer%receive)
        buffer%protect = .true.
        buffer%request_number=0
        buffer%vec_num=0
    end subroutine init_buffer !}}}


    subroutine destroy_buffer(buffer) !{{{
        use mpi
        implicit none
        class (mini_buffer) :: buffer
        integer::mpinfo

        if (buffer%protect) then
            write(*,*) 'CRITICAL!!! in: destroy_buffer - buffer write protected'
            !   write(9,*) "destroy_buffer critical error: my buffer was write protected, num:",buffer%vec_num
            !   call myflush(9)
            call MPI_ABORT(MPI_COMM_WORLD,-1,mpinfo) 
        endif

        if (associated(buffer%send)) deallocate(buffer%send)
        if (associated(buffer%receive)) deallocate(buffer%receive)
#ifdef NOCOLLECTIVE
        if (allocated(buffer%request)) deallocate(buffer%request)
#endif
    end subroutine destroy_buffer !}}}

    subroutine allocate_buffer(buffer,mini_parallel,ierr) !{{{
        use mini_parallel_data_mod
        implicit none
        class (mini_buffer) :: buffer
        type (mini_parallel_data), intent (in) :: mini_parallel
        integer, intent(inout) :: ierr
        integer :: alccheck

        allocate(buffer%send (mini_parallel%maxcomm),stat=alccheck)
        if (alccheck /=0) ierr=ierr+110
        allocate(buffer%receive (mini_parallel%maxcomm2),stat=alccheck)
        if (alccheck /=0) ierr=ierr+120
#ifdef NOCOLLECTIVE
        allocate(buffer%request (2*mini_parallel%group_size), stat=alccheck)
        if (alccheck /=0) ierr=ierr+130
#endif
        !now you can write to this buffer
        call buffer%wrelease()
        !buffer%protect = .false.

    end subroutine allocate_buffer !}}}

    subroutine setup_buffer(buffer,mini_parallel,ierr) !{{{
        use mini_parallel_data_mod
        use mpi, ONLY:MPI_REQUEST_NULL
        implicit none
        class (mini_buffer) :: buffer
        type (mini_parallel_data), intent (in) :: mini_parallel
        integer, intent(inout) :: ierr

        !basic protection against 1 PE jobs
        if (mini_parallel%group_size == 1) return

        call buffer%alloc(mini_parallel,ierr)

        buffer%request = MPI_REQUEST_NULL

    end subroutine setup_buffer !}}}

    subroutine buff_write_protect(buffer) !{{{
        implicit none
        class (mini_buffer) :: buffer
        buffer%protect = .TRUE.
    end subroutine buff_write_protect!}}}

    subroutine buff_write_release(buffer) !{{{
        implicit none
        class (mini_buffer) :: buffer
        buffer%protect = .FALSE.
    end subroutine buff_write_release !}}}

    subroutine buff_pack(buffer,mini_parallel,p,current,tobuff) !{{{
        use mini_constants
        use mini_parallel_data_mod
        use mpi
        implicit none
        !
        !  Input/Output variables:
        !
        !
        class (mini_buffer) :: buffer
        type (mini_parallel_data), intent(in) :: mini_parallel
        !
        integer, intent(in) :: current !the column of p that is going to be calculated on
        integer, intent(in) :: tobuff !the column of p that we are exchanging data for
        !
        real(dp), intent(in) :: p(mini_parallel%ldn,mini_parallel%num_wvfunc) ! p is p
        !
        !
        ! Work variables
        !
        integer :: inode,jst,jelem,j_comm,mpinfo

        !basic protection against 1 PE jobs
        if (mini_parallel%group_size == 1) return

        if (buffer%protect) then
            if( mini_parallel%verbosity == DEBUGEACH) then
                write(9,*) "buff_pack error: my buffer was write protected (",buffer%vec_num, ")"
            endif
            call MPI_ABORT(MPI_COMM_WORLD,-1,mpinfo) 
        endif
        !Data packing - but only if needed based on the distance between the current
        ! and next vectors in the buffer space

        !the check is done inside this routine so that future openmp block is well defined
        !with no explicit branches
        if (tobuff > mini_parallel%num_wvfunc) then
            if( mini_parallel%verbosity == DEBUGEACH) then
                write(9,*) "buff_pack: exceeded num_wvfunc, I refuse to pack this buffer (",tobuff,")"
            endif
            return
        endif

        do inode = 0,mini_parallel%countcomm-1
        jst   = mini_parallel%jsendp(mini_parallel%sources(inode))-1
        jelem = mini_parallel%jsendp(mini_parallel%sources(inode)+1) - jst - 1
        ! Manual packing
        do j_comm = 1, jelem
        buffer%send(j_comm+mini_parallel%sdispls(inode)) = &
            p(mini_parallel%senrows(j_comm+jst),tobuff)
        enddo
        enddo
        !
        ! Add information about the tobuff that you just buffered
        ! currently only within blksize. 
        ! TODO/maybe : higher level awareness?
        if (current > tobuff) then !ugly hack for labeling
            buffer%vec_num=current
        else
            buffer%vec_num=tobuff
        endif

    end subroutine buff_pack !}}}

    subroutine buff_unpack(buffer,mini_parallel,workvec) !{{{
        use mini_constants
        use mini_parallel_data_mod
        implicit none
        !
        !  Input/Output variables:
        !
        !
        class (mini_buffer) :: buffer
        type (mini_parallel_data), intent(in) :: mini_parallel
        ! sorry, still huge vector since comm neigh can
        real(dp), intent(inout) :: workvec(mini_parallel%workvec_edge+1)
        !
        ! Work variables
        !
        integer :: inode,jst,jelem,j_comm
        !
        !basic protection against 1 PE jobs
        if (mini_parallel%group_size == 1) return
        do inode = 0,mini_parallel%countcomm-1
        jst   = mini_parallel%irecvp(mini_parallel%destinations(inode))-1
        jelem = mini_parallel%irecvp(mini_parallel%destinations(inode)+1) - jst - 1
        ! Manual unpacking
        ! AJB: Must make sure that this is mapped to the correct memory package (future)
!!$OMP DO SCHEDULE(RUNTIME)
        do j_comm = 1, jelem
        workvec(mini_parallel%irecvp(mini_parallel%destinations(inode))+j_comm-1)=&
            buffer%receive(mini_parallel%rdispls(inode) +j_comm)
        enddo
!!$OMP END DO
        enddo

        ! done unpacking
    end subroutine buff_unpack !}}}

    !========================Communication part======================!
    subroutine buff_comm_non_block(buffer,mini_parallel) !{{{
        use mini_constants
        use mini_parallel_data_mod
        use mpi
        implicit none

        class (mini_buffer) :: buffer
        type (mini_parallel_data), intent(in) :: mini_parallel

        integer :: mpinfo
        !this ifdef could be switched to runtime check
#ifdef NOCOLLECTIVE
        integer :: msgtype
        integer :: icomm,inode
        integer :: req(0:2*mini_parallel%group_size)
#else
        integer :: nonblock_req
#endif

        !basic protection against 1 PE jobs
        if (mini_parallel%group_size == 1) return

#ifndef NOCOLLECTIVE
        ! ineighbor alltoallv
        call MPI_Ineighbor_alltoallv(buffer%send,mini_parallel%sendcounts,mini_parallel%sdispls,MPI_DOUBLE_PRECISION,&
            buffer%receive,mini_parallel%recvcounts,mini_parallel%rdispls,MPI_DOUBLE_PRECISION,&
            mini_parallel%group_comm_topo,nonblock_req,mpinfo)
        ! copy req to buffer struct
        buffer%request=nonblock_req
        buffer%request_number=1
        call buffer%wprotect()
#else
        ! isend/irecv

        icomm = 0
        req(:) = 0
        ! manually post receives
        do inode = 0,mini_parallel%countcomm-1
        msgtype = 777 + mini_parallel%destinations(inode)
        call MPI_IRECV(buffer%receive(1+mini_parallel%rdispls(inode)), mini_parallel%recvcounts(inode),MPI_DOUBLE_PRECISION, &
            mini_parallel%destinations(inode), msgtype, mini_parallel%group_comm_topo,  req(icomm), mpinfo )
        icomm = icomm + 1
        enddo
        ! manually post sends 
        do inode = 0,mini_parallel%countcomm-1
        msgtype = 777 + mini_parallel%group_iam
        call MPI_ISEND(buffer%send(1+mini_parallel%sdispls(inode)), mini_parallel%sendcounts(inode),MPI_DOUBLE_PRECISION, &
            mini_parallel%sources(inode), msgtype, mini_parallel%group_comm_topo,  req(icomm), mpinfo )
        icomm = icomm + 1
        enddo
        ! copy req to buffer struct
        buffer%request=req
        buffer%request_number=icomm
        call buffer%wprotect()
#endif
    end subroutine buff_comm_non_block !}}}

    subroutine buff_comm_block(buffer,mini_parallel) !{{{
        use mini_constants
        use mini_parallel_data_mod
        use mpi
        implicit none

        class (mini_buffer) :: buffer
        type (mini_parallel_data), intent(in) :: mini_parallel

        integer :: mpinfo
        !this ifdef could be switched to runtime check
#ifdef NOCOLLECTIVE
        integer :: msgtype
        integer :: icomm,inode
#endif

        !basic protection against 1 PE jobs
        if (mini_parallel%group_size == 1) return

#ifdef NOCOLLECTIVE
        ! send/recv?? I'm not implementing this.
        write(*,*) 'CRITICAL!!! Blocking communication without neighbor_alltoallv is not implemented!'
        call MPI_ABORT(MPI_COMM_WORLD,-1,mpinfo) 
#else
        ! neighbor alltoallv
        call MPI_neighbor_alltoallv(buffer%send,mini_parallel%sendcounts,mini_parallel%sdispls,MPI_DOUBLE_PRECISION,&
            buffer%receive,mini_parallel%recvcounts,mini_parallel%rdispls,MPI_DOUBLE_PRECISION,&
            mini_parallel%group_comm_topo,mpinfo)
        !        call buffer%wprotect()
#endif
    end subroutine buff_comm_block !}}}

    subroutine buff_test_non_block(buffer)
        use mini_constants
        use mpi
        implicit none
        class (mini_buffer) :: buffer

        !  logical, intent(out) :: answer
        logical :: answer
        integer :: comm_status(MPI_STATUS_SIZE,buffer%request_number)
        integer :: mpinfo

        if (buffer%protect) then
            ! is it bad practice to call testall of 1 request?

            call MPI_TESTALL(buffer%request_number, buffer%request, answer, comm_status, mpinfo)
            if (answer) then
                call buffer%wrelease()
            endif
        endif

    end subroutine buff_test_non_block !}}}

    subroutine buff_wait_non_block(buffer) !,vec_num) !{{{
        use mini_constants
        use mpi
        implicit none
        ! does mpi test if check_mode is test
        ! if done removes write protect from buffer
        class (mini_buffer) :: buffer
        ! integer, intent(in) :: vec_num !compare this number with the one in the buffer
        !
        ! work variables:
        !
        integer :: comm_status(MPI_STATUS_SIZE,buffer%request_number)
        integer :: mpinfo

        ! if(buffer%vec_num /= vec_num ) then
        !   if( mini_parallel%verbosity == DEBUGEACH) then
        !     write(9,*) "buff_wait_non_block: I got a vector that I did not expect"
        !     write(9,*) "buff_wait_non_block: got -",buffer%vec_num,"while expecting",vec_num
        !    endif
        ! endif

        ! is it bad practice to call waitall of 1 request?
        call MPI_WAITALL(buffer%request_number, buffer%request, comm_status, mpinfo)
        call buffer%wrelease() 

    end subroutine buff_wait_non_block !}}}

    subroutine buff_nothing(buffer)
        !does nothing, just releases the write protect if it is on
        use mini_constants
        use mpi
        implicit none
        class (mini_buffer) :: buffer
        if (buffer%protect) then
            call buffer%wrelease()
        endif
    end subroutine

end module mini_buffers_mod
