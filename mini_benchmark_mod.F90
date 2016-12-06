! proper header
! description

module mini_benchmark_mod

contains

    subroutine mini_run_benchmark(n_outer,walltime,ierr) !{{{
        use mini_constants
        use mini_common_mod
#ifdef OMP
        use omp_lib
#endif
        implicit none
#ifdef ITAC
        include 'VT.inc'
        include 'mini_vt.inc'
#endif

        integer, intent(in) :: n_outer
        real(dp), intent(out) :: walltime(2)
        integer, intent(inout) :: ierr

        integer :: ii,mvierr
        

        real(dp) :: timestamp0, timestamp1

        ierr = 0
        ii = 0 
        walltime = zero

#ifdef ITAC
!if (mod(ii,2)) then
! we would like this to be inside parallel region but I think it does not work
call VTTIMESYNC(mvierr)
!endif
#endif

#ifdef OMP
!$OMP PARALLEL DEFAULT(SHARED) FIRSTPRIVATE(ii,mvierr)
#endif

#ifdef ITAC
    call VTBEGIN( vt_outside, mvierr)
#endif

     do ii=1,n_outer
!$OMP MASTER
        call my_wtime(timestamp0)
!$OMP END MASTER
        
            call outside_kernel()

!$OMP MASTER
        call my_wtime(timestamp1)
        walltime(1)=walltime(1)+(timestamp1-timestamp0)
!$OMP END MASTER

#ifdef ITAC
    call VTEND( vt_outside, mvierr)
#endif
! removing this barrier
! !$OMP BARRIER

        call kernel()


#ifdef ITAC
    call VTBEGIN( vt_outside, mvierr)
#endif

!$OMP MASTER
        call my_wtime(timestamp0)
        walltime(2)=walltime(2)+(timestamp0-timestamp1)
!$OMP END MASTER

     end do
#ifdef ITAC
    call VTEND( vt_outside, mvierr)
#endif
#ifdef OMP
!$OMP END PARALLEL
#endif

     end subroutine mini_run_benchmark !}}}

    subroutine outside_kernel() !{{{
        use mini_constants
        use mini_global_data_mod
        implicit none
        integer :: kk,jj,ii,mvierr,dummy
        real(dp) :: tempmax(mini_parallel%num_wvfunc)

        tempmax=zero

        !$OMP DO COLLAPSE (2) SCHEDULE(STATIC)
        do kk= 1,OUTSIDE_WORK
            do ii = 1,mini_parallel%num_wvfunc
                    tempmax(ii) = max(one/maxval(wvfuncs(1:mini_parallel%mydim,ii)),tempmax(ii))
            end do
        end do
        !$OMP END DO 

            !doing some allreduce to simulate non kernel communication - and force MPI sync
        !$OMP MASTER
        call psum(tempmax,mini_parallel%num_wvfunc,mini_parallel%group_size,mini_parallel%group_comm_topo)
        !$OMP END MASTER
        !$OMP BARRIER

        !$OMP DO SIMD COLLAPSE (2) SCHEDULE(RUNTIME)
            do ii = 1,mini_parallel%num_wvfunc
              do jj = 1,mini_parallel%mydim
                    wvfuncs(jj,ii) = wvfuncs(jj,ii)*tempmax(ii)
              end do
            end do
        !$OMP END DO SIMD NOWAIT


    end subroutine outside_kernel !}}}

    subroutine kernel() !{{{
        use mini_constants
        use mini_global_data_mod
        use mini_buffers_mod
#ifdef OMP
        use omp_lib
#endif
        implicit none
#ifdef ITAC
        include 'VT.inc'
        include 'mini_vt.inc'
#endif

        integer :: ii,mvierr
        integer :: num_prebuffer,prebuff_index,k
        integer :: current_buffer,next_buffer,num_buffers
        real(dp) :: tempvec(mini_parallel%mydim)

     if( mini_parallel%group_size > 1 ) then !prebuffer if more than one PE !{{{
#ifdef ITAC
    call VTBEGIN( vt_prebuffer, mvierr)
#endif

!$OMP MASTER
#ifdef BLOCKING
        num_buffers=1
        num_prebuffer=0
#else
        num_buffers=mini_parallel%blksize
        num_prebuffer=min(num_buffers-1,mini_parallel%num_wvfunc)
        num_prebuffer=max(1,num_prebuffer)

        do k=1,num_prebuffer
            prebuff_index=mod(k,num_buffers)+1
            call mini_buffers(prebuff_index)%packing(mini_parallel,wvfuncs,k,k)
! write(9,*) 'packed prebuffer - ', mod(k,num_buffers)+1
! write(9,*) 'with wavefunction  - ', k 
            call mini_buffers(prebuff_index)%buffcomm(mini_parallel)
        end do
#endif
#ifdef DEBUG
write(9,*) 'packed buffers up to - ', mod(k,num_buffers)+1
write(9,*) 'total prebuffering - ', num_prebuffer
#endif
!$OMP END MASTER

#ifdef ITAC
    call VTEND( vt_prebuffer, mvierr)
#endif
      endif !prebuffering !}}}

    do ii = 1,mini_parallel%num_wvfunc  
    !{{{ Buffering
    ! buffer the next vector if not prebufferd
     if( mini_parallel%group_size > 1 ) then 

!$OMP MASTER
            current_buffer=mod(ii,num_buffers)+1
#ifdef DEBUG
write(9,*) 'current_buffer',current_buffer
write(9,*) 'current wavefunction ',ii
#endif
!$OMP END MASTER

        if(ii+num_buffers-1>mini_parallel%num_wvfunc) then

!$OMP MASTER
            !no need to buffer anymore, test all buffers?
!            call mini_buffers(current_buffer)%testing()
            next_buffer = -1
!$OMP END MASTER

        else

#ifdef ITAC
    call VTBEGIN( vt_buffer, mvierr)
#endif

!$OMP MASTER
            next_buffer=mod(ii+num_buffers+num_prebuffer,num_buffers)+1

#ifdef DEBUG
write(9,*) 'packing next buffer ',next_buffer
write(9,*) 'with wave function ',ii+num_prebuffer
#endif
            call mini_buffers(next_buffer)%packing(mini_parallel,wvfuncs,ii,ii+num_buffers-1)
            !what if I test for my buffer before sending more stuff?
!            call mini_buffers(current_buffer)%testing()
#ifndef BLOCKING
            call mini_buffers(next_buffer)%buffcomm(mini_parallel)
#endif
!$OMP END MASTER
        endif
#ifdef ITAC
    call VTEND( vt_buffer, mvierr)
#endif
    endif
 !}}}
 !{{{ Diag and waiting for MPI
#ifdef ITAC
    call VTBEGIN( vt_diag, mvierr)
#endif
    
    ! Do some work that does not need communication
    ! static schedule is not recommended here!
    ! and start writing to tempvec
    ! also copies data to shared-memory workvec
    ! assuming here all workvec is on the same NUMA node so first touch is uneeded
    ! will have conseqences when testing on MCDRAM
    call mock_diagonal(ii,tempvec)

#ifdef ITAC
    call VTEND( vt_diag, mvierr)
#endif

    
!$OMP MASTER
    if (mini_parallel%group_size > 1) then
        ! let's just do old school manual progression here. commenting out:
        ! call mini_buffers(current_buffer)%testing()
        ! if (next_buffer>0) then
        !     call mini_buffers(next_buffer)%testing()
        ! endif
#ifdef BLOCKING
        call mini_buffers(next_buffer)%buffcomm(mini_parallel)
#else
        call mini_buffers(current_buffer)%waiting()
#endif
        ! since eventually there will not be an unpacking event the context is
        ! left as 'idle' in itac
        call mini_buffers(current_buffer)%unpacking(mini_parallel,workvec)
    endif 
!$OMP END MASTER

!}}}
!{{{ Inner (no need for MPI) grid points
#ifdef ITAC
    call VTBEGIN( vt_mocklap_inner, mvierr)
#endif

  call mock_laplacian(tempvec,1,mini_parallel%inter1)

#ifdef ITAC
    call VTEND( vt_mocklap_inner, mvierr)
#endif
! I dont know if this really helps or not: right now I think it does not
! !$OMP MASTER
!    if (mini_parallel%group_size > 1) then
!        !        if (mini_parallel%inter1 > mini_parallel%mydim/2 ) then
!         if (next_buffer>0) then
!             call mini_buffers(next_buffer)%testing()
!         endif
!        !        endif
!    endif 
!!$OMP END MASTER

! wait for master to return from mpi:
!$OMP BARRIER !need the barrier here right now because mock_lap has NOWAIT
!}}}


!{{{ Outer grid points
#ifdef ITAC
    call VTBEGIN( vt_mocklap_outer, mvierr)
#endif

    call mock_laplacian(tempvec,mini_parallel%inter1+1,mini_parallel%mydim)
#ifdef ITAC
    call VTEND( vt_mocklap_outer, mvierr)
#endif

!$OMP BARRIER !need the barrier here right now because mock_lap has NOWAIT
!}}}


    end do !ii=1,mini_parallel%num_wvfunc 

    end subroutine kernel !}}}

    subroutine mock_diagonal(wvidx,outvec) !{{{
        use mini_constants
        use mini_global_data_mod
        implicit none

        integer, intent(in) :: wvidx
        real(dp), intent(out) :: outvec(mini_parallel%mydim)

        integer :: ii

        ! Now also copies to workvec
!$OMP DO SIMD SCHEDULE(RUNTIME)    
        do ii = 1, mini_parallel%mydim
            outvec(ii) = mockpot(ii) * wvfuncs(ii,wvidx)
            workvec(ii) = wvfuncs(ii,wvidx)
        enddo
!$OMP END DO SIMD NOWAIT


    end subroutine mock_diagonal !}}}

    subroutine copy_to_workvec(wvidx) !{{{
        use mini_constants
        use mini_global_data_mod
        implicit none

        integer, intent(in) :: wvidx
        integer :: ii

!$OMP SINGLE
            workvec(mini_parallel%workvec_edge+1) = zero
!$OMP END SINGLE NOWAIT

!$OMP DO SIMD SCHEDULE(RUNTIME)    
        do ii = 1, mini_parallel%mydim
        enddo
!$OMP END DO SIMD NOWAIT

    end subroutine copy_to_workvec !}}}

    ! subroutine mock_outer_product()
    ! end subroutine mock_outer_product

    subroutine mock_laplacian(outvec,first_ind,last_ind) !{{{
        use mini_constants
        use mini_global_data_mod
        implicit none

        real(dp), intent(inout) :: outvec(mini_parallel%mydim)
        integer, intent(in) :: first_ind, last_ind

        integer :: ii,irow,shell_index,mm
        integer :: current_shell_shift
        real(dp) :: tmp_shell,shellmemb(LAPDIR)
!$OMP DO &
!$OMP& SCHEDULE(RUNTIME)    & 
!$OMP& PRIVATE(ii,irow,shell_index,shellmemb,tmp_shell,current_shell_shift,mm) 
        do ii = first_ind, last_ind
            irow = mini_parallel%pint(ii)
            tmp_shell = zero
            do shell_index = 1,mini_grid%norder
                do mm = 1,LAPDIR
                current_shell_shift = (shell_index-1)*LAPDIR
                shellmemb(mm) = workvec(mini_parallel%neibs(mm+current_shell_shift,irow))
                shellmemb(mm) = shellmemb(mm) * mini_grid%coeff2_1D(mm+current_shell_shift)
                end do
                tmp_shell = tmp_shell+sum(shellmemb(:))
            end do
            outvec(irow)=outvec(irow)+tmp_shell
        end do
!$OMP END DO NOWAIT
        
    end subroutine mock_laplacian !}}}

    subroutine mini_run_nested_benchmark(n_outer,walltime,ierr) !{{{
        use mini_constants
        use mini_common_mod
#ifdef OMP
        use omp_lib
#endif
        implicit none
#ifdef ITAC
        include 'VT.inc'
        include 'mini_vt.inc'
#endif

        integer, intent(in) :: n_outer
        real(dp), intent(out) :: walltime(2)
        integer, intent(inout) :: ierr

        integer :: ii,mvierr
        

        real(dp) :: timestamp0, timestamp1

        ierr = 0
        ii = 0 
        walltime = zero

#ifdef ITAC
!if (mod(ii,2)) then
! we would like this to be inside parallel region but I think it does not work
call VTTIMESYNC(mvierr)
!endif
#endif

#ifdef OMP
!$OMP PARALLEL DEFAULT(SHARED) FIRSTPRIVATE(ii,mvierr)
#endif

#ifdef ITAC
    call VTBEGIN( vt_outside, mvierr)
#endif

     do ii=1,n_outer
!$OMP MASTER
        call my_wtime(timestamp0)
!$OMP END MASTER
        
            call outside_kernel()

!$OMP MASTER
        call my_wtime(timestamp1)
        walltime(1)=walltime(1)+(timestamp1-timestamp0)
!$OMP END MASTER

#ifdef ITAC
    call VTEND( vt_outside, mvierr)
#endif
! removing this barrier
! !$OMP BARRIER

        call kernel()


#ifdef ITAC
    call VTBEGIN( vt_outside, mvierr)
#endif

!$OMP MASTER
        call my_wtime(timestamp0)
        walltime(2)=walltime(2)+(timestamp0-timestamp1)
!$OMP END MASTER

     end do
#ifdef ITAC
    call VTEND( vt_outside, mvierr)
#endif
#ifdef OMP
!$OMP END PARALLEL
#endif

     end subroutine mini_run_nested_benchmark !}}}

    subroutine outside_nested_kernel() !{{{
        use mini_constants
        use mini_global_data_mod
        implicit none
        integer :: kk,jj,ii,mvierr,dummy
        real(dp) :: tempmax(mini_parallel%num_wvfunc)

        tempmax=zero

        !$OMP DO COLLAPSE (2) SCHEDULE(STATIC)
        do kk= 1,OUTSIDE_WORK
            do ii = 1,mini_parallel%num_wvfunc
                    tempmax(ii) = max(one/maxval(wvfuncs(1:mini_parallel%mydim,ii)),tempmax(ii))
            end do
        end do
        !$OMP END DO 

            !doing some allreduce to simulate non kernel communication - and force MPI sync
        !$OMP MASTER
        call psum(tempmax,mini_parallel%num_wvfunc,mini_parallel%group_size,mini_parallel%group_comm_topo)
        !$OMP END MASTER
        !$OMP BARRIER

        !$OMP DO SIMD COLLAPSE (2) SCHEDULE(RUNTIME)
            do ii = 1,mini_parallel%num_wvfunc
              do jj = 1,mini_parallel%mydim
                    wvfuncs(jj,ii) = wvfuncs(jj,ii)*tempmax(ii)
              end do
            end do
        !$OMP END DO SIMD NOWAIT


    end subroutine outside_nested_kernel !}}}

    subroutine nested_kernel() !{{{
        use mini_constants
        use mini_global_data_mod
        use mini_buffers_mod
#ifdef OMP
        use omp_lib
#endif
        implicit none
#ifdef ITAC
        include 'VT.inc'
        include 'mini_vt.inc'
#endif

        integer :: ii,mvierr
        integer :: num_prebuffer,prebuff_index,k
        integer :: current_buffer,next_buffer,num_buffers
        real(dp) :: tempvec(mini_parallel%mydim)

     if( mini_parallel%group_size > 1 ) then !prebuffer if more than one PE !{{{
#ifdef ITAC
    call VTBEGIN( vt_prebuffer, mvierr)
#endif

!$OMP MASTER
#ifdef BLOCKING
        num_buffers=1
        num_prebuffer=0
#else
        num_buffers=mini_parallel%blksize
        num_prebuffer=min(num_buffers-1,mini_parallel%num_wvfunc)
        num_prebuffer=max(1,num_prebuffer)

        do k=1,num_prebuffer
            prebuff_index=mod(k,num_buffers)+1
            call mini_buffers(prebuff_index)%packing(mini_parallel,wvfuncs,k,k)
! write(9,*) 'packed prebuffer - ', mod(k,num_buffers)+1
! write(9,*) 'with wavefunction  - ', k 
            call mini_buffers(prebuff_index)%buffcomm(mini_parallel)
        end do
#endif
#ifdef DEBUG
write(9,*) 'packed buffers up to - ', mod(k,num_buffers)+1
write(9,*) 'total prebuffering - ', num_prebuffer
#endif
!$OMP END MASTER

#ifdef ITAC
    call VTEND( vt_prebuffer, mvierr)
#endif
      endif !prebuffering !}}}

    do ii = 1,mini_parallel%num_wvfunc  
    !{{{ Buffering
    ! buffer the next vector if not prebufferd
     if( mini_parallel%group_size > 1 ) then 

!$OMP MASTER
            current_buffer=mod(ii,num_buffers)+1
#ifdef DEBUG
write(9,*) 'current_buffer',current_buffer
write(9,*) 'current wavefunction ',ii
#endif
!$OMP END MASTER

        if(ii+num_buffers-1>mini_parallel%num_wvfunc) then

!$OMP MASTER
            !no need to buffer anymore, test all buffers?
!            call mini_buffers(current_buffer)%testing()
            next_buffer = -1
!$OMP END MASTER

        else

#ifdef ITAC
    call VTBEGIN( vt_buffer, mvierr)
#endif

!$OMP MASTER
            next_buffer=mod(ii+num_buffers+num_prebuffer,num_buffers)+1

#ifdef DEBUG
write(9,*) 'packing next buffer ',next_buffer
write(9,*) 'with wave function ',ii+num_prebuffer
#endif
            call mini_buffers(next_buffer)%packing(mini_parallel,wvfuncs,ii,ii+num_buffers-1)
            !what if I test for my buffer before sending more stuff?
!            call mini_buffers(current_buffer)%testing()
#ifndef BLOCKING
            call mini_buffers(next_buffer)%buffcomm(mini_parallel)
#endif
!$OMP END MASTER
        endif
#ifdef ITAC
    call VTEND( vt_buffer, mvierr)
#endif
    endif
 !}}}
 !{{{ Diag and waiting for MPI
#ifdef ITAC
    call VTBEGIN( vt_diag, mvierr)
#endif
    
    ! Do some work that does not need communication
    ! static schedule is not recommended here!
    ! and start writing to tempvec
    ! also copies data to shared-memory workvec
    ! assuming here all workvec is on the same NUMA node so first touch is uneeded
    ! will have conseqences when testing on MCDRAM
    call mock_diagonal(ii,tempvec)

#ifdef ITAC
    call VTEND( vt_diag, mvierr)
#endif

    
!$OMP MASTER
    if (mini_parallel%group_size > 1) then
        ! let's just do old school manual progression here. commenting out:
        ! call mini_buffers(current_buffer)%testing()
        ! if (next_buffer>0) then
        !     call mini_buffers(next_buffer)%testing()
        ! endif
#ifdef BLOCKING
        call mini_buffers(next_buffer)%buffcomm(mini_parallel)
#else
        call mini_buffers(current_buffer)%waiting()
#endif
        ! since eventually there will not be an unpacking event the context is
        ! left as 'idle' in itac
        call mini_buffers(current_buffer)%unpacking(mini_parallel,workvec)
    endif 
!$OMP END MASTER

!}}}
!{{{ Inner (no need for MPI) grid points
#ifdef ITAC
    call VTBEGIN( vt_mocklap_inner, mvierr)
#endif

  call mock_laplacian(tempvec,1,mini_parallel%inter1)

#ifdef ITAC
    call VTEND( vt_mocklap_inner, mvierr)
#endif
! I dont know if this really helps or not: right now I think it does not
! !$OMP MASTER
!    if (mini_parallel%group_size > 1) then
!        !        if (mini_parallel%inter1 > mini_parallel%mydim/2 ) then
!         if (next_buffer>0) then
!             call mini_buffers(next_buffer)%testing()
!         endif
!        !        endif
!    endif 
!!$OMP END MASTER

!}}}
!{{{ Outer grid points
! wait for master to return from mpi:
!$OMP BARRIER !need the barrier here right now because mock_lap has NOWAIT

#ifdef ITAC
    call VTBEGIN( vt_mocklap_outer, mvierr)
#endif

    call mock_laplacian(tempvec,mini_parallel%inter1+1,mini_parallel%mydim)


#ifdef ITAC
    call VTEND( vt_mocklap_outer, mvierr)
#endif
!}}}

!$OMP BARRIER !need the barrier here right now because mock_lap has NOWAIT

    end do !ii=1,mini_parallel%num_wvfunc 

    end subroutine nested_kernel !}}}


end module mini_benchmark_mod
