! proper header
!> @description

module mini_global_data_mod
    use mini_constants
    use mini_parallel_data_mod
    use mini_grid_mod
    use mini_buffers_mod
    
    implicit none
    
    !parallel data
    type (mini_parallel_data) :: mini_parallel

    !the grid
    type (mini_grid_data) :: mini_grid

    !buffer - hard coded !!!
    type (mini_buffer) :: mini_buffers(NUM_BUFFER)

    !work vector (to be mapped to HBM!)
    real(dp), dimension(:), pointer :: workvec

    !mock potential vector (to be mapped to HBM?)
    real(dp), dimension(:), pointer :: mockpot

    !wavefunction vectors
    real(dp), dimension(:,:), pointer :: wvfuncs

    contains

    !> @description initialize the derived types used as global data
    !! these are:
    !! grid
    !! parallel
    !! buffers
    !! also these arrays:
    !! workvec
    !! mockpot
    !! wavefunctions
    !! also init MPI (collective call) and trace analyzer if enabled
    subroutine mini_global_data_init(ierr) !{{{
            implicit none

            integer, intent (out) :: ierr

            ierr = 0
            ! init mpi environment
            !
            call mini_parallel%init()
            call mini_parallel%comm_init(ierr)
            !ITAC API, does nothing if ITAC not enabled
            call mini_parallel%itac_init(ierr)

            nullify(workvec)
            nullify(mockpot)
            nullify(wvfuncs)

            call mini_grid%init()

            call mini_global_set_defaults()

        end subroutine mini_global_data_init !}}}


        subroutine mini_global_data_destroy() !{{{
            implicit none

            call mini_parallel%destroy()
            call mini_grid%destroy()
            call mini_destroy_all_buffers()

            if (associated(workvec)) then
                deallocate(workvec)
            endif

            if (associated(mockpot)) then
                deallocate(mockpot)
            endif

            if (associated(wvfuncs)) then
                deallocate(wvfuncs)
            endif

        end subroutine mini_global_data_destroy!}}}

        subroutine mini_init_all_buffers()
            implicit none
            integer :: num_buffers,ii

            num_buffers = mini_parallel%blksize 

            do ii=1,num_buffers !init each of the buffers
               call mini_buffers(ii)%init()
            end do
        end subroutine mini_init_all_buffers

        subroutine mini_destroy_all_buffers()
            implicit none
            integer :: num_buffers,ii

            num_buffers = mini_parallel%blksize

            do ii=1,num_buffers
               call mini_buffers(ii)%destroy()
            end do
        end subroutine mini_destroy_all_buffers

        subroutine mini_global_set_defaults() !{{{
            implicit none

            mini_grid%norder = DEF_ORDER
            mini_grid%stepin = DEF_STEP
            mini_grid%rmax = DEF_RMAX

            mini_grid%d_shape_param = mini_grid%rmax
            mini_grid%i_shape_param = 0
            mini_grid%domain_shape = 0

            mini_parallel%blksize = size(mini_buffers)
            mini_parallel%num_wvfunc = NUM_WAVEFUN
            mini_parallel%num_outer = NUM_OUTER
            mini_parallel%verbosity = MINIMAL

            mini_grid%shift = zero
            !mini_grid%shift = half

        end subroutine mini_global_set_defaults !}}}

        !> @description This is the main exit subroutine of miniparsec
        !! it also serves as an error handler. If ierr>0 it calls MPI_FINALIZE
        subroutine mini_exit_err(ierr) !{{{

            use mpi
            implicit none
            integer, intent (inout) :: ierr
            integer :: mpinfo
            ! whether or not a file is open
            logical :: lopen

! #ifdef DEBUG
!             write(*,*) 'calling exit err with ierr=',ierr
! #endif
            
            call MPI_allREDUCE(MPI_IN_PLACE,ierr,1,MPI_INTEGER,MPI_MAX,mini_parallel%comm,mpinfo) !AJB: mpi max should work

            if(ierr /= 0) then !we are quitting

                !if any text files are open let's close them first since it should induce a flush !{{{
                inquire(unit=9,opened=lopen) !{{{
                if(lopen) then
                    if (ierr > 0) then
                        write(9,*) ' Abort! Error code ierr = ',ierr
                    elseif (ierr == -666) then !special case:program end
                        write(9,*) ' Program End, finalizing '
                    endif
                    close(unit=9)
                else
                    if (ierr > 0) then
                        write(*,*) ' PROC ID:',mini_parallel%iam ,'Abort! Error code ierr = ',ierr
                    endif
                endif !}}}

                if (mini_parallel%iammaster) then !{{{
                    inquire(unit=7,opened=lopen)
                    if (lopen) then
                        write(7,*) ' Program End, finalizing '
                        close(unit=7)
                    endif
                endif !}}} 
                !}}}

                call mini_global_data_destroy()

                call MPI_FINALIZE(mpinfo)
                stop
            endif

        end subroutine mini_exit_err!}}}

end module mini_global_data_mod

