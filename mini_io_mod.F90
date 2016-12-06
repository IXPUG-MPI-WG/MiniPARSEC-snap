! proper header
! description

module mini_io_mod
    use mini_constants

    implicit none

    type cmd_args !{{{

        real(dp) :: r_max
        real(dp) :: grid_spacing
        integer  :: verbosity
        integer  :: norder
        logical  :: ok
        !! TODO: type of communication

    contains

        procedure :: init => cmd_args_init
        procedure :: parse => parse_args

    end type cmd_args !}}}

contains

    subroutine  cmd_args_init(cmd_arguments) !{{{
        implicit none

        class(cmd_args) :: cmd_arguments

        cmd_arguments%r_max = mone
        cmd_arguments%grid_spacing = mone
        cmd_arguments%verbosity = -1
        cmd_arguments%norder = 12
        cmd_arguments%ok = .FALSE.

    end subroutine  cmd_args_init !}}}

    subroutine parse_args(cmd_arguments) !{{{
        implicit none

        class(cmd_args) :: cmd_arguments

        integer :: cmd_count

        character(len=10) :: arg

        cmd_count = command_argument_count()

        select case(cmd_count)

        case (0)
            ! all default run
            cmd_arguments%ok = .TRUE.

        case (2)
            call get_command_argument(1,arg)
            read(arg,*) cmd_arguments%r_max
            call get_command_argument(2,arg)
            read(arg,*) cmd_arguments%grid_spacing

            cmd_arguments%ok = .TRUE.

        case (3)
            call get_command_argument(1,arg)
            read(arg,*) cmd_arguments%r_max
            call get_command_argument(2,arg)
            read(arg,*) cmd_arguments%grid_spacing
            call get_command_argument(3,arg)
            read(arg,*) cmd_arguments%verbosity

            cmd_arguments%ok = .TRUE.

        case (4)
            call get_command_argument(1,arg)
            read(arg,*) cmd_arguments%r_max
            call get_command_argument(2,arg)
            read(arg,*) cmd_arguments%grid_spacing
            call get_command_argument(3,arg)
            read(arg,*) cmd_arguments%verbosity
            call get_command_argument(4,arg)
            read(arg,*) cmd_arguments%norder


            cmd_arguments%ok = .TRUE.

        case default
            ! something is not quite right
            write(*,*) 'Are you sure about this command line input?'
            write(*,*) 'Usage:: miniparsec rmax(real) spacing(real) [verbosity(int) norder(int)]'
            cmd_arguments%ok = .FALSE.
        end select

#ifdef DEBUG
       !always force multifile output in debug compile
       cmd_arguments%verbosity = DEBUGEACH
#endif
    end subroutine parse_args !}}}

    subroutine mini_io_set_verbosity(verbosity) !{{{
        use mini_constants
        use mini_global_data_mod
#ifdef OMP
        ! access openmp functions
        use omp_lib          
#endif
        implicit none

        integer, intent(in) :: verbosity

        integer :: i
        character(len=4) :: idstring

        if (mini_parallel%iammaster) then

        if (verbosity == NOTXT) then
            write(*,*) 'Running miniparsec with no text output, have a nice day!'
        endif

        if (verbosity > NOTXT) then
            open(7,file='miniparsec.out',form='formatted')  
            write(7,*)
            write(7,*) ('=',i=1,65)
            write(7,*)
            write(7,*) ' MINI PARSEC BENCHMARK 0.3.2 '
#ifdef DEBUG
            write(7,*) ' '
            write(7,*) ' !!!!!!  DEBUG BUILD !!!!!!'
            write(7,*) ' '
#else
            write(7,*) ' '
            write(7,*) " Physicists stick to FORTRAN, so that they can share each other's programs, bugs included"
            write(7,*) ' '
#endif

#ifdef NOCOLLECTIVE
            write(7,*) ' Kernel comm. is MPI2 IRECV/ISEND'
#else
            write(7,*) ' Kernel comm. is MPI3 neighborhood collectives'
#endif
            write(7,*) ' '
#ifdef BLOCKING
            write(7,*) ' NOTE: Blocking communication activated. '
            write(7,*) ' '
#endif
            write(7,*) ' MPI: I am registering ', mini_parallel%procs_num,' PEs '
#ifdef OMP
            write(7,*) ' OMP enabled. Max. Threads/PE: ',OMP_GET_MAX_THREADS()
#endif
            write(7,*)
            write(7,*) ('=',i=1,65)
        endif

        endif


        if (verbosity >= DEBUGEACH) then
            if (mini_parallel%iammaster) then
                write(7,*) 'WARNING, High verbosity mode engaged. Timing data unreliable'
            endif

            write(idstring,'(I4.4)') mini_parallel%iam !limited for 10k processes
            open(9,file='proc.'//idstring,form='formatted',status='unknown')
            write(9,*) 'WARNING, High verbosity mode engaged. Timing data unreliable'
            write(9,*) 'Processor No ',mini_parallel%iam,'  has started'
            write(9,*) 'Number of processors: ',mini_parallel%procs_num
#ifdef OMP
            write(9,*) 'Number of threads/processor: ',OMP_GET_MAX_THREADS()

            write(9,*) ('=',i=1,65)
            write(9,*) ' DEBUG: Runtime OMP thread affinity set to ',OMP_GET_PROC_BIND()
            write(9,*) ' where...'
            write(9,*)  "  false state stat ="  , omp_proc_bind_false
            write(9,*)  "  true  state stat ="   , omp_proc_bind_true
            write(9,*) ' alternatively...'
            write(9,*)  "  master   binding =" , omp_proc_bind_master
            write(9,*)  "  close    binding ="  , omp_proc_bind_close
            write(9,*)  "  spread   binding =" , omp_proc_bind_spread
            write(9,*) ('=',i=1,65)
#else
            write(9,*) 'Number of threads/processor: N/A (No openMP) '
#endif
        endif
    end subroutine mini_io_set_verbosity !}}}

    subroutine mini_input_read(ierr)
        use mini_global_data_mod

        integer, intent(out) :: ierr

        type (cmd_args) :: cmd_arguments

        call cmd_arguments%init()

        ierr = 0

            call cmd_arguments%parse()
            if (cmd_arguments%ok) then

                if (cmd_arguments%r_max > zero) then
                    call mini_grid%set_rmax(cmd_arguments%r_max)
! #ifdef DEBUG
!                      write(*,*) 'called set_rmax.'
!                      write(*,*) 'rmax is now', mini_grid%rmax
! #endif
                else
                    call mini_grid%set_rmax(mini_grid%rmax)
                endif

                if (cmd_arguments%grid_spacing > zero) then
                    call mini_grid%set_step(cmd_arguments%grid_spacing)
! #ifdef DEBUG
!                      write(*,*) 'called set_step.'
!                      write(*,*) 'stepin is now', mini_grid%stepin
!                      write(*,*) 'n1 is now', mini_grid%n1
! #endif
                else
! #ifdef DEBUG
!                      write(*,*) 'grid spacing is set to default: ',mini_grid%stepin
! #endif
                    call mini_grid%set_step(mini_grid%stepin)
                endif

                if (cmd_arguments%norder > 0) then
                    !note that we set it to 1/2 of the input value e.g 12->6
                        if ( mod(cmd_arguments%norder,2) /=0) then
                            !            write(*,*) "norder bad"
                        cmd_arguments%norder=cmd_arguments%norder+1
                        endif
                    call mini_grid%set_norder(cmd_arguments%norder)
                endif

                if (cmd_arguments%verbosity >= NOTXT ) then
                    mini_parallel%verbosity = cmd_arguments%verbosity
                    !else - use the hard-coded default
! #ifdef DEBUG
!                     write(*,*) 'master: my verbosity level=',mini_parallel%verbosity
! #endif
                endif
            else
                ! ERROR CODE 70 :: Command line arguments could not be read
#ifdef DEBUG
                    write(*,*) 'bad argumentttssss'
#endif
                ierr = ierr + 70
            endif


        ! rest of the input parameters are BCASTED in mini_setup
    end subroutine mini_input_read

    subroutine mini_finalize_pre(setup_time) !{{{
        use mini_constants
        use mini_global_data_mod

        implicit none

        real(dp), intent(in) :: setup_time

        if (mini_parallel%verbosity > NOTXT) then
         if (mini_parallel%iammaster) then
            write(7,*) ''
            write(7,*) 'Finished setup, commencing benchmark'
            write(7,*) 'Setup time was:', setup_time, 'seconds'
            write(7,*) ''
            write(7,*) '-------------------------'
            write(7,*) 'REQUESTED PARAMETERS ARE:'
            write(7,*) '-------------------------'

            write(7,*) 'Input grid spacing in AU:', mini_grid%stepin
            write(7,*) 'Actual grid spacing in AU:', mini_grid%step(1)
            write(7,*) 'Specified domain radius in AU', mini_grid%rmax
            write(7,*) '(Half-) Order of Finite Difference stencil', mini_grid%norder/2

#ifndef BLOCKING
            write(7,*) 'Number of buffers', mini_parallel%blksize
#else
            write(7,*) ''
#endif
            write(7,*) 'Number of wave functions', mini_parallel%num_wvfunc
            write(7,*) 'Number of outer loop iterations', mini_parallel%num_outer
#ifdef DEBUG
            write(7,*) 'Domain shape ' , mini_grid%domain_shape
            write(7,*) 'Domain shape d(1) param', mini_grid%d_shape_param(1)
            write(7,*) 'Domain shape i param' , mini_grid%i_shape_param
#endif
         endif
        endif

        if(mini_parallel%verbosity == DEBUGEACH) then
            write(9,*) ''
            write(9,*) 'Finished setup, commencing benchmark'
            write(9,*) 'my Setup time was:', setup_time, 'seconds'
            write(9,*) '------------------------'
            write(9,*) 'PARAMETERS IN MY MEMORY:'
            write(9,*) '------------------------'

            write(9,*) 'Input grid spacing in AU:', mini_grid%stepin
            write(9,*) 'Actual grid spacing in AU:', mini_grid%step(1)
            write(9,*) 'Specified domain radius in AU', mini_grid%rmax
            write(9,*) '(Half-) Order of Finite Difference stencil', mini_grid%norder

            write(9,*) 'Number of buffers', mini_parallel%blksize
            write(9,*) 'Number of wave functions', mini_parallel%num_wvfunc
            write(9,*) 'Number of outer loop iterations', mini_parallel%num_outer
#ifdef DEBUG
            write(9,*) 'Domain shape ' , mini_grid%domain_shape
            write(9,*) 'Domain shape d(1) param', mini_grid%d_shape_param(1)
            write(9,*) 'Domain shape i param' , mini_grid%i_shape_param
#endif
            write(9,*) 'MORE USEFUL DATA FOR THIS PE'
        endif

    end subroutine mini_finalize_pre !}}}

    subroutine mini_end_credits(benchmark_walltime,benchmark_time)
        use mini_constants
        use mini_global_data_mod

        implicit none

        real(dp), intent(in) :: benchmark_walltime(2), benchmark_time

         if (mini_parallel%iammaster) then
        if (mini_parallel%verbosity == NOTXT) then
            write(*,*) 'miniparsec run with no output now done, finishing up...'
            write(7,*) 'Time per benchmark loop:', benchmark_time/mini_parallel%num_outer, 'seconds'
        endif

        if (mini_parallel%verbosity > NOTXT) then
            write(7,*) ''
            write(7,*) 'Finished benchmark'
            write(7,*) 'My watch shows that the benchmark took', benchmark_time, 'seconds'
            write(7,*) ''
            write(7,*) 'Time per benchmark loop:', benchmark_time/mini_parallel%num_outer, 'seconds'
            write(7,*) 'Time per benchmark loop per vector:', benchmark_time/mini_parallel%num_outer/mini_parallel%num_wvfunc, 'seconds'
            write(7,*) ''
            write(7,*) 'Time per benchmark loop (method 2):', (benchmark_walltime(1)+benchmark_walltime(2))/mini_parallel%num_outer, 'seconds'
            write(7,*) 'Walltime for just "outside_kernel"/loop', benchmark_walltime(1)/mini_parallel%num_outer, 'seconds'
            write(7,*) 'Walltime for just "kernel"/loop', benchmark_walltime(2)/mini_parallel%num_outer, 'seconds'
            write(7,*) ''
            write(7,*) 'MORE USEFUL DATA'
        endif
         endif

        if(mini_parallel%verbosity == DEBUGEACH) then
            write(9,*) ''
            write(9,*) 'Finished benchmark'
            write(9,*) 'I think the benchmark took', benchmark_time, 'seconds'
            write(9,*) 'USEFUL DATA FOR THIS PE'
            write(9,*) 'USEFUL DATA FOR THIS PE'
        endif

    end subroutine 



end module mini_io_mod
