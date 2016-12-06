!============================================================================
!
!                          +-+-+-+-+-+-+-+-+-+-+-+-+
!                            MINI PARSEC - 0.1.x
!                          +-+-+-+-+-+-+-+-+-+-+-+-+
!
! PARESC is a real-space, high-order finite difference electronic structure
! code.  MINI PARSEC is an MPI3 finite difference benchmark that is closely
! related to the main computational kernel of PARSEC, matvec.
!
! Copyright? (C) 2005,2015? Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! St, Fifth Floor, Boston, MA 02110-1301 USA
!-------------------------------------------------------------------------------

program miniparsec
    use mini_constants !> Definition of constants
    use mini_common_mod, ONLY:my_wtime !> For timing
    use mini_global_data_mod  !> Common data in terms of scope, each PE has different data!
    use mini_benchmark_mod !> Main benchmark subroutines
    use mini_io_mod !> Main input/output
    use mpi !> need here because using BCAST             

    implicit none
#ifdef ITAC
    include 'VT.inc' !> trace-analyzer definitions
    include 'mini_vt.inc' !> sadly one needs global variables for the ITAC api
#endif 
    ! other variable definitions
    real(dp) :: walltime_bench(2), tstamp0, tstamp1 !> time stamps for my_wtime
    integer :: mpinfo,ierr !> misc flags
#ifdef ITAC
    integer :: vt_ierr !> ITAC flag
#endif

#ifdef ITAC
    vt_ierr = 0
!    call VTTRACEOFF(vt_ierr) !turn off tracing, not sure it actually works
#endif

    ierr = 0 
    call mini_global_data_init(ierr)
    call mini_exit_err(ierr)
    ! ===============================================================
    ! Start the master clock?
    ! ===============================================================

    call my_wtime(tstamp0)
    if (mini_parallel%iammaster)  then
        call mini_input_read(ierr)
    endif
    ! see if everything went okay for master process
    call mini_exit_err(ierr)
    ! then share the verbosity value and open files if necessary.
    call MPI_BCAST(mini_parallel%verbosity,1,MPI_INTEGER,mini_parallel%masterid,mini_parallel%comm,mpinfo)
    call mini_io_set_verbosity(mini_parallel%verbosity)

    ! set up, includes sharing of data from master process and 
    ! actual set up of run parameters for all PEs
    ! and filling of all important vectors with random numbers
    call mini_setup(ierr)
#ifdef DEBUG
                     ! write(*,*) 'back from mini_setup'
                     ! write(*,*) 'I am', mini_parallel%iam
                     ! write(*,*) 'rmax is now', mini_grid%rmax
                     ! write(*,*) 'stepin is now', mini_grid%stepin
#endif
    call mini_exit_err(ierr)

    ! if needed, generate report on what has been done and what is going to happen
    call my_wtime(tstamp1)
    call mini_finalize_pre(tstamp1-tstamp0)

#ifdef ITAC
    ! Tracing should begin here
!    call VTTRACEON(vt_ierr) 
!    call VTBEGIN(vt_outside,vt_ierr)
#endif

    ! perform the benchmark
    call mini_run_benchmark(mini_parallel%num_outer,walltime_bench,ierr)
    call mini_exit_err(ierr)

#ifdef ITAC
!    call VTEND(vt_outside,vt_ierr)
#endif

    call my_wtime(tstamp0)
    ! generate small report on the benchmark run
    call mini_end_credits(walltime_bench,tstamp0-tstamp1)

    ! Done. mini_exit_err does destruction work
    ierr = -666
    call mini_exit_err(ierr)

    !FAIL SAFE
!    call MPI_ABORT(MPI_COMM_WORLD,-1,mpinfo)

end program miniparsec
