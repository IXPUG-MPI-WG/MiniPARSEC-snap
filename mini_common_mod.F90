! proper header
! description

module mini_common_mod

contains

    subroutine my_wtime(time)
        use mini_constants
        use mpi
        implicit none

        real(dp), intent(out) :: time
        time = MPI_WTIME()
    end subroutine my_wtime

    subroutine my_random_array(array, ldn, ndim, blksize) !,info)

        use mini_constants
        implicit none
        !
        !  Input/Output variables:
        !
        !  ldn: leading dimension of array
        !  ndim: working dimension of array
        !  blksize: width of array
        integer, intent(in)  :: ldn, ndim, blksize
        !  integer, intent(out) :: info

        real(dp), intent(inout) :: array(ldn, blksize)

        !
        !  Work variables:
        !
        !  idist (input) specifies the distribution of the random numbers
        !      idist = 1:  uniform (0,1)
        !      idist = 2:  uniform (-1,1)
        !      idist = 3:  normal  (0,1)
        integer, save :: idist = 2 

        !  iseed (input/output) INTEGER array, dimension (4)
        !  On input, the seed of the random number generator; 
        !  the array elements must be between 0 and 4095, 
        !  and ISEED(4) must be odd. On output, the seed is updated
        !              .
        integer, save :: iseed(4)
        data iseed/1,1000,2000,4095/
        !  data iseed/1,3000,700,4095/

        integer :: jcol
        !
        !  External subroutines:
        external  dlarnv

        !-------------------------------------------------------------------
        !!! the following are not necessary at all, since idist is set to 2 
        !
        ! if (idist /= 1 .and. idist /= 2  .and. idist /= 3) then
        !    write (9,*) "ERROR: idist input to random_array is illegal"
        !    write (9,*) "idist must be 1 or 2 or 3"
        !    ierr = idist
        !    return
        ! elseif (mod(iseed(4), 2) == 0) then
        !    write (9,*) "ERROR: iseed(4) must be odd for random_array()"
        !    ierr = -1
        !    return
        ! endif
        ! AJB: uncomment this to find a bug later.
        !info=666

        do jcol = 1, blksize
        call dlarnv(idist, iseed, ndim, array(1,jcol))
        enddo

    end subroutine my_random_array
    !
    !=====================================================================


end module mini_common_mod
