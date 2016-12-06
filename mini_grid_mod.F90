! proper header
! description

module mini_grid_mod
    use mini_constants

    type mini_grid_data

        ! boundary sphere size (atomic units) - for confined systems
        real(dp) :: rmax
        ! user provided  grid spacing (atomic units) - h 
        real(dp) :: stepin
        ! step = program corrected h in each direction
        real(dp) :: step(3)
        ! hcub = h**3, hcub2 = 2/(h**3),h_2 = h**2
        real(dp) :: hcub, hcub2, h_2
        ! total number of effective grid points
        integer :: ndim
        ! total number of grid points in irreducible wedge
        integer :: nwedge
        ! number of neighbors used on one side in numerical derivative
        integer :: norder

        ! nxmax-nxmin+1 -- number of grid points along x-axis (including
        ! points for derivative); same for y and z directions
        integer :: nxmin, nxmax
        integer :: nymin, nymax
        integer :: nzmin, nzmax
        integer :: nxyz !
        ! n1,n2,n3 = number of grid points along x-axis, y-axis, and
        ! z-axis respectively (the ones that are actually used)
        integer :: n1, n2, n3

        ! 3d to 1d conversion table for grid point indexing - gives the
        ! 1d index based on the 3-d position (in h units, relative to the
        ! origin). Returns value of ndim+1 if point is outside the boundary sphere
        ! indexg(i,j,k) = ndim+1 for points outside the current boundary sphere
        ! indexg(i,j,k) = 1 to ndim for points inside the current boundary sphere
        integer, dimension (:,:,:), pointer :: indexg
        ! indexw(i,j,k) : has the same meaning as indexg, but it is
        ! defined only in the irreducible wedge, after considering
        ! symmetry operations in the reduced, Abelian subgroup
        integer, dimension (:,:,:), pointer :: indexw
        ! 1d to 3d conversion tables for grid point indexing - giving INDEX
        ! These arrays retrieve the three 3d indices for the ith 1d grid point
        ! 1d to 3d conversion tables for grid point indexing - giving
        ! cartesian coordinates (for irreducible wedge, replace fx with kx etc):
        ! xx = (shift(1) + kx)*h, yy=(shift(2) + ky)*h, zz=(shift(3) + kz)*h
        integer, dimension (:), pointer :: kx, ky, kz
        integer, dimension (:), pointer :: fx, fy, fz
        ! shift of grid points
        real(dp), dimension(3) :: shift
        ! irreducible wedge parameters
        integer, dimension (:), pointer :: rindex

        ! number of additional directions for calculating derivatives:
        integer lap_dir_num

        ! inverse of the normalized lattice vectors matrix.
        real(dp), dimension(3,3) :: grad_bvec_norm

        ! using mock laplacian, coefficents are random
        ! without central term
        real(dp), dimension (:), pointer :: coeff2_1D  !(0:LAPDIR*norder)


        ! work array for grid partitioning
        integer, dimension(:), pointer :: ist

        ! flag for domain shape in cluster BCs
        ! possible shapes:
        ! 0 - sphere, this is the default
        ! 1 - ellipsoid
        ! 2 - cylindrical
        ! 3 - box
        ! more options can be added.

        integer :: domain_shape

        ! Parameters relevant to the cluster shape
        ! The size of the array and the meaning of its parameters
        ! depend on the cluster shape.
        !
        ! Sphere: d_shape_param(1) = radius
        ! Ellipse: d_shape_param(1) = x radius
        !          d_shape_param(2) = y radius
        !          d_shape_param(3) = z radius
        ! Cylinder: d_shape_param(1) = radius
        !           d_shape_param(2) = length (centered at 0)
        !           i_shape_param(1) = orientation
        !                              (x = 1, y = 2, z = 3)
        ! Box: d_shape_param(1) = x length (centered at 0)
        !      d_shape_param(2) = y length (centered at 0)
        !      d_shape_param(3) = z length (centered at 0)

        real(dp), dimension(3)  :: d_shape_param
        integer, dimension(1) :: i_shape_param
        ! IMPORTANT: if the dimesions of these shape_param arrays
        ! are changed, update the broadcast statements in init_var
        !

    contains

        procedure :: init => mini_grid_init 
        procedure :: destroy => mini_grid_destroy 
        procedure :: set_rmax => set_rmax 
        procedure :: set_norder => set_norder
        procedure :: set_step => set_step 
        procedure :: set_index => set_index
        procedure :: set_wedge => set_wedge 
        procedure :: set_ist => set_ist 
        procedure :: set_ndim => set_ndim 
        procedure :: setup_mock_laplacian => mock_lap_setup


    end type mini_grid_data

contains

    subroutine mini_grid_init(mini_grid)
        implicit none
        class (mini_grid_data) :: mini_grid

        integer, parameter :: int_def = -1

        nullify(mini_grid%indexg) 
        nullify(mini_grid%indexw) 
        nullify(mini_grid%kx) 
        nullify(mini_grid%ky) 
        nullify(mini_grid%kz) 
        nullify(mini_grid%fx) 
        nullify(mini_grid%fy) 
        nullify(mini_grid%fz) 
        nullify(mini_grid%rindex) 
        nullify(mini_grid%ist) 

        mini_grid%step(:)=zero

        mini_grid%nxmin = int_def
        mini_grid%nxmax = int_def
        mini_grid%nymin = int_def
        mini_grid%nymax = int_def
        mini_grid%nzmin = int_def
        mini_grid%nzmax = int_def

        mini_grid%lap_dir_num = 0
    end subroutine mini_grid_init

    subroutine mini_grid_destroy (mini_grid)
        implicit none
        class (mini_grid_data) :: mini_grid

        if (associated (mini_grid%indexg)) deallocate (mini_grid%indexg) 
        if (associated (mini_grid%indexw)) deallocate (mini_grid%indexw)

        if (associated (mini_grid%kx)) deallocate (mini_grid%kx) 
        if (associated (mini_grid%ky)) deallocate (mini_grid%ky) 
        if (associated (mini_grid%kz)) deallocate (mini_grid%kz)

        if (associated (mini_grid%fx)) deallocate (mini_grid%fx)
        if (associated (mini_grid%fy)) deallocate (mini_grid%fy)
        if (associated (mini_grid%fz)) deallocate (mini_grid%fz)

        if (associated (mini_grid%rindex)) deallocate (mini_grid%rindex)


        if (associated(mini_grid%ist)) deallocate(mini_grid%ist)

    end subroutine mini_grid_destroy

    subroutine set_rmax(mini_grid,rmax)
        implicit none

        class (mini_grid_data):: mini_grid
        real(dp), intent (in) :: rmax

        mini_grid%rmax = rmax
        mini_grid%d_shape_param = rmax

    end subroutine set_rmax

    subroutine set_norder(mini_grid,norder)
        implicit none

        class (mini_grid_data):: mini_grid
        integer, intent (in) :: norder

        !in PARSEC the norder is half the stencil
        mini_grid%norder = norder/2

    end subroutine set_norder

    subroutine set_step(mini_grid,step)
        implicit none

        class (mini_grid_data):: mini_grid
        real(dp), intent (in) :: step

        mini_grid%stepin = step
        mini_grid%step(:) = step
        mini_grid%n1 = int(two*mini_grid%rmax/step) + 2
        mini_grid%n2 = mini_grid%n1
        mini_grid%n3 = mini_grid%n1
        !TODO: print initial cube grid dimensions

        mini_grid%hcub = step**3
        mini_grid%hcub2 = two/(step**3)
    end subroutine set_step

    subroutine set_index(mini_grid,nxmin,nxmax,nymin,nymax,nzmin,nzmax,ierr)
        implicit none

        class (mini_grid_data):: mini_grid
        integer, intent (in) :: nxmin,nxmax,nymin,nymax,nzmin,nzmax
        integer, intent (inout) :: ierr

        integer :: alccheck, nxyz

        mini_grid%nxmax = nxmax
        mini_grid%nxmin = nxmin
        mini_grid%nymax = nymax
        mini_grid%nymin = nymin
        mini_grid%nzmax = nzmax
        mini_grid%nzmin = nzmin

        nxyz = (nxmax-nxmin+1)*(nymax-nymin+1)*(nzmax-nzmin+1)
        mini_grid%nxyz = nxyz

        allocate (mini_grid%indexg (nxmin:nxmax,nymin:nymax,nzmin:nzmax),stat=alccheck)
        if (alccheck /= 0) then 
            ierr = ierr + 2000
            ! ERROR CODE 2000 :: MEMORY ALLOCATION FAILURE FOR INDEXG
            return
        endif
        mini_grid%indexg(:,:,:) = 0

        allocate (mini_grid%indexw (nxmin:nxmax,nymin:nymax,nzmin:nzmax),stat=alccheck)
        if (alccheck /= 0) then 
            ierr = ierr + 2001
            ! ERROR CODE 2001 :: MEMORY ALLOCATION FAILURE FOR INDEXW
            return
        endif
        mini_grid%indexw(:,:,:) = 0

    end subroutine set_index

    subroutine set_wedge(mini_grid,nnodes,ierr)
        implicit none
        class (mini_grid_data):: mini_grid
        integer, intent (in) :: nnodes
        integer, intent (inout) :: ierr
        integer :: alccheck(3),ii

        alccheck(:) = 0
        allocate (mini_grid%kx (mini_grid%nwedge),stat=alccheck(1)) 
        allocate (mini_grid%ky (mini_grid%nwedge),stat=alccheck(2)) 
        allocate (mini_grid%kz (mini_grid%nwedge),stat=alccheck(3)) 
        do ii=1,3
            if (alccheck(ii) /= 0) then
                if (ierr < 2000) ierr=ierr+2000
                ierr = ierr+10
                ! ERRORCODE 20X0 :: MEMORY ALLOCATION FAILURE FOR 3D-1D CONVERSION MAP
            endif
        end do
        allocate(mini_grid%ist(0:nnodes))
    end subroutine set_wedge

    subroutine set_ist(mini_grid,nnodes)
        implicit none
        class (mini_grid_data):: mini_grid
        integer, intent (in) :: nnodes
        integer ii, isize, remainder

        mini_grid%ist(:) = 0

        isize = mini_grid%nwedge / nnodes
        remainder = mini_grid%nwedge - isize*nnodes
        ! distribute (mini_grid%nwedge) points so that the first PEs receive
        ! (isize+1) points and the last ones receive (isize)
        mini_grid%ist(0) = 1
        do ii = 1, remainder
            mini_grid%ist(ii) = mini_grid%ist(ii-1) + isize + 1
        enddo
        do ii = remainder+1,nnodes-1
            mini_grid%ist(ii) = mini_grid%ist(ii-1) + isize
        enddo
            mini_grid%ist(nnodes) = mini_grid%nwedge + 1

    end subroutine set_ist

    subroutine set_ndim(mini_grid,ndim,ierr)
        implicit none
        class (mini_grid_data):: mini_grid
        integer, intent (in) :: ndim
        integer, intent (inout) :: ierr
        integer :: alccheck(5),ii

        alccheck(:) = 0
        mini_grid%ndim = ndim
        !TODO...
        allocate (mini_grid%rindex (ndim+1),stat=alccheck(1))

        allocate (mini_grid%fx(ndim),stat=alccheck(3))
        allocate (mini_grid%fy(ndim),stat=alccheck(4))
        allocate (mini_grid%fz(ndim),stat=alccheck(5))
        do ii=1,5
            if (alccheck(ii) /= 0) then
                if (ierr < 2000) ierr=ierr+2000
                ierr = ierr+1
                ! ERRORCODE 200X :: MEMORY ALLOCATION FAILURE FOR SPARSE REPRESENTATION
            endif
        end do
    end subroutine set_ndim


    subroutine mock_lap_setup(mini_grid)
        implicit none
        class (mini_grid_data):: mini_grid

        if (associated(mini_grid%coeff2_1D)) deallocate(mini_grid%coeff2_1D)
        call mock_lap_coeffs(mini_grid%norder,mini_grid%coeff2_1D)

    end subroutine mock_lap_setup

    ! or function
    subroutine mock_lap_coeffs(norder,coeff_pntr)
        use mini_constants
        use mini_common_mod
        implicit none

        integer, intent(in) :: norder
        real(dp), dimension (:), pointer,intent(out) :: coeff_pntr

        integer :: info


        allocate(coeff_pntr(LAPDIR*norder))
        coeff_pntr(:) = zero
        call my_random_array(coeff_pntr,LAPDIR*norder,LAPDIR*norder,1)

    end subroutine mock_lap_coeffs


end module mini_grid_mod
