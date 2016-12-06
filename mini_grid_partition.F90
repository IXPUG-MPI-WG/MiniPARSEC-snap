!proper header, this is based on parsec etc.
subroutine mini_grid_partition(ierr)
    use mini_constants
    use mini_global_data_mod
    implicit none

    integer, intent (inout) :: ierr

    integer, dimension(1):: ierrvec
    !  cutoff radius used in the determination of grid points
    real(dp) :: rcut
    !  number of neighbors used for derivation (on one side)
    integer nord2
    !  actual size of hamiltonian, dimension of grid arrays
    integer ndim,ldn
    !  bounds of grid along each direction
    integer nxmin,nxmax,nymin,nymax,nzmin,nzmax
    !  variable for distance between points, maximum distance
    real(dp) :: dist,dmax
    !  temporary 1-d index of each point
    integer nef, neibs_num
    !  wrapping indexes
    integer iwrap, iwrap2, iwrap3
    !  counters
    integer ii,jj,i,j,k,j00,k00,ny1,ny2,nz1,nz2,iold
    !  variables for the irreducible wedge
    integer itrans, ir(3), itest, nwedge
    !  Cartesian coordinates of a grid point and its equivalent one
    real(dp) :: rr(3), rw(3), rstep
    !  arrays for irreducible wedge
    integer, allocatable :: ninv(:),rinv(:)

    !  number of procs available
    integer nnodes
    !  processor rank
    integer inode
    !  Indexing of distributed blocks
    integer irows(0:mini_parallel%group_size)
    !
    !  External function:
    !
    integer, external :: Part
    logical, external :: outside_domain

    ierr = 0
    nnodes = mini_parallel%group_size
    nord2 = mini_grid%norder !problem here, see below
    !
    !  PART1: Setup the GRID
    !
    !  determine minimal number of mesh points (per each axis, per
    !  each side of origin) that would encompass a sphere of size rmax
    !  (confined system), or a box of side rmax (periodic system)
    !  inside which is a grid with a spacing of h. Note that for a
    !  periodic system rmax is the WHOLE side of the box, not half of
    !  it, so it's really comparable to the diameter of the sphere!

    nxmin = -mini_grid%n1/2
    nxmax = mini_grid%n1 + nxmin -1
    nymin = -mini_grid%n2/2
    nymax = mini_grid%n2 + nymin -1
    nzmin = -mini_grid%n3/2
    nzmax = mini_grid%n3 + nzmin -1

    ! We will be creating indexg and indexw with padding (2*nord2):
    nxmax = nxmax + nord2
    nxmin = nxmin - nord2
    nymax = nymax + nord2
    nymin = nymin - nord2
    nzmax = nzmax + nord2
    nzmin = nzmin - nord2

    call mini_grid%set_index(nxmin,nxmax,nymin,nymax,nzmin,nzmax,ierr)
#ifdef DEBUG
    if (mini_parallel%verbosity > NOTXT) then
        write(7,*) 'DEBUG : nxyz after set_index are:',mini_grid%nxyz
        write(7,*) 'DEBUG : nxmax after set_index are:',mini_grid%nxmax
    endif
#endif
    if (ierr > 0 ) return

    ! Resort to the orignal number without padding
    nxmin = -mini_grid%n1/2
    nxmax = mini_grid%n1 + nxmin -1
    nymin = -mini_grid%n2/2
    nymax = mini_grid%n2 + nymin -1
    nzmin = -mini_grid%n3/2
    nzmax = mini_grid%n3 + nzmin -1

    !
    !  Assign spatial coordinates (x,y,z) to each array point (i,j,k).
    !  To each point, assign a 1d index (the counter is
    !  nef). Build the 3d->1d conversion arrays, kx, ky, kz. The
    !  actual position (atomic units, cartesian coordinates, relative
    !  to the origin) of each grid point is: xx=h*(kx+shift),
    !  yy=h*(ky+shift), zz=h*(kz+shift).
    !
    !  start by initializing the 3d -> 1d conversion index
    !TODO this updates indexg in the wrong order. should change loop to k,j,i

    rcut = mini_grid%rmax * mini_grid%rmax
    nef = 0
    do i = nxmax, nxmin, -1
        rr(1) = (mini_grid%shift(1) + i)*mini_grid%step(1)
        do j = nymax, nymin, -1
            rr(2) = (mini_grid%shift(2) + j)*mini_grid%step(2)
            do k = nzmax, nzmin, -1
                rr(3) = (mini_grid%shift(3) + k)*mini_grid%step(3)
            ! cluster calculation, check the point against the domain
            if(outside_domain(mini_grid, rr)) cycle

            nef = nef + 1
            mini_grid%indexg(i,j,k) = nef
            enddo
        enddo
    enddo
    !
    !  The total number of points inside the sphere is the size of the
    !  Hamiltonian of the system, ndim. Allocate arrays accordingly.
    !
    ndim = nef
    call mini_grid%set_ndim(ndim,ierr)
#ifdef DEBUG
    if (mini_parallel%verbosity > NOTXT) then
        write(7,*) 'DEBUG : ndim after set_ndim is:',ndim

    endif
#endif
    if (ierr > 0 ) return

    !
    !  Define the set of grid points in the irreducible wedge, IW.
    !  The IW contains all grid points that can reconstruct the
    !  original grid by the application of symmetry operations, and it
    !  is such that no two points in the IW are equivalent to each
    !  other by the application of any symmetry operation.
    !
    !  Information about the wedge is contained in arrays indexw, rindex:
    !  the n-th point in the IW has coordinates i,j,k so that
    !  indexw(i,j,k)=n. If indexw(i,j,k)=0, this point is not in the IW
    !
    !  if rindex(n)=m and rtrans(n)=i, the n-th point in original grid
    !  is equivalent to the m-th point in IW upon symmetry operation i.
    !
    !  dmax stores the maximum distance between a grid point and its
    !  equivalent in the wedge after a symmetry operation. Ideally,
    !  these two point should coincide, but numerical roundoff or the
    !  fractional translation symm%tnp (see subroutine symop) may interfere.
    dmax = zero
    nwedge = 0

    !TODO this updates indexg in the wrong order. should change loop to k,j,i
    do i = nxmax, nxmin, -1
      rr(1) = (mini_grid%shift(1) + i)*mini_grid%step(1)
      do j = nymax, nymin, -1
        rr(2) = (mini_grid%shift(2) + j)*mini_grid%step(2)
        do k = nzmax, nzmin, -1
          rr(3) = (mini_grid%shift(3) + k)*mini_grid%step(3)
          if (mini_grid%indexg(i,j,k) == 0) cycle

            nwedge = nwedge + 1
            mini_grid%indexw(i,j,k) = nwedge
            itest = nwedge
            mini_grid%rindex(mini_grid%indexg(i,j,k)) = itest
        enddo
      enddo
    enddo

    mini_grid%nwedge = nwedge

    call mini_grid%set_wedge(mini_parallel%procs_num,ierr)
    if (ierr > 0 ) return
 
    !AJB: TODO: ldn should be divisible by the AVX register size (MKL: just by 16)
    !AJB: TODO (2) : MKL says ldn should NOT be divisible by 2048
    ldn = nwedge/nnodes + nnodes
    mini_parallel%ldn = ldn
    mini_parallel%ndim = ndim
    mini_parallel%nwedge = nwedge

    if (mini_parallel%verbosity == DEBUGEACH) then
          write(9,*)
          write(9,*) 'MASTER PROC - Setup messages:'
          write(9,*) '---------------'
          write(9,*) 'The effective matrix size is', nwedge
          ! write(9,*) 'The effective matrix size is', ndim
          ! write(9,*) ' reduced size is ',nwedge
          write(9,*) ' the ldn length is ',ldn
          write(9,*)
          write(9,15) LAPDIR*nord2
    endif

    if (mini_parallel%verbosity > NOTXT) then
        if (mini_parallel%iammaster) then
            write(7,*)
            write(7,*) 'Setup messages:'
            write(7,*) '---------------'
            write(7,*) 'The effective matrix size is', nwedge
            ! write(7,*) 'The effective matrix size is', ndim
            ! write(7,*) ' reduced size is ',nwedge
            write(7,*)
            write(7,15) LAPDIR*nord2
        endif
    endif

15 format(1x,'There are',1x,i2,1x, 'laplacian-related non-diagonal elements per row') 

  !
  !  PART2: Partitioning The Domain
  !
  !  At this point we are ready to do any kind of partitioning
  !  geometrical or based on the nwedge-array.
  !  Part(i,j,k,..) = processor number that the (i,j,k) resides on
  !

  irows(:) = 0

  ! only master process here, rest do it later 
  call mini_grid%set_ist(nnodes)

  do i = nxmax, nxmin, -1
     do j = nymax, nymin, -1
        do k = nzmax, nzmin, -1
           if (mini_grid%indexg(i,j,k) /= 0) then
              ii = mini_grid%rindex(mini_grid%indexg(i,j,k))
              mini_grid%fx(mini_grid%indexg(i,j,k)) = i
              mini_grid%fy(mini_grid%indexg(i,j,k)) = j
              mini_grid%fz(mini_grid%indexg(i,j,k)) = k
           endif
           if (mini_grid%indexw(i,j,k) /= 0) then
              irows(Part(mini_grid,nnodes,i,j,k,.true.)) &
                   =irows(Part(mini_grid,nnodes,i,j,k,.true.)) + 1
           endif

        enddo
     enddo
  enddo

  !  Transform irows to show where the first row of each node
  !  begins. I.e., 3 nodes: first has 10 rows, second has 5 rows,
  !  third 6 rows irows(0) = 1, irows(1) = 11, irows(2) = 16,
  !  irows(3) = 22.
  !  Used for permuting the rows later.
  !
  irows(nnodes) = nwedge + 1
  do i = nnodes-1, 0, -1
     irows(i) = irows(i+1) - irows(i)
  enddo

  do i = nxmin,nxmax
   do j = nymin, nymax
    do k = nzmin, nzmax
      if (mini_grid%indexw(i,j,k) /= 0) then
          iold  = mini_grid%indexw(i,j,k)
          inode = Part(mini_grid,nnodes,i,j,k,.true.)
          nef   = irows( inode )
          mini_grid%indexw(i,j,k) = nef
          irows( inode )= irows( inode ) + 1
          !AJB: this is the only place where k[xyz] gets updated
          mini_grid%kx(nef)  = i
          mini_grid%ky(nef)  = j
          mini_grid%kz(nef)  = k
      endif
     enddo
    enddo
  enddo
! #ifdef DEBUG
! write(9,*) 'irows in partition', irows
! write(9,*) 'nnodes', nnodes
! #endif

  !
  !  Restoring irows so that they point at the first row of each
  !  processor
  !
  irows(nnodes) = mini_grid%nwedge + 1
  do i = nnodes -1, 1,-1
     irows(i) = irows(i-1)
  enddo
  irows(0) = 1
  !
  !  Transfer irows to parallel structure
  !
  mini_parallel%irows = irows

    ! Set the indexes for grid points outside the domain to be nwedge+1 (it was
    ! zero before). This is relevant along non-periodic directions only.
    !
    ! AJB: TODO: since the size of the workverctor would vary according to the actual data received,
    ! eventually we would like to have the "out of scope" cell to be index 0 
    !
    ! AJB2: Actually, the workvector(s!) should have dimension mydim, and the recieved data should be kept
    ! seperately, so "out of scope" has two meanings. one "out of scope" is explicitly nwedge+1,
    ! and a local "out of scope" meaning  ">parallel%mydim". This should be setup later.

  do k = nzmin - nord2,nzmax + nord2
     do j = nymin - nord2,nymax + nord2
        do i = nxmin - nord2,nxmax + nord2
           if (mini_grid%indexw(i,j,k) == 0) mini_grid%indexw(i,j,k) = mini_grid%nwedge + 1
           if (mini_grid%indexg(i,j,k) == 0) mini_grid%indexg(i,j,k) = mini_grid%ndim + 1
        enddo
     enddo
  enddo

  mini_grid%rindex(mini_grid%ndim+1) = mini_grid%nwedge + 1

end subroutine mini_grid_partition
!===============================================================
!
!  Simple partitioning routine.
!  Split the nwedge-long vector into nnodes equal parts.
!
!---------------------------------------------------------------
integer function Part(mini_grid,nnodes,i,j,k,lprint)

  use mini_constants
  use mini_grid_mod
  implicit none
  !
  !  Input/Output variables:
  !
  !AJB: why inout?am I missing something here?
  type (mini_grid_data), intent(inout) :: mini_grid
  integer, intent(in) :: nnodes
  integer, intent(in) :: i,j,k
  logical, intent(in) :: lprint

  integer nef, ii
!---------------------------------------------------------------

  !  Give the processor number for the current i,j,k
  !  For this routine partition according to rows:
  !  (ijk) -> Indexw(ijk)
  !
  nef = mini_grid%indexw(i,j,k)
  do ii = 0, nnodes-1
     if (nef >= mini_grid%ist(ii) .and. nef < mini_grid%ist(ii+1)) then
        Part = ii
        return
     endif
  enddo

  Part = -1

  !if (lprint) write(7,*) 'The row nef:',nef &
  !     ,'is outside bounds nwedge:',mini_grid%nwedge

end function Part
!===============================================================
!
!  Determines if the point rr is outside the computational
!  domain.  The grid structure is needed since it contains
!  the type of domain shape and the necessary parameters.
!
!---------------------------------------------------------------
logical function outside_domain(mini_grid,rr)

  use mini_constants
  use mini_grid_mod
  implicit none
  !
  !  Input/Output variables:
  !
  type (mini_grid_data), intent(in) :: mini_grid
  real(dp), intent(in) :: rr(3)
  !
  !  Work variables:
  !
  real(dp) :: rad2, rr_rad2, coord, polar1, polar2
  !---------------------------------------------------------------

  ! generic return -> point rr is outside the domain
  outside_domain = .true.
  
  select case(mini_grid%domain_shape)
  case (0)
     ! Sphere
     ! points inside the sphere obey 
     ! x**2 + y**2 + z**2 <= r**2

     ! grid%d_shape_param(1) is the radius
     rad2 = mini_grid%d_shape_param(1)**2

     rr_rad2 = dot_product(rr, rr)

     outside_domain = (rr_rad2 > rad2)
  case (1)
     ! Ellipsoid
     ! points inside the sphere obey 
     ! (x/r_x)**2 + (y/r_y)**2 + (z/r_z)**2 <= 1

     outside_domain = &
       ((rr(1) / mini_grid%d_shape_param(1))**2 + &
        (rr(2) / mini_grid%d_shape_param(2))**2 + &
        (rr(3) / mini_grid%d_shape_param(3))**2     > 1.0d0)
  case (2)
     ! Cylinder
     ! points inside the cylinder obey
     ! x**2 + y**2 <= r**2 AND |z| <= length
     ! for a cylinder oriented along the z-axis

     ! Get the length coordinate and the polar coordinates
     select case(mini_grid%i_shape_param(1))
         case (1)
            ! oriented along x-axis
            coord = rr(1)
            polar1 = rr(2)
            polar2 = rr(3)
         case (2)
            ! oriented along y-axis
            coord = rr(2)
            polar1 = rr(1)
            polar2 = rr(3)
         case (3)
            ! oriented along z-axis
            coord = rr(3)
            polar1 = rr(2)
            polar2 = rr(1)
     end select

     outside_domain = &
       ( (polar1**2 + polar2**2 >  mini_grid%d_shape_param(1)**2   ) &
                  .or. coord**2 > (mini_grid%d_shape_param(2)/2.0d0)**2 )
  case (3)
     ! Box
     ! If the coordinate is outside the length/2 (in magnitude)...

     outside_domain = (rr(1)**2 > (mini_grid%d_shape_param(1)/2.0d0)**2 .or. &
        rr(2)**2 > (mini_grid%d_shape_param(2)/2.0d0)**2 .or. &
        rr(3)**2 > (mini_grid%d_shape_param(3)/2.0d0)**2)
  end select

end function outside_domain
!===============================================================
