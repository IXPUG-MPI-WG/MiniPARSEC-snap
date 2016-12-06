!===============================================================
!
!  Copyright (C) 2005 Finite Difference Research Group
!  This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
!  NUMERICAL CONSTANTS
!
!---------------------------------------------------------------
module mini_constants

  implicit none
  !
  !  Definitions of double precision real/complex types:
  !
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: dpc = kind((1.0d0,1.0d0))
  !  Machine precision for 64bit 2**-52
  real(dp), parameter :: mach_eps    = 2.220446049250313080847263336182d-16
  real(dp), parameter :: twomach_eps = 4.440892098500626161694526672363d-16
  !
  ! SCENARIO PARAMETERS - CHANGE HERE
  !
  !number of comm. buffers to populate (ring buffer, min 2)
  integer, parameter :: NUM_BUFFER = 2
  !number of internal kernel runs, increase to show inter-node sync problems
  integer, parameter :: NUM_WAVEFUN = 10
  !number of outer iterations, increase for statistics
  integer, parameter :: NUM_OUTER = 40
  !number of iterations outside of kernel (generally don't touch)
  integer, parameter :: OUTSIDE_WORK = 1
  !CLI variables - defaults
  real(dp), parameter :: DEF_RMAX = 45.0d0
  real(dp), parameter :: DEF_STEP = 0.3d0
  integer, parameter :: DEF_ORDER = 6
  !
  !  Numerical constants:
  !
  real(dp), parameter :: zero = 0.d0, one = 1.d0, &
       two = 2.d0, three = 3.d0, four = 4.d0, five = 5.d0, &
       six = 6.d0, seven = 7.d0, eight = 8.d0, nine = 9.d0, mone = -1.d0, &
       half = 0.5d0, third = 1.d0/3.d0, mhalf = -0.5d0, &
       root3 = 1.73205080756887729352744634151d0, &
       pi = 3.1415926535897932384626433832795d0, &
       twopi = two * pi
  complex(dpc), parameter :: zzero = (0.d0,0.d0), &
       zone = (1.d0,0.d0), zi = (0.d0,1.d0)
  !
  !  Physical constants:
  !
  !  convert energy units from rydbergs to eV
  real(dp), parameter :: rydberg = 13.6058d0
  !  convert length units from atomic units (bohrs) to angstroms
  real(dp), parameter :: angs = 0.529177d0, cubangs = angs**3
  !  convert dipole units from au to debyes
  real(dp), parameter :: debye = 2.54d0
  !  convert temperature from kelvins to rydbergs
  real(dp), parameter :: tempkry = 6.33327186d-06
  !  convert time from a.u. to ps
  real(dp), parameter :: timeaups = 2.4188843d-05
  !
  !
  !PARAMETERS
  !
  integer, parameter :: NOTXT = 0, MINIMAL = 1, DEBUGEACH = 2
  integer, parameter :: LAPDIR = 6

end module mini_constants
!===============================================================
