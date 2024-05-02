module rolls_grid_mod

use shr_kind_mod,     only : r8 => shr_kind_r8

   implicit none

   integer, parameter :: nsrmx  = 1
   integer, parameter :: nsrmy  = 1

   integer, parameter :: nxsrm = 513
   integer, parameter :: nzsrm = 101

   real(KIND=r8), parameter :: stepx = 30.0
   real(KIND=r8), parameter :: stepz = 30.0

end module rolls_grid_mod
