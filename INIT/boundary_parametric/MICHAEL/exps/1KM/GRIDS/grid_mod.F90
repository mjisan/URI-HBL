module grid_mod

use shr_kind_mod,  only : r4 => shr_kind_r4
use shr_kind_mod,  only : r8 => shr_kind_r8


   integer, parameter :: nx = 1801
   integer, parameter :: ny = 1801
   integer, parameter :: nz = 101

   integer, parameter :: npx = 1
   integer, parameter :: npy = 1

   integer, parameter :: icen = nx/2 + 1
   integer, parameter :: jcen = ny/2 + 1

   integer, parameter :: rpts = 151

   real(KIND=r4), parameter :: stepx = 1.0E03 
   real(KIND=r4), parameter :: stepy = 1.0E03 

end module grid_mod
