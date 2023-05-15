module constants_mod

use shr_kind_mod,  only : r4 => shr_kind_r4
use shr_kind_mod,  only : r8 => shr_kind_r8

   implicit none
!   real(KIND=r8), parameter :: fcor    =  6.575702e-05            ! coriolis parameter (rad/sec)
   real(KIND=r8), parameter :: fcor    =  6.38e-05                 ! coriolis parameter (rad/sec)   
   real(KIND=r8), parameter :: pie     = 4.0_r8*ATAN(1.0_r8)       ! pie
   real(KIND=r8), parameter :: rearth  = 6.37122e+06               ! radius of earth (m)
   real(KIND=r8), parameter :: grav    = 9.80616e+00               ! acceleration of gravity (m/s^2)
   real(KIND=r8), parameter :: rhosw   = (4.1_r8/3.99_r8)*1.0e+03  ! density of salt water (kg/m3)
   real(KIND=r8), parameter :: rhoair  = 1.15e+00                  ! density of air(kg/m^3)

   real(KIND=r8), parameter :: mdmv    = 0.622                     ! molar weight for dry air/
                                                                   ! mean molar weight for moist air
   real(KIND=r8), parameter :: epsilon = (1.0_r8-mdmv)/mdmv 

   real(KIND=r8), parameter :: deg2rad = pie/180.0_r8              ! degrees to radians

end module constants_mod


