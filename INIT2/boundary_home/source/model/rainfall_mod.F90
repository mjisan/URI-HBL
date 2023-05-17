module rainfall_mod

use shr_kind_mod, only : r8 => shr_kind_r8

use grid_mod, only : nx, ny
use time_mod, only : delt

implicit none

private

public :: rainfall

real(KIND=r8), parameter :: a1 = -1.10, a2 = -1.60, a3 =  64.5, a4 = 150.0
real(KIND=r8), parameter :: b1 =  3.96, b2 =  4.80, b3 = -13.0, b4 = -16.0

logical :: first_time = .true.

contains

subroutine rainfall (vmax, rmax, rad, rain_acc)

  real(KIND=r8),                  intent(in)    :: vmax
  real(KIND=r8),                  intent(in)    :: rmax
  real(KIND=r8), dimension(nx,ny),intent(in)    :: rad 
  real(KIND=r8), dimension(nx,ny),intent(inout) :: rain_acc

  real(KIND=r8) :: rain, vmax1, rmax1, rad1
  real(KIND=r8) :: t0, tm, rm, re, u

  integer :: i, j

  if ( first_time ) then
    call rainfall_init (rain_acc)
    first_time = .false.
  endif

  do j=1,ny
    do i=1,nx
      vmax1 = vmax*1.94384     ! knots
      rmax1 = rmax*1.0E-03     ! km
      rad1  = rad(i,j)*1.0E-03 ! km

      u  = 1.0 + (vmax1 - 35.0)/33.0
      t0 = a1 + b1*u
      tm = a2 + b2*u
      rm = a3 + b3*u
      re = a4 + b4*u

      if (rad1 .lt. rm) then
        rain = t0 + (tm -t0)*(rad1/rm)
      else
        rain = tm*exp(-1.0*(rad1 - rm)/re)
      endif 

      rain_acc(i,j) = rain_acc(i,j) + delt*rain

    enddo
  enddo

end subroutine rainfall

subroutine rainfall_init ( rain_acc )

  real(KIND=r8), dimension(nx,ny),intent(inout) :: rain_acc

  integer :: i, j

  do j=1,ny 
    do i=1,nx 
      rain_acc(i,j) = 0.0
    enddo
  enddo

end subroutine rainfall_init

end module rainfall_mod
