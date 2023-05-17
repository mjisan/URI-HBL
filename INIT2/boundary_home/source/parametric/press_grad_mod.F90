module press_grad_mod

use shr_kind_mod,  only : r4 => shr_kind_r4
use shr_kind_mod,  only : r8 => shr_kind_r8

use grid_mod,      only : nx, ny
use grid_mod,      only : stepx, stepy
use constants_mod, only : fcor

implicit none

private

#include 'grid_params_r4.h'

public :: press_grad

contains

subroutine press_grad (um, vm, pgfx, pgfy)

  real(KIND=r4), dimension(nx,ny),intent(in)  :: um
  real(KIND=r4), dimension(nx,ny),intent(in)  :: vm
  real(KIND=r4), dimension(nx,ny),intent(out) :: pgfx
  real(KIND=r4), dimension(nx,ny),intent(out) :: pgfy
  
  real(KIND=r4), dimension(nx,ny) :: coriolis

  real(KIND=r4) :: difu, difv
  integer :: i, j
  
  pgfx = 0.0_r4
  pgfy = 0.0_r4

  do j=1,ny
    do i=1,nx
      coriolis(i,j) = fcor
    enddo
  enddo

  do j=2,ny-1
   do i=2,nx-1
      if (um(i,j)  >= 0.0_r4) then
         difu = vm(i,j)-vm(i-1,j)
      else
         difu = vm(i+1,j)-vm(i,j)
      endif
      if (vm(i,j) >= 0.0_r4) then
         difv = vm(i,j)-vm(i,j-1)
      else
         difv = vm(i,j+1)-vm(i,j)
      endif
      pgfy(i,j) = um(i,j)*difu*stepxinv +    &
                  vm(i,j)*difv*stepyinv + coriolis(i,j)*um(i,j)

      if (um(i,j) >= 0.0_r4) then
         difu = um(i,j)-um(i-1,j)
      else
         difu = um(i+1,j)-um(i,j)
      endif
      if (vm(i,j) >= 0.0_r4) then
         difv = um(i,j)-um(i,j-1)
      else
         difv = um(i,j+1)-um(i,j)
      endif
      pgfx(i,j) = um(i,j)*difu*stepxinv +    &
                  vm(i,j)*difv*stepyinv - coriolis(i,j)*vm(i,j)
   enddo
  enddo

end subroutine press_grad

end module press_grad_mod
