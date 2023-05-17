module mpp_io_boundary_model_mod

use shr_kind_mod,          only : r8 => shr_kind_r8
use shr_kind_mod,          only : r4 => shr_kind_r4

use mpp_io_boundary_model_spatial_mod,  only : write_boundary_model_spatial
use mpp_io_boundary_model_temporal_mod, only : write_boundary_model_temporal
   implicit none
   private

   public :: write_boundary_model

contains
  
subroutine write_boundary_model ( um, vm, pgfx, pgfy, wm, tur, znot, mask, &
                  rain_acc, wsmax, x_center, y_center, x_corner, y_corner, &
                                      recs, nx, ny, nz, ngx, ngy, npx, npy )

   integer, intent(in) :: nx, ny, nz
   integer, intent(in) :: ngx, ngy
   integer, intent(in) :: npx, npy

   real(KIND=r8), dimension(:,:,:), intent(in)    :: um, vm
   real(KIND=r8), dimension(:,:),   intent(inout) :: pgfx, pgfy
   real(KIND=r8), dimension(:,:,:), intent(inout) :: wm, tur
   real(KIND=r8), dimension(:,:),   intent(inout) :: znot, mask
   real(KIND=r8), dimension(:,:),   intent(inout) :: rain_acc
   real(KIND=r8),                   intent(in)    :: wsmax
   real(KIND=r8),                   intent(in)    :: x_center, y_center
   real(KIND=r8), dimension(:),     intent(in)    :: x_corner, y_corner
   integer,                         intent(in)    :: recs
 
   call write_boundary_model_spatial ( um, vm, pgfx, pgfy, wm, tur, znot, mask, &
                                 rain_acc, recs, nx, ny, nz, ngx, ngy, npx, npy )

   call write_boundary_model_temporal ( wsmax, x_center, y_center, &
                                          x_corner, y_corner, recs )
 
end subroutine write_boundary_model

end module mpp_io_boundary_model_mod
