module mpp_grid_mod

   use shr_kind_mod, only : r8 => shr_kind_r8

   implicit none
   private

   include 'mpif.h'

   public :: mpp_get_global_indices

contains

subroutine mpp_get_global_indices ( xbeg_data, xend_data, ybeg_data, yend_data, &
                                    xbeg_comp, xend_comp, ybeg_comp, yend_comp, &
                                                 id, nx, ny, ngx, ngy, npx, npy )

!  computes beginning and ending global indices on data domain for each mpi-rank

   integer, intent(out) :: xbeg_data, xend_data, ybeg_data, yend_data
   integer, intent(out) :: xbeg_comp, xend_comp, ybeg_comp, yend_comp
   integer, intent(in)  :: id, nx, ny, ngx, ngy, npx, npy

   integer :: xid, yid

!  computes beginning and ending global indices on data domain for each mpi-rank

   if ( id <= npx ) then
     xid = id
     yid = 1
   else
     yid = (id-1)/npx + 1
     xid = id - (yid-1)*npx
   endif

   xbeg_data = (nx-2)*(xid-1) + 1
   xend_data = xbeg_data + nx - 1

   ybeg_data = (ny-2)*(yid-1) + 1
   yend_data = ybeg_data + ny - 1

!  computes beginning and ending global indices on compute domain for each mpi-rank

   xbeg_comp = xbeg_data + 1
   xend_comp = xend_data - 1
   ybeg_comp = ybeg_data + 1
   yend_comp = yend_data - 1

   if ( xbeg_data == 1   ) xbeg_comp = xbeg_data
   if ( xend_data == ngx ) xend_comp = xend_data
   if ( ybeg_data == 1   ) ybeg_comp = ybeg_data
   if ( yend_data == ngy ) yend_comp = yend_data

end subroutine mpp_get_global_indices

end module mpp_grid_mod
