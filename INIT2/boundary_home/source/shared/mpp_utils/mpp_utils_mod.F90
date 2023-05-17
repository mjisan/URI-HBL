module mpp_utils_mod

use shr_kind_mod,  only : r4 => shr_kind_r4
use shr_kind_mod,  only : r8 => shr_kind_r8

use mpp_setup_mod, only : size, master, rank
use mpp_io_mod,    only : stdout, error_handler
use mpp_setup_mod, only : npx, npy

   implicit none
   include 'mpif.h'

   private

   interface mpp_max
      module procedure mpp_max_3d
      module procedure mpp_max_2d
      module procedure mpp_max_1d
      module procedure mpp_max_0d
   end interface

   interface mpp_min
      module procedure mpp_min_3d
      module procedure mpp_min_2d
      module procedure mpp_min_1d
      module procedure mpp_min_0d
   end interface

   interface mpp_bcast
      module procedure mpp_bcast_0d
   end interface

   public :: mpp_max
   public :: mpp_min
   public :: mpp_bcast

contains

subroutine mpp_max_3d (var_local, max_var_global)

   real(KIND=r8), dimension(:,:,:), intent(in)  :: var_local
   real(KIND=r8),                   intent(out) :: max_var_global

   integer       :: ierr
   real(KIND=r8) :: var_max_local
   
   var_max_local = MAXVAL(var_local)

   if (npx*npy == 1 ) then
     max_var_global = var_max_local
   else
     call MPI_REDUCE (var_max_local, max_var_global, 1, MPI_DOUBLE_PRECISION, &
                                         MPI_MAX, master, MPI_COMM_WORLD, ierr)
   endif

end subroutine mpp_max_3d

subroutine mpp_max_2d (var_local, max_var_global)

   real(KIND=r8), dimension(:,:), intent(in)  :: var_local
   real(KIND=r8),                 intent(out) :: max_var_global

   integer       :: ierr
   real(KIND=r8) :: var_max_local
   
   var_max_local = MAXVAL(var_local)

   if (npx*npy == 1 ) then
     max_var_global = var_max_local
   else
     call MPI_REDUCE (var_max_local, max_var_global, 1, MPI_DOUBLE_PRECISION, &
                                         MPI_MAX, master, MPI_COMM_WORLD, ierr)
   endif

end subroutine mpp_max_2d

subroutine mpp_max_1d (var_local, max_var_global)

   real(KIND=r8), dimension(:), intent(in)  :: var_local
   real(KIND=r8),               intent(out) :: max_var_global

   integer       :: ierr
   real(KIND=r8) :: var_max_local
   
   var_max_local = MAXVAL(var_local)

   if (npx*npy == 1 ) then
     max_var_global = var_max_local
   else
     call MPI_REDUCE (var_max_local, max_var_global, 1, MPI_DOUBLE_PRECISION, &
                                         MPI_MAX, master, MPI_COMM_WORLD, ierr)
   endif

end subroutine mpp_max_1d
 
subroutine mpp_max_0d (var_local, max_var_global)

   real(KIND=r8), intent(in)  :: var_local
   real(KIND=r8), intent(out) :: max_var_global

   integer       :: ierr
   real(KIND=r8) :: var_max_local
   
   var_max_local = var_local

   if (npx*npy == 1 ) then
     max_var_global = var_max_local
   else
     call MPI_REDUCE (var_max_local, max_var_global, 1, MPI_DOUBLE_PRECISION, &
                                         MPI_MAX, master, MPI_COMM_WORLD, ierr)
   endif

end subroutine mpp_max_0d

subroutine mpp_min_3d (var_local, min_var_global)

   real(KIND=r8), dimension(:,:,:), intent(in)  :: var_local
   real(KIND=r8),                   intent(out) :: min_var_global

   integer       :: ierr
   real(KIND=r8) :: var_min_local
   
   var_min_local = MINVAL(var_local)

   if (npx*npy == 1 ) then
     min_var_global = var_min_local
   else
     call MPI_REDUCE (var_min_local, min_var_global, 1, MPI_DOUBLE_PRECISION, &
                                         MPI_MIN, master, MPI_COMM_WORLD, ierr)
   endif

end subroutine mpp_min_3d

subroutine mpp_min_2d (var_local, min_var_global)

   real(KIND=r8), dimension(:,:), intent(in)  :: var_local
   real(KIND=r8),                 intent(out) :: min_var_global

   integer       :: ierr
   real(KIND=r8) :: var_min_local
   
   var_min_local = MINVAL(var_local)

   if (npx*npy == 1 ) then
     min_var_global = var_min_local
   else
     call MPI_REDUCE (var_min_local, min_var_global, 1, MPI_DOUBLE_PRECISION, &
                                         MPI_MIN, master, MPI_COMM_WORLD, ierr)
   endif

end subroutine mpp_min_2d

subroutine mpp_min_1d (var_local, min_var_global)

   real(KIND=r8), dimension(:), intent(in)  :: var_local
   real(KIND=r8),               intent(out) :: min_var_global

   integer       :: ierr
   real(KIND=r8) :: var_min_local
   
   var_min_local = MINVAL(var_local)

   if (npx*npy == 1 ) then
     min_var_global = var_min_local
   else
     call MPI_REDUCE (var_min_local, min_var_global, 1, MPI_DOUBLE_PRECISION, &
                                         MPI_MIN, master, MPI_COMM_WORLD, ierr)
   endif

end subroutine mpp_min_1d

subroutine mpp_min_0d (var_local, min_var_global)

   real(KIND=r8), intent(in)  :: var_local
   real(KIND=r8), intent(out) :: min_var_global

   integer       :: ierr
   real(KIND=r8) :: var_min_local
   
   var_min_local = var_local

   if (npx*npy == 1 ) then
     min_var_global = var_min_local
   else
     call MPI_REDUCE (var_min_local, min_var_global, 1, MPI_DOUBLE_PRECISION, &
                                         MPI_MIN, master, MPI_COMM_WORLD, ierr)
   endif

end subroutine mpp_min_0d

subroutine mpp_bcast_0d (var_local)

   real(KIND=r8), intent(in)  :: var_local
   integer :: ierr

   call MPI_BCAST ( var_local, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr) 

end subroutine mpp_bcast_0d

end module mpp_utils_mod

