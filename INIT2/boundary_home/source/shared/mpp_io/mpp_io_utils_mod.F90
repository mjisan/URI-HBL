module mpp_io_utils_mod

use shr_kind_mod,     only : r8 => shr_kind_r8
use shr_kind_mod,     only : r4 => shr_kind_r4

use mpp_setup_mod,    only : size, master, rank
use mpp_io_mod,       only : stdout, error_handler

   implicit none
   private

   include 'netcdf.inc'

   public :: check_file_dims, get_file_dims

   integer :: x_dimid, y_dimid, z_dimid, rec_dimid
   integer :: y_varid, x_varid, z_varid, rec_varid
   integer :: retval, ncid

contains

subroutine check_file_dims (file_name, x_size, x_name, y_size, y_name, z_size, z_name, &
                                                                     rec_size, rec_name)

   character(LEN=*),           intent(in) :: file_name

   integer,                    intent(in) :: x_size
   character(LEN=*),           intent(in) :: x_name
   integer,                    intent(in) :: y_size
   character(LEN=*),           intent(in) :: y_name

   integer,          optional, intent(in) :: z_size
   character(LEN=*), optional, intent(in) :: z_name
   integer,          optional, intent(in) :: rec_size
   character(LEN=*), optional, intent(in) :: rec_name

   integer :: x_dim, y_dim, z_dim, rec_dim

!  open netCDF file

   retval = nf_open(file_name, nf_nowrite, ncid)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'check_file_dims/nf_open')

!  inquiry into information about the dimensions

   if ( present(z_name)) then
   call get_file_dims ( file_name, x_dim, x_name, y_dim, y_name, z_dim=z_dim, &
                                                                z_name=z_name )
   endif

   if ( present(rec_name)) then
   call get_file_dims ( file_name, x_dim, x_name, y_dim, y_name, rec_dim=rec_dim, &
                                                                rec_name=rec_name )
  endif
   
   if ( present(z_name) .and. present(rec_name)) then
   call get_file_dims ( file_name, x_dim, x_name, y_dim, y_name, &
                                     z_dim=z_dim, z_name=z_name, &
                              rec_dim=rec_dim, rec_name=rec_name )
  endif

  if ( x_size /= x_dim ) then
     if ( master == rank ) then
       write(stdout,1000) x_size, x_dim
     endif
1000 format(2x, 'x_size = ', i4, 2x, 'x_dim = ', i4 )

     call flush(stdout); call error_handler ( 'ERROR', 'check_file_dims/x_size .NE. x_dim' )
   endif

   if ( y_size /= y_dim ) then
     if ( master == rank ) then
       write(stdout,2000) y_size, y_dim
     endif
2000 format(2x, 'y_size = ', i4, 2x, 'y_dim = ', i4 )

     call flush(stdout); call error_handler ( 'ERROR', 'check_file_dims/y_size .NE. y_dim' )
   endif

   if ( present(z_name) ) then
     if ( z_size /= z_dim ) then
       if ( master == rank ) then
         write(stdout,3000) z_size, z_dim
       endif
3000 format(2x, 'z_size = ', i4, 2x, 'z_dim = ', i4 )

       call flush(stdout); call error_handler ( 'ERROR', 'check_file_dims/z_size .NE. z_dim' )
     endif
   endif

   if ( present(rec_name) ) then
     if ( rec_size > rec_dim ) then
       if ( master == rank ) then
         write(stdout,4000) rec_size, rec_dim
       endif
4000 format(2x, 'rec_size = ', i4, 2x, 'rec_dim = ', i4 )

       call error_handler ( 'ERROR', 'check_file_dims/rec_size .NE. rec_dim' )
     endif
   endif

end subroutine check_file_dims

subroutine get_file_dims ( file_name, x_dim, x_name, y_dim, y_name, z_dim, z_name, &
                                                                 rec_dim, rec_name )

  character(LEN=*),           intent(in)  :: file_name

  integer,                    intent(out) :: x_dim 
  character(LEN=*),           intent(in)  :: x_name
  integer,                    intent(out) :: y_dim 
  character(LEN=*),           intent(in)  :: y_name

  integer,          optional, intent(out) :: z_dim 
  character(LEN=*), optional, intent(in)  :: z_name
  integer,          optional, intent(out) :: rec_dim 
  character(LEN=*), optional, intent(in)  :: rec_name

!  open netCDF file

   retval = nf_open(file_name, nf_nowrite, ncid)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'check_file_dims/nf_open')

!  inquiry into information about the dimensions

   retval = nf_inq_dimid(ncid, x_name, x_dimid)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'check_file_dims/nf_inq_dimid/x_name' )

   retval = nf_inq_dimlen(ncid, x_dimid, x_dim)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'check_file_dims/nf_inq_dimlen/x_dimid' )

   retval = nf_inq_dimid(ncid, y_name, y_dimid)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'check_file_dims/nf_inq_dimid/y_name' )

   retval = nf_inq_dimlen(ncid, y_dimid, y_dim) 
   if (retval .NE. 0) call error_handler ( 'ERROR', 'check_file_dims/nf_inq_dimlen/y_dimid' )

   if ( present (z_name) ) then
     retval = nf_inq_dimid(ncid, z_name, z_dimid)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'check_file_dims/nf_inq_dimid/z_name' )

     retval = nf_inq_dimlen(ncid, z_dimid, z_dim) 
     if (retval .NE. 0) call error_handler ( 'ERROR', 'check_file_dims/nf_inq_dimlen/z_dimid' )
   endif

   if ( present (rec_name) ) then
     retval = nf_inq_dimid(ncid, rec_name, rec_dimid)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'check_file_dims/nf_inq_dimid/rec_name' )

     retval = nf_inq_dimlen(ncid, rec_dimid, rec_dim) 
   endif

end subroutine get_file_dims

end module mpp_io_utils_mod
