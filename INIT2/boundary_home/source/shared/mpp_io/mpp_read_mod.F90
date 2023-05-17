module mpp_read_mod

  use shr_kind_mod,  only : r8 => shr_kind_r8
  use shr_kind_mod,  only : r4 => shr_kind_r4
  use mpp_setup_mod, only : size, master, rank
  use mpp_setup_mod, only : mpp_init, mpp_end
  use mpp_io_mod,    only : stdout, error_handler
  use netcdf


  implicit none
  private

  include 'netcdf.inc'

  public :: mpp_read_axes, mpp_read_data

  interface mpp_read_data
    module procedure mpp_read_data_2d
    module procedure mpp_read_data_0d
  end interface

contains

subroutine mpp_read_axes ( lons, lats, data_file, im, jm )

   real(KIND=r4),    dimension(jm),   intent(out) :: lats
   real(KIND=r4),    dimension(im),   intent(out) :: lons
   character(LEN=*),                  intent(in)  :: data_file
   integer,                           intent(in)  :: im, jm
  
   integer :: lons_varid, lats_varid
   integer :: ncid, retval

!  open netCDF file
   retval = nf90_open ( data_file, nf90_nowrite, ncid )
   if (retval .NE. 0) call error_handler ( 'ERROR', 'mpp_read/nf90_open/data_file')

!  get the varids of the latitude and longitude coordinate variables
   retval = nf90_inq_varid (ncid, 'lat', lats_varid)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'mpp_read/nf90_inq_varid/lat_varid' )

   retval = nf90_inq_varid (ncid, 'lon', lons_varid)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'mpp_read/nf90_inq_varid/lon_varid' )

!  read the latitude and longitude data
   retval = nf90_get_var (ncid, lats_varid, lats)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'mpp_read/nf90_get_var_real/lats' )

   retval = nf90_get_var (ncid, lons_varid, lons)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'mpp_read/nf90_get_var_real/lons' )

!  close netcdf file
   retval = nf90_close (ncid)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'mpp_read/nf90_close' )

end subroutine mpp_read_axes

subroutine mpp_read_data_0d ( data_out, data_name, data_file, nrec )

   real(KIND=r4),    intent(out) :: data_out
   character(LEN=*), intent(in)  :: data_name
   character(LEN=*), intent(in)  :: data_file
   integer,          intent(in)  :: nrec
  
   integer :: data_varid
   integer :: ncid, retval
   integer :: start(1), count(1) 

!  open netCDF file
   retval = nf90_open ( data_file, nf90_nowrite, ncid )
   if (retval .NE. 0) call error_handler ( 'ERROR', 'mpp_read/nf90_open/data_file')

!  get the varid of data
   retval = nf90_inq_varid (ncid, data_name, data_varid)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'mpp_read/nf90_inq_varid/data_varid' )

   start(1) = nrec
   count(1) = 1

!  read data
   retval = nf_get_vara_REAL (ncid, data_varid, start, count, data_out)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'mpp_read/nf_get_var_REAL/data_out' )

!  close netcdf file
   retval = nf90_close (ncid)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'mpp_read/nf90_close' )

end subroutine mpp_read_data_0d

subroutine mpp_read_data_2d ( data_out, data_name, data_file, im, jm, nrec )

   real(KIND=r4),    dimension(im,jm),intent(out) :: data_out
   character(LEN=*),                  intent(in)  :: data_name
   character(LEN=*),                  intent(in)  :: data_file
   integer,                           intent(in)  :: im, jm
   integer, optional,                 intent(in)  :: nrec
  
   integer :: data_varid
   integer :: ncid, retval
   integer, allocatable, dimension(:) :: start
   integer, allocatable, dimension(:) :: count

!  open netCDF file
   retval = nf90_open ( data_file, nf90_nowrite, ncid )
   if (retval .NE. 0) call error_handler ( 'ERROR', 'mpp_read/nf90_open/data_file')

!  get the varid of data
   retval = nf90_inq_varid (ncid, data_name, data_varid)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'mpp_read/nf90_inq_varid/data_varid' )

   if ( present(nrec) ) then
     allocate( start(3) )
     allocate( count(3) ) 
     start(1) = 1 
     start(2) = 1 
     start(3) = nrec
     count(1) = im
     count(2) = jm
     count(3) = 1
   else
     allocate( start(2) )
     allocate( count(2) )
     start(1) = 1
     start(2) = 1
     count(1) = im
     count(2) = jm
   endif

!  read data
   retval = nf90_get_var (ncid, data_varid, data_out, start, count)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'mpp_read/nf90_get_varid/data_out' )

!  close netcdf file
   retval = nf90_close (ncid)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'mpp_read/nf90_close' )

end subroutine mpp_read_data_2d

end module mpp_read_mod
