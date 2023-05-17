module mpp_io_boundary_swath_mod

use shr_kind_mod,     only : r4 => shr_kind_r4

use mpp_setup_mod,    only : size, master, rank
use mpp_io_mod,       only : get_unit_num, error_handler
use mpp_io_mod,       only : stdout
use mpp_io_utils_mod, only : get_file_dims

   implicit none

   private

   include 'netcdf.inc'

   public :: write_boundary_swath

   character(LEN=48)  ::  read_swath_file ! file name associated with read
   character(LEN=48)  :: write_swath_file ! file name associated with write

   integer, parameter :: ndims = 2 ! number of dimensions in netCDF file

   integer :: start(ndims), count(ndims)

   character(LEN=32) :: var_units  
   character(LEN=32) :: lons_units = 'degrees'  
   character(LEN=32) :: lats_units = 'degrees'  

   character(LEN=32), parameter :: units = 'units'

!  create netCDF files, variables and dimensions and get back an ID for each
   integer, dimension(ndims) :: dimids
   integer :: ncid, lons_dimid, lats_dimid

   integer :: lons_varid
   integer :: lats_varid
   integer :: uvbot_varid
   integer :: rainfall_varid

!  loop indexes, and error handling
   integer :: i, j, k, retval, rec, id

contains

subroutine write_boundary_swath ( uvbot, rainfall, lons, lats, nx, ny )

   real(KIND=r4), dimension(nx,ny), intent(inout) :: uvbot
   real(KIND=r4), dimension(nx,ny), intent(inout) :: rainfall
   real(KIND=r4), dimension(nx)   , intent(in) :: lons
   real(KIND=r4), dimension(ny)   , intent(in) :: lats
   integer,                         intent(in) :: nx
   integer,                         intent(in) :: ny

   namelist /swath_write/ write_swath_file

!  loop indexes, and error handling
   integer :: i, j, k, retval, rec, id, io, io_unit

!  set swath file name
   call get_unit_num ( io_unit )

   open (io_unit,FILE='input.nml',STATUS='UNKNOWN')
   read (io_unit,NML=swath_write)
   close(io_unit,STATUS='KEEP')

   if ( rank == master ) then
     write(stdout, 1000) write_swath_file
1000 format( 2x, 'name of swath file written in: write_boundary_swath = ', a32 )
   endif

!  create netCDF file, overwrite if already exists
   call write_boundary_swath_state(1)

!  define the dimensions for each axis and get back ID for each
   retval = nf_def_dim(ncid, 'lons', nx, lons_dimid  )
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_swath/nf_def_dim/lons_dimid' )

   retval = nf_def_dim(ncid, 'lats', ny, lats_dimid  )
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_swath/nf_def_dim/lats_dimid' )

!  define coordinate variables
   retval = nf_def_var(ncid, 'lons', NF_REAL, 1, lons_dimid, lons_varid)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_swath/nf_def_var/lons_varid' )

   retval = nf_def_var(ncid, 'lats', NF_REAL, 1, lats_dimid, lats_varid)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_swath/nf_def_var/lats_varid' )

!  assign units attributes to coordinate variable data
   retval = nf_put_att_text(ncid, lons_varid, UNITS, len(TRIM(lons_units)), TRIM(lons_units))
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_swath/nf_put_att_text/lons_units' )

   retval = nf_put_att_text(ncid, lats_varid, UNITS, len(TRIM(lats_units)), TRIM(lats_units))
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_swath/nf_put_att_text/lats_units' )

!  the dimids array is used to pass the ID of the dimensions of the variables. note that in fortran
!  arrays are stored in column-major
   dimids(1) = lons_dimid
   dimids(2) = lats_dimid

!  define the attributes for var
   retval = nf_def_var(ncid, 'uvbot', NF_REAL, ndims, dimids, uvbot_varid )
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_swath/nf_def_var/uvbot_varid' )

   retval = nf_def_var(ncid, 'rainfall', NF_REAL, ndims, dimids, rainfall_varid )
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_swath/nf_def_var/rainfall_varid' )

!  define the units for var
   var_units = 'degrees'
   retval = nf_put_att_text(ncid, lons_varid, UNITS, len(TRIM(var_units)), var_units)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_swath/nf_put_att_text/var_units' )
  
   var_units = 'degrees'
   retval = nf_put_att_text(ncid, lats_varid, UNITS, len(TRIM(var_units)), var_units)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_swath/nf_put_att_text/var_units' )
 
   var_units = 'm/s' 
   retval = nf_put_att_text(ncid, uvbot_varid, UNITS, len(TRIM(var_units)), var_units)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_swath/nf_put_att_text/var_units' )
  
   var_units = 'm' 
   retval = nf_put_att_text(ncid, rainfall_varid, UNITS, len(TRIM(var_units)), var_units)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_swath/nf_put_att_text/var_units' )
  
!  end define mode; this tells netCDF we are done defining metadata
   retval = nf_enddef(ncid)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_swath/nf_enddef' )

!  put axes data
   retval = nf_put_var_REAL(ncid, lons_varid, lons)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_swath/nf_put_var_REAL/lons_axis' )

   retval = nf_put_var_REAL(ncid, lats_varid, lats)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_swath/nf_put_var_REAL/lats_axis' )

!  put var data
   count(1) = nx
   count(2) = ny
   start(1) = 1
   start(2) = 1

   retval = nf_put_vara_REAL(ncid, uvbot_varid, start, count, uvbot)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_swath/nf_put_var_REAL/uvbot' )

   retval = nf_put_vara_REAL(ncid, rainfall_varid, start, count, rainfall)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_swath/nf_put_var_REAL/rainfall')

   call write_boundary_swath_state(3)

end subroutine write_boundary_swath

subroutine write_boundary_swath_state(case)
   integer, intent(in) :: case
   integer :: retval

   select CASE(case)
     CASE(1)
!    create netCDF file, overwrite if already exists
     retval = nf_create(write_swath_file, nf_clobber, ncid)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_swath_state/nf_create' )

     CASE(2)
!    open netCDF file
     retval = nf_open(write_swath_file, nf_write, ncid)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_swath_state/nf_open' )

     CASE(3)
!    close netCDF file
     retval = nf_close(ncid)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_swath_state/nf_close' )

    end select

end subroutine write_boundary_swath_state

end module mpp_io_boundary_swath_mod
