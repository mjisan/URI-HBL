module mpp_io_boundary_model_temporal_mod

use shr_kind_mod,          only : r8 => shr_kind_r8
use shr_kind_mod,          only : r4 => shr_kind_r4
use mpp_setup_mod,         only : size, master, rank
use mpp_io_mod,            only : get_unit_num, error_handler
use mpp_io_mod,            only : stdout, read_diag_table

   implicit none
   private

   include 'netcdf.inc'

   public :: write_boundary_model_temporal

   logical :: first_time = .true.

   character(LEN=16), allocatable, dimension(:) :: diag_table_list
   integer :: num_diag_table

   character(LEN=64)  :: file_name ! file name associated with diags

   character(LEN=32) :: var_units  ! units of variables

   character(LEN=32), parameter :: units   = 'units'

   character(LEN=32), parameter :: rec_name  = 't' ! name of time axis

!  when we create netCDF files, variables and dimensions, we get back an ID for each

   integer :: ncid, rec_dimid
   
   integer :: wsmax_varid
   integer :: x_center_varid, y_center_varid
   integer :: x_corner_west_varid
   integer :: x_corner_east_varid
   integer :: y_corner_south_varid
   integer :: y_corner_north_varid
  
contains
  
subroutine write_boundary_model_temporal ( wsmax, x_center, y_center, &
                                             x_corner, y_corner, recs )

   real(KIND=r8),               intent(in) :: wsmax
   real(KIND=r8),               intent(in) :: x_center, y_center
   real(KIND=r8), dimension(:), intent(in) :: x_corner, y_corner
   integer,                     intent(in) :: recs
  
   real(KIND=r4) :: wsmax_r4
   real(KIND=r4) :: x_center_r4
   real(KIND=r4) :: y_center_r4
   real(KIND=r4) :: x_corner_r4
   real(KIND=r4) :: y_corner_r4

!  loop indexes, and error handling
   integer :: i, j, k, retval, id, io_unit

   character(LEN=64) :: diag_file='OUTPUT/diagnostics_temporal.nc'

   if ( rank == master ) then
     if ( first_time == .true. ) then

       write(stdout, 1000) diag_file
1000 format( 2x, 'name of diag file, written in: write_boundary_model_temporal = ', a64 )

       file_name = diag_file

!      read diag_table
       call read_diag_table ( diag_table_list, num_diag_table )

!      create netCDF file, overwrite if already exists
       call write_boundary_model_state(1)

!      define the dimensions for each axis and get back ID for each

       retval = nf_def_dim(ncid, rec_name, NF_UNLIMITED, rec_dimid)
       if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_model_temporal/nf_def_dim/rec_dimid' )

!      define the attributes for var
       do i=1,num_diag_table

         if ( TRIM(diag_table_list(i)) == 'wsmax') then
           retval = nf_def_var(ncid, 'wsmax', NF_REAL, 1, rec_dimid, wsmax_varid)
           if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_model_temporal/nf_def_var/wsmax_varid' )
         endif

         if ( TRIM(diag_table_list(i)) == 'xy_center') then
           retval = nf_def_var(ncid, 'x_center', NF_REAL, 1, rec_dimid, x_center_varid)
           if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_model_temporal/nf_def_var/x_center_varid' )

           retval = nf_def_var(ncid, 'y_center', NF_REAL, 1, rec_dimid, y_center_varid)
           if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_model_temporal/nf_def_var/y_center_varid' )
         endif

         if ( TRIM(diag_table_list(i)) == 'xy_corner') then
           retval = nf_def_var(ncid, 'x_corner_west', NF_REAL, 1, rec_dimid, x_corner_west_varid)
           if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_model_temporal/nf_def_var/x_corner_varid' )

           retval = nf_def_var(ncid, 'x_corner_east', NF_REAL, 1, rec_dimid, x_corner_east_varid)
           if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_model_temporal/nf_def_var/x_corner_varid' )

           retval = nf_def_var(ncid, 'y_corner_south', NF_REAL, 1, rec_dimid, y_corner_south_varid)
           if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_model_temporal/nf_def_var/y_corner_varid' )

           retval = nf_def_var(ncid, 'y_corner_north', NF_REAL, 1, rec_dimid, y_corner_north_varid)
           if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_model_temporal/nf_def_var/y_corner_north_varid' )
         endif

         if ( TRIM(diag_table_list(i)) == 'wsmax') then
           var_units = ''
           retval = nf_put_att_text(ncid, wsmax_varid, UNITS, len(TRIM(var_units)), var_units)
           if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_model_temporal/nf_put_att_text/var_units' )
         endif

         if ( TRIM(diag_table_list(i)) == 'xy_center') then
           var_units = 'degrees'
           retval = nf_put_att_text(ncid, x_center_varid, UNITS, len(TRIM(var_units)), var_units)
           if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_model_temporal/nf_put_att_text/var_units' )

           var_units = 'degrees'
           retval = nf_put_att_text(ncid, y_center_varid, UNITS, len(TRIM(var_units)), var_units)
           if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_model_temporal/nf_put_att_text/var_units' )
         endif

         if ( TRIM(diag_table_list(i)) == 'xy_corner') then
           var_units = 'degrees'
           retval = nf_put_att_text(ncid, x_corner_west_varid, UNITS, len(TRIM(var_units)), var_units)
           if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_model_temporal/nf_put_att_text/var_units' )

           var_units = 'degrees'
           retval = nf_put_att_text(ncid, x_corner_east_varid, UNITS, len(TRIM(var_units)), var_units)
           if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_model_temporal/nf_put_att_text/var_units' )

           var_units = 'degrees'
           retval = nf_put_att_text(ncid, y_corner_south_varid, UNITS, len(TRIM(var_units)), var_units)
           if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_model_temporal/nf_put_att_text/var_units' )

           var_units = 'degrees'
           retval = nf_put_att_text(ncid, y_corner_north_varid, UNITS, len(TRIM(var_units)), var_units)
           if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_model_temporal/nf_put_att_text/var_units' )
         endif

       enddo

!      end define mode; this tells netCDF we are done defining metadata
       retval = nf_enddef(ncid)
       if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_model_temporal/nf_enddef' )

       first_time = .false.
     else
!    open netCDF file
       call write_boundary_model_state(2)
     endif

!    put var data 
     do i=1,num_diag_table

       if ( TRIM(diag_table_list(i)) == 'wsmax') then
         wsmax_r4 = wsmax
         retval = nf_put_var1_REAL(ncid, wsmax_varid, recs, wsmax_r4)
         if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_model_temporal/nf_put_var_REAL/wsmax' )
       endif

       if ( TRIM(diag_table_list(i)) == 'xy_center') then
         x_center_r4 = x_center
         retval = nf_put_var1_REAL(ncid, x_center_varid, recs, x_center_r4)
         if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_model_temporal/nf_put_var_REAL/x_center' )

         y_center_r4 = y_center
         retval = nf_put_var1_REAL(ncid, y_center_varid, recs, y_center_r4)
         if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_model_temporal/nf_put_var_REAL/y_center' )
       endif

       if ( TRIM(diag_table_list(i)) == 'xy_corner') then
         x_corner_r4 = x_corner(1)
         retval = nf_put_var1_REAL(ncid, x_corner_west_varid, recs, x_corner_r4)
         if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundarx_model/nf_put_var_REAL/x_corner_west' )

         x_corner_r4 = x_corner(2)
         retval = nf_put_var1_REAL(ncid, x_corner_east_varid, recs, x_corner_r4)
         if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundarx_model/nf_put_var_REAL/x_corner_east' )

         y_corner_r4 = y_corner(1)
         retval = nf_put_var1_REAL(ncid, y_corner_south_varid, recs, y_corner_r4)
         if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_model_temporal/nf_put_var_REAL/y_corner_south' )

         y_corner_r4 = y_corner(2)
         retval = nf_put_var1_REAL(ncid, y_corner_north_varid, recs, y_corner_r4)
         if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_model_temporal/nf_put_var_REAL/y_corner' )
       endif

     enddo
  
     call write_boundary_model_state(3)

    endif

end subroutine write_boundary_model_temporal

subroutine write_boundary_model_state(case)
   integer, intent(in) :: case
   integer :: retval

   select CASE(case)
     CASE(1)
!    create netCDF file, overwrite if already exists
     retval = nf_create(file_name, nf_clobber, ncid)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_model_state/nf_create' )

     CASE(2) 
!    open netCDF file
     retval = nf_open(file_name, nf_write, ncid)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_model_state/nf_open' )

     CASE(3) 
!    close netCDF file
     retval = nf_close(ncid)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_model_state/nf_close' )

    end select

end subroutine write_boundary_model_state

end module mpp_io_boundary_model_temporal_mod
