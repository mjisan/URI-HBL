module mpp_io_boundary_parametric_mod

use shr_kind_mod,     only : r8 => shr_kind_r8
use shr_kind_mod,     only : r4 => shr_kind_r4

use mpp_setup_mod,    only : size, master, rank
use mpp_grid_mod,     only : mpp_get_global_indices
use mpp_io_mod,       only : get_unit_num, add_rank_to_file_name, error_handler
use mpp_io_mod,       only : stdout
use mpp_io_utils_mod, only : check_file_dims

   implicit none
   private

   include 'netcdf.inc'

   public :: read_boundary_parametric, write_boundary_parametric

   logical :: first_time = .true.

   character(LEN=48)  :: param_file ! file name associated with param

   integer, parameter :: ndims = 3  ! number of dimensions in netCDF file

   integer :: start(ndims), count(ndims)

   character(LEN=32) :: var_units  ! units of variables

   character(LEN=32), parameter :: units = 'units'
   character(LEN=32), parameter :: times = 'time-units'
   character(LEN=32)            :: x_units = 'unit-less'
   character(LEN=32)            :: y_units = 'unit-less'
   character(LEN=48)            :: time_units
   real(KIND=r4)                :: time_level

   character(LEN=32), parameter :: x_name    = 'x' ! name of x axis
   character(LEN=32), parameter :: y_name    = 'y' ! name of y axis
   character(LEN=32), parameter :: rec_name  = 't' ! name of time axis

!  when we create netCDF files, variables and dimensions, we get back an ID for
!  each
   integer, dimension(ndims) :: dimids
   integer                   :: ncid, x_dimid, y_dimid, rec_dimid

   integer :: y_varid, x_varid, rec_varid

   integer :: x_center_varid
   integer :: y_center_varid
   integer :: x_corner_west_varid
   integer :: x_corner_east_varid
   integer :: y_corner_south_varid
   integer :: y_corner_north_varid
   integer :: utx_varid
   integer :: uty_varid
   integer :: rmax_varid
   integer :: wsmax_varid
   integer :: ubot_varid
   integer :: vbot_varid
   integer :: utop_varid 
   integer :: vtop_varid
   integer :: ug_varid
   integer :: vg_varid
   integer :: pgfx_varid
   integer :: pgfy_varid
   integer :: topog_varid
   integer :: land_mask_varid
   integer :: land_rough_varid
   integer :: pres_sea_varid

!  loop indexes, and error handling
   integer :: i, j, k, retval, rec, id

contains

subroutine read_boundary_parametric ( ubot, vbot, utop, vtop, pgfx, pgfy, &
                                            topog, land_mask, land_rough, &
                                  x_center, y_center, x_corner, y_corner, &
                                    nx, ny, nz, ngx, ngy, npx, npy, nrecs )

   real(KIND=r8), dimension(:,:,:),intent(out) :: ubot 
   real(KIND=r8), dimension(:,:,:),intent(out) :: vbot
   real(KIND=r8), dimension(:,:,:),intent(out) :: utop
   real(KIND=r8), dimension(:,:,:),intent(out) :: vtop
   real(KIND=r8), dimension(:,:,:),intent(out) :: pgfx
   real(KIND=r8), dimension(:,:,:),intent(out) :: pgfy
   real(KIND=r8), dimension(:,:,:),intent(out) :: topog
   real(KIND=r8), dimension(:,:,:),intent(out) :: land_mask
   real(KIND=r8), dimension(:,:,:),intent(out) :: land_rough
   real(KIND=r8), dimension(:),    intent(out) :: x_center
   real(KIND=r8), dimension(:),    intent(out) :: y_center
   real(KIND=r8), dimension(:,:),  intent(out) :: x_corner
   real(KIND=r8), dimension(:,:),  intent(out) :: y_corner

   integer,                        intent(in)  :: nx, ny, nz
   integer,                        intent(in)  :: ngx, ngy
   integer,                        intent(in)  :: npx, npy 
   integer,                        intent(in)  :: nrecs
  
!  local variables
   real(KIND=r4), dimension(nx,ny) :: var_read_r4
   real(KIND=r4)                   :: x_center_r4
   real(KIND=r4)                   :: y_center_r4
   real(KIND=r4)                   :: x_corner_r4
   real(KIND=r4)                   :: y_corner_r4

   integer :: ncid, rec, retval, id, io_unit
   integer :: xbeg_data, xend_data, ybeg_data, yend_data
   integer :: xbeg_comp, xend_comp, ybeg_comp, yend_comp
   integer :: xdim, ydim

   namelist /parametric/ param_file

!  set parametric file names
   call get_unit_num ( io_unit )

   open (io_unit,FILE='input.nml',STATUS='UNKNOWN')
   read (io_unit,NML=parametric)
   close(io_unit,STATUS='KEEP')

   if ( rank == master ) then
     write(stdout, 1000) param_file
1000 format( 2x, 'name of parametric file read in: read_boundary_parametric = ', a32 )
   endif

!  inquiry into information about the varables dimensions in the file
   call check_file_dims ( param_file, ngx, x_name, ngy, y_name, rec_size=nrecs, rec_name=rec_name )

!  open netCDF file
   retval = nf_open(param_file, nf_nowrite, ncid)

! get the varid of the var
   retval = nf_inq_varid(ncid, 'ubot', ubot_varid    )
   if (retval .NE. 0) call error_handler ( 'ERROR', 'read_boundary_parametric/nf_inq_varid/ubot_varid' )

   retval = nf_inq_varid(ncid, 'vbot', vbot_varid    )
   if (retval .NE. 0) call error_handler ( 'ERROR', 'read_boundary_parametric/nf_inq_varid/vbot_varid' )

   retval = nf_inq_varid(ncid, 'utop', utop_varid   )
   if (retval .NE. 0) call error_handler ( 'ERROR', 'read_boundary_parametric/nf_inq_varid/utop_varid' )

   retval = nf_inq_varid(ncid, 'vtop', vtop_varid    )
   if (retval .NE. 0) call error_handler ( 'ERROR', 'read_boundary_parametric/nf_inq_varid/vtop_varid' )

   retval = nf_inq_varid(ncid, 'pgfx', pgfx_varid  )
   if (retval .NE. 0) call error_handler ( 'ERROR', 'read_boundary_parametric/nf_inq_varid/pgfx_varid' )

   retval = nf_inq_varid(ncid, 'pgfy', pgfy_varid  )
   if (retval .NE. 0) call error_handler ( 'ERROR', 'read_boundary_parametric/nf_inq_varid/pgfy_varid' )

   retval = nf_inq_varid(ncid, 'topog', topog_varid )  
   if (retval .NE. 0) call error_handler ( 'ERROR', 'read_boundary_parametric/nf_inq_varid/topog_varid' )

   retval = nf_inq_varid(ncid, 'land_mask', land_mask_varid )  
   if (retval .NE. 0) call error_handler ( 'ERROR', 'read_boundary_parametric/nf_inq_varid/land_mask_varid' )

   retval = nf_inq_varid(ncid, 'land_rough', land_rough_varid)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'read_boundary_parametric/nf_inq_varid/land_rough_varid' )

   if ( rank == master ) then

     retval = nf_inq_varid(ncid, 'x_center', x_center_varid)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'read_boundary_parametric/nf_inq_varid/x_center_varid' )

     retval = nf_inq_varid(ncid, 'y_center', y_center_varid)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'read_boundary_parametric/nf_inq_varid/y_center_varid' )

     retval = nf_inq_varid(ncid, 'x_corner_west', x_corner_west_varid)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'read_boundary_parametric/nf_inq_varid/x_corner_west_varid' )

     retval = nf_inq_varid(ncid, 'x_corner_east', x_corner_east_varid)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'read_boundary_parametric/nf_inq_varid/x_corner_east_varid' )

     retval = nf_inq_varid(ncid, 'y_corner_south', y_corner_south_varid)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'read_boundary_parametric/nf_inq_varid/y_corner_south_varid' )

     retval = nf_inq_varid(ncid, 'y_corner_north', y_corner_north_varid)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'read_boundary_parametric/nf_inq_varid/y_corner_north_varid' )

   endif

   id = rank + 1
   call mpp_get_global_indices ( xbeg_data, xend_data, ybeg_data, yend_data, &
                                 xbeg_comp, xend_comp, ybeg_comp, yend_comp, &
                                              id, nx, ny, ngx, ngy, npx, npy ) 

   count(1) = nx
   count(2) = ny
   count(3) = 1 
   start(1) = xbeg_data
   start(2) = ybeg_data

!  get var data 
   do rec=1,nrecs
     start(3) = rec
     retval = nf_get_vara_REAL(ncid, ubot_varid, start, count, var_read_r4)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'read_boundary_parametric/nf_get_vara_REAL/ubot' )
     ubot(:,:,rec) = var_read_r4(:,:)

     retval = nf_get_vara_REAL(ncid, vbot_varid, start, count, var_read_r4) 
     if (retval .NE. 0) call error_handler ( 'ERROR', 'read_boundary_parametric/nf_get_vara_REAL/vbot' )
     vbot(:,:,rec) = var_read_r4(:,:)

     retval = nf_get_vara_REAL(ncid, utop_varid, start, count, var_read_r4)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'read_boundary_parametric/nf_get_vara_REAL/utop' )
     utop(:,:,rec) = var_read_r4(:,:)

     retval = nf_get_vara_REAL(ncid, vtop_varid, start, count, var_read_r4)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'read_boundary_parametric/nf_get_vara_REAL/vtop' )
     vtop(:,:,rec) = var_read_r4(:,:)

     retval = nf_get_vara_REAL(ncid, pgfx_varid, start, count, var_read_r4) 
     if (retval .NE. 0) call error_handler ( 'ERROR', 'read_boundary_parametric/nf_get_vara_REAL/pgfx' )
     pgfx(:,:,rec) = var_read_r4(:,:)

     retval = nf_get_vara_REAL(ncid, pgfy_varid, start, count, var_read_r4)  
     if (retval .NE. 0) call error_handler ( 'ERROR', 'read_boundary_parametric/nf_get_vara_REAL/pgfy' )
     pgfy(:,:,rec) = var_read_r4(:,:)

     retval = nf_get_vara_REAL(ncid, topog_varid, start, count, var_read_r4)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'read_boundary_parametric/nf_get_vara_REAL/topog' )
     topog(:,:,rec) = var_read_r4(:,:)

     retval = nf_get_vara_REAL(ncid, land_mask_varid, start, count, var_read_r4)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'read_boundary_parametric/nf_get_vara_REAL/land_mask' )
     land_mask(:,:,rec) = var_read_r4(:,:)

     retval = nf_get_vara_REAL(ncid, land_rough_varid, start, count, var_read_r4)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'read_boundary_parametric/nf_get_vara_REAL/land_rough' )
     land_rough(:,:,rec) = var_read_r4(:,:)

     if ( rank == master ) then

       retval = nf_get_var1_REAL(ncid, x_center_varid, rec, x_center_r4)
       if (retval .NE. 0) call error_handler ( 'ERROR', 'read_boundary_parametric/nf_get_vara_REAL/x_center' )
       x_center(rec) = x_center_r4

       retval = nf_get_var1_REAL(ncid, y_center_varid, rec, y_center_r4)
       if (retval .NE. 0) call error_handler ( 'ERROR', 'read_boundary_parametric/nf_get_vara_REAL/y_center' )
       y_center(rec) = y_center_r4

       retval = nf_get_var1_REAL(ncid, x_corner_west_varid, rec, x_corner_r4)
       if (retval .NE. 0) call error_handler ( 'ERROR', 'read_boundary_parametric/nf_get_vara_REAL/x_corner' )
       x_corner(rec,1) = x_corner_r4

       retval = nf_get_var1_REAL(ncid, x_corner_east_varid, rec, x_corner_r4)
       if (retval .NE. 0) call error_handler ( 'ERROR', 'read_boundary_parametric/nf_get_vara_REAL/x_corner' )
       x_corner(rec,2) = x_corner_r4

       retval = nf_get_var1_REAL(ncid, y_corner_south_varid, rec, y_corner_r4)
       if (retval .NE. 0) call error_handler ( 'ERROR', 'read_boundary_parametric/nf_get_vara_REAL/y_corner' )
       y_corner(rec,1) = y_corner_r4

       retval = nf_get_var1_REAL(ncid, y_corner_north_varid, rec, y_corner_r4)
       if (retval .NE. 0) call error_handler ( 'ERROR', 'read_boundary_parametric/nf_get_vara_REAL/y_corner' )
       y_corner(rec,2) = y_corner_r4

     endif


   enddo

!  close netcdf file
   retval = nf_close(ncid)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'read_boundary_parametric/nf_close' )

end subroutine read_boundary_parametric

subroutine write_boundary_parametric ( ubot, vbot, utop, vtop, pgfx, pgfy, &
                                   topog, land_mask, land_rough, pres_sea, &
            x_center, y_center, x_corner, y_corner, utx, uty, rmax, wsmax, &
                                            time, start_time, recs, nx, ny )

   real(KIND=r4) , dimension(:,:),intent(in) ::ubot 
   real(KIND=r4) , dimension(:,:),intent(in) :: vbot 
   real(KIND=r4) , dimension(:,:),intent(in) :: utop
   real(KIND=r4) , dimension(:,:),intent(in) :: vtop
   real(KIND=r4) , dimension(:,:),intent(in) :: pgfx
   real(KIND=r4) , dimension(:,:),intent(in) :: pgfy
   real(KIND=r4) , dimension(:,:),intent(in) :: topog
   real(KIND=r4) , dimension(:,:),intent(in) :: land_mask
   real(KIND=r4) , dimension(:,:),intent(in) :: land_rough
   real(KIND=r4) , dimension(:,:),intent(in) :: pres_sea
   real(kind=r4) ,                intent(in) :: x_center
   real(kind=r4) ,                intent(in) :: y_center
   real(kind=r4) , dimension(:),  intent(in) :: x_corner
   real(kind=r4) , dimension(:),  intent(in) :: y_corner
   real(kind=r4) ,                intent(in) :: utx
   real(kind=r4) ,                intent(in) :: uty
   real(kind=r4) ,                intent(in) :: rmax 
   real(kind=r4) ,                intent(in) :: wsmax 
   real(KIND=r4) ,                intent(in) :: time
   character(LEN=*),              intent(in) :: start_time
   integer,                       intent(in) :: recs
   integer,                       intent(in) :: nx
   integer,                       intent(in) :: ny

   namelist /parametric/ param_file

   real(KIND=r4), dimension(nx+1,ny+1) :: var_write

   real(kind=r8), dimension(nx+1) :: xaxis
   real(kind=r8), dimension(ny+1) :: yaxis

!  loop indexes, and error handling
   integer :: i, j, k, retval, rec, id, io, io_unit

   if ( first_time == .true. ) then

!  set parametric file name
   call get_unit_num ( io_unit )

   open (io_unit,FILE='input.nml',STATUS='UNKNOWN')
   read (io_unit,NML=parametric)
   close(io_unit,STATUS='KEEP')

   if ( rank == master ) then
     write(stdout, 1000) param_file
1000 format( 2x, 'name of parametric file written in: read_boundary_parametric = ', a32 )
   endif

!    define axis 
     do j = 1,ny+1
       yaxis(j) = j
     end do
     do i = 1,nx+1
       xaxis(i) = i
     end do

!    create netCDF file, overwrite if already exists
     call write_boundary_parametric_state(1)

!    define the dimensions for each axis and get back ID for each
     retval = nf_def_dim(ncid, y_name, ny+1, y_dimid  )
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_def_dim/y_dimid' )

     retval = nf_def_dim(ncid, x_name, nx+1, x_dimid  )
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_def_dim/x_dimid' )

     retval = nf_def_dim(ncid, rec_name, NF_UNLIMITED, rec_dimid)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_def_dim/rec_dimid' )

!    define the coordinate variables
     retval = nf_def_var(ncid, y_name, NF_DOUBLE, 1, y_dimid, y_varid)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_def_var/y_varid' )

     retval = nf_def_var(ncid, x_name, NF_DOUBLE, 1, x_dimid, x_varid)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_def_var/x_varid' )

!    assign units attributes to coordinate variable data

     retval = nf_put_att_text(ncid, x_varid, UNITS, len(TRIM(x_units)), TRIM(x_units))
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_hbl_param/nf_put_att_text/x_units' )

     retval = nf_put_att_text(ncid, y_varid, UNITS, len(TRIM(y_units)), TRIM(y_units))
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_hbl_param/nf_put_att_text/y_units' )

     time_units = 'start time: ' // TRIM(start_time)

     retval = nf_put_att_text(ncid, rec_varid, TIMES, len(TRIM(time_units)), TRIM(time_units))
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_att_text/time_units' )

!    the dimids array is used to pass the ID of the dimensions of the variables. note that in fortran
!    arrays are stored in column-major
     dimids(1) = x_dimid
     dimids(2) = y_dimid
     dimids(3) = rec_dimid

!    define the attributes for var

     retval = nf_def_var(ncid, 'time', NF_REAL, 1, rec_dimid, rec_varid )
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_hbl_param/nf_def_var/rec_varid' )

     retval = nf_def_var(ncid, 'x_center', NF_REAL, 1, rec_dimid, x_center_varid )
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_def_var/x_center_varid' )

     retval = nf_def_var(ncid, 'y_center', NF_REAL, 1, rec_dimid, y_center_varid )
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_def_var/y_center_varid' )

     retval = nf_def_var(ncid, 'x_corner_west', NF_REAL, 1, rec_dimid, x_corner_west_varid )
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_def_var/x_corner_west_varid' )

     retval = nf_def_var(ncid, 'x_corner_east', NF_REAL, 1, rec_dimid, x_corner_east_varid )
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_def_var/x_corner_east_varid' )

     retval = nf_def_var(ncid, 'y_corner_south', NF_REAL, 1, rec_dimid, y_corner_south_varid )
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_def_var/y_corner_south_varid' )

     retval = nf_def_var(ncid, 'y_corner_north', NF_REAL, 1, rec_dimid, y_corner_north_varid )
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_def_var/y_corner_varid' )

     retval = nf_def_var(ncid, 'utx', NF_REAL, 1, rec_dimid, utx_varid )
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_def_var/utx_varid' )

     retval = nf_def_var(ncid, 'uty', NF_REAL, 1, rec_dimid, uty_varid )
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_def_var/uty_varid' )

     retval = nf_def_var(ncid, 'rmax', NF_REAL, 1, rec_dimid, rmax_varid )
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_def_var/rmax_varid' )

     retval = nf_def_var(ncid, 'wsmax', NF_REAL, 1, rec_dimid, wsmax_varid )
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_def_var/wsmax_varid' )

     retval = nf_def_var(ncid, 'ubot', NF_REAL, ndims, dimids, ubot_varid )
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_def_var/ubot_varid' )

     retval = nf_def_var(ncid, 'vbot', NF_REAL, ndims, dimids, vbot_varid )
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_def_var/vbot_varid' )

     retval = nf_def_var(ncid, 'utop', NF_REAL, ndims, dimids, utop_varid )
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_def_var/utop_varid' )

     retval = nf_def_var(ncid, 'vtop', NF_REAL, ndims, dimids, vtop_varid )
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_def_var/vtop_varid' )

     retval = nf_def_var(ncid, 'pgfx', NF_REAL, ndims, dimids, pgfx_varid )
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_def_var/pgfx_varid' )

     retval = nf_def_var(ncid, 'pgfy', NF_REAL, ndims, dimids, pgfy_varid )
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_def_var/pgfy_varid' )

     retval = nf_def_var(ncid, 'topog' , NF_REAL, ndims, dimids, topog_varid )
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_def_var/topog_varid' )

     retval = nf_def_var(ncid, 'land_mask', NF_REAL, ndims, dimids, land_mask_varid )
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_def_var/land_mask_varid' )

     retval = nf_def_var(ncid, 'land_rough', NF_REAL, ndims, dimids, land_rough_varid)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_def_var/land_rough_varid' )

     retval = nf_def_var(ncid, 'pres_sea', NF_REAL, ndims, dimids, pres_sea_varid)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_def_var/pres_sea_varid' )

     var_units = 'hours'
     retval = nf_put_att_text(ncid, rec_varid, UNITS, len(TRIM(var_units)), TRIM(var_units))
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_hbl_param/nf_put_att_text/var_units' )

     var_units = 'degrees'
     retval = nf_put_att_text(ncid, x_center_varid, UNITS, len(TRIM(var_units)), var_units)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_att_text/var_units' )
  
     var_units = 'degrees'
     retval = nf_put_att_text(ncid, y_center_varid, UNITS, len(TRIM(var_units)), var_units)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_att_text/var_units' )
 
     var_units = 'degrees'
     retval = nf_put_att_text(ncid, x_corner_west_varid, UNITS, len(TRIM(var_units)), var_units)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_att_text/var_units' )
  
     var_units = 'degrees'
     retval = nf_put_att_text(ncid, x_corner_east_varid, UNITS, len(TRIM(var_units)), var_units)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_att_text/var_units' )
  
     var_units = 'degrees'
     retval = nf_put_att_text(ncid, y_corner_south_varid, UNITS, len(TRIM(var_units)), var_units)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_att_text/var_units' )
 
     var_units = 'degrees'
     retval = nf_put_att_text(ncid, y_corner_north_varid, UNITS, len(TRIM(var_units)), var_units)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_att_text/var_units' )
 
     var_units = '' 
     retval = nf_put_att_text(ncid, utx_varid, UNITS, len(TRIM(var_units)), var_units)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_att_text/var_units' )
  
     var_units = '' 
     retval = nf_put_att_text(ncid, uty_varid, UNITS, len(TRIM(var_units)), var_units)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_att_text/var_units' )
  
     var_units = '' 
     retval = nf_put_att_text(ncid, rmax_varid, UNITS, len(TRIM(var_units)), var_units)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_att_text/var_units' )
  
     var_units = '' 
     retval = nf_put_att_text(ncid, wsmax_varid, UNITS, len(TRIM(var_units)), var_units)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_att_text/var_units' )
  
     var_units = 'm/s' 
     retval = nf_put_att_text(ncid, ubot_varid, UNITS, len(TRIM(var_units)), var_units)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_att_text/var_units' )
  
     var_units = 'm/s' 
     retval = nf_put_att_text(ncid, vbot_varid, UNITS, len(TRIM(var_units)), var_units)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_att_text/var_units' )
  
     var_units = 'm/s' 
     retval = nf_put_att_text(ncid, utop_varid, UNITS, len(TRIM(var_units)), var_units)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_att_text/var_units' )
  
     var_units = 'm/s' 
     retval = nf_put_att_text(ncid, vtop_varid, UNITS, len(TRIM(var_units)), var_units)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_att_text/var_units' )
  
     retval = nf_put_att_text(ncid, pgfx_varid, UNITS, len(TRIM(var_units)), var_units)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_att_text/var_units' )

     var_units = '' 
     retval = nf_put_att_text(ncid, pgfy_varid, UNITS, len(TRIM(var_units)), var_units)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_att_text/var_units' )

     var_units = 'm' 
     retval = nf_put_att_text(ncid, topog_varid, UNITS, len(TRIM(var_units)), var_units)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_att_text/var_units' )

     var_units = '' 
     retval = nf_put_att_text(ncid, land_mask_varid, UNITS, len(TRIM(var_units)), var_units)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_att_text/var_units' )

     var_units = '' 
     retval = nf_put_att_text(ncid, land_rough_varid, UNITS, len(TRIM(var_units)), var_units)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_att_text/var_units' )

     var_units = 'kPa' 
     retval = nf_put_att_text(ncid, pres_sea_varid, UNITS, len(TRIM(var_units)), var_units)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_att_text/var_units' )

!    end define mode; this tells netCDF we are done defining metadata
     retval = nf_enddef(ncid)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_enddef' )

!    put axis data
     retval = nf_put_var_DOUBLE(ncid, y_varid, yaxis)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_var_DOUBLE/yaxis' )

     retval = nf_put_var_DOUBLE(ncid, x_varid, xaxis)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_var_DOUBLE/xaxis' )

     first_time = .false.
   else
!  open netCDF file
     call write_boundary_parametric_state(2)
   endif

   count(1) = nx+1
   count(2) = ny+1
   count(3) = 1 
   start(1) = 1 
   start(2) = 1 
   start(3) = recs

!  put var data 
   retval = nf_put_var1_REAL(ncid, rec_varid, recs, time)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_hbl_param/nf_put_var_REAL/time' )

   retval = nf_put_var1_REAL(ncid, x_center_varid, recs, x_center)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_var_REAL/x_center' )

   retval = nf_put_var1_REAL(ncid, y_center_varid, recs, y_center)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_var_REAL/y_center' )

   retval = nf_put_var1_REAL(ncid, x_corner_west_varid, recs, x_corner(1))
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_var_REAL/x_corner(1)' )

   retval = nf_put_var1_REAL(ncid, x_corner_east_varid, recs, x_corner(2))
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_var_REAL/x_corner(2)' )

   retval = nf_put_var1_REAL(ncid, y_corner_south_varid, recs, y_corner(1))
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_var_REAL/y_corner(1)' )

   retval = nf_put_var1_REAL(ncid, y_corner_north_varid, recs, y_corner(2))
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_var_REAL/y_corner(2)' )

   retval = nf_put_vara_REAL(ncid, utx_varid, recs, 1, utx)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_var_REAL/utx' )

   retval = nf_put_vara_REAL(ncid, uty_varid, recs, 1, uty)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_var_REAL/uty' )

   retval = nf_put_vara_REAL(ncid, rmax_varid, recs, 1, rmax)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_var_REAL/rmax' )

   retval = nf_put_vara_REAL(ncid, wsmax_varid, recs, 1, wsmax)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_var_REAL/wsmax' )

   call extra_column_row ( ubot, var_write, nx, ny )
   retval = nf_put_vara_REAL(ncid, ubot_varid, start, count, var_write)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_var_REAL/ubot' )

   call extra_column_row ( vbot, var_write, nx, ny )
   retval = nf_put_vara_REAL(ncid, vbot_varid, start, count, var_write)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_var_REAL/vbot' )

   call extra_column_row ( utop, var_write, nx, ny )
   retval = nf_put_vara_REAL(ncid, utop_varid, start, count, var_write)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_var_REAL/utop' )

   call extra_column_row ( vtop, var_write, nx, ny )
   retval = nf_put_vara_REAL(ncid, vtop_varid, start, count, var_write)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_var_REAL/vtop' )

   call extra_column_row ( pgfx, var_write, nx, ny )
   retval = nf_put_vara_REAL(ncid, pgfx_varid, start, count, var_write)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_var_REAL/pgfx' )

   call extra_column_row ( pgfy, var_write, nx, ny )
   retval = nf_put_vara_REAL(ncid, pgfy_varid, start, count, var_write)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_var_REAL/pgfy' )

   call extra_column_row ( topog, var_write, nx, ny )
   retval = nf_put_vara_REAL(ncid, topog_varid, start, count, var_write)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_var_REAL/topog')

   call extra_column_row ( land_mask, var_write, nx, ny )
   retval = nf_put_vara_REAL(ncid, land_mask_varid, start, count, var_write)
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_var_REAL/land_mask' )

   call extra_column_row ( land_rough, var_write, nx, ny )
   retval = nf_put_vara_REAL(ncid, land_rough_varid, start, count, var_write) 
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_var_REAL/land_rough' )

   call extra_column_row ( pres_sea, var_write, nx, ny )
   retval = nf_put_vara_REAL(ncid, pres_sea_varid, start, count, var_write) 
   if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric/nf_put_var_REAL/pres_sea' )

   call write_boundary_parametric_state(3)

end subroutine write_boundary_parametric

subroutine write_boundary_parametric_state(case)
   integer, intent(in) :: case
   integer :: retval

   select CASE(case)
     CASE(1)
!    create netCDF file, overwrite if already exists
     retval = nf_create(param_file, nf_clobber, ncid)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric_state/nf_create' )

     CASE(2)
!    open netCDF file
     retval = nf_open(param_file, nf_write, ncid)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric_state/nf_open' )

     CASE(3)
!    close netCDF file
     retval = nf_close(ncid)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_boundary_parametric_state/nf_close' )

    end select

end subroutine write_boundary_parametric_state

subroutine extra_column_row ( var_write_in, var_write_out, nx, ny )

   real(KIND=r4) , dimension(nx,ny),     intent(in)  :: var_write_in
   real(KIND=r4) , dimension(nx+1,ny+1), intent(out) :: var_write_out
   integer,                              intent(in)  :: nx
   integer,                              intent(in)  :: ny

   integer :: i, j

   do j=1,ny
     do i=1,nx
       var_write_out(i,j) = var_write_in(i,j)
     enddo
   enddo

   do i=1,nx
     var_write_out(i,ny+1) = var_write_in(i,ny)
   enddo

   do j=1,ny
     var_write_out(nx+1,j) =  var_write_in(nx,j)
   enddo

   var_write_out(nx+1,ny+1) = var_write_in(nx,ny)

end subroutine extra_column_row

end module mpp_io_boundary_parametric_mod
