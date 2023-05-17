module land_roughness_mod

use shr_kind_mod,        only : r4 => shr_kind_r4
use shr_kind_mod,        only : r8 => shr_kind_r8
use mpp_setup_mod,       only : master, rank
use mpp_io_mod,          only : stdout, error_handler
use mpp_io_utils_mod,    only : get_file_dims
use mpp_io_mod,          only : get_unit_num

use bilinear_interp_mod, only : bilinear_interp, pair
use mpp_read_mod,        only : mpp_read_axes, mpp_read_data
use grid_mod,            only : nx, ny
use constants_mod,       only : rearth, deg2rad

implicit none

private

!#define DEBUG 1

include 'netcdf.inc'

logical :: first_time_land_rough=.true.

real(KIND=r4) :: const_land_rough=0.03_r4

integer :: im, jm

integer, parameter :: flagint  = -9999
real(KIND=r4)      :: flagreal = -9999_r4

real(KIND=r4)      :: land_rough_fill = 0.03_r4

real(KIND=r4), allocatable,  dimension(:,:) :: land_rough_data

real(KIND=r4), allocatable, dimension(:) :: lons
real(KIND=r4), allocatable, dimension(:) :: lats
real(KIND=r4), allocatable, dimension(:) :: xdist
real(KIND=r4), allocatable, dimension(:) :: ydist

type(pair) :: ig1, jg1

public :: land_roughness

contains

subroutine land_roughness (land_rough, land_mask, xgrid, ygrid, xcen, ycen, &
                                                                 storm_name )

   real(KIND=r4), dimension(nx,ny), intent(out) :: land_rough
   real(KIND=r4), dimension(nx,ny), intent(in)  :: land_mask
   real(KIND=r4), dimension(nx),    intent(in)  :: xgrid
   real(KIND=r4), dimension(ny),    intent(in)  :: ygrid
   real(KIND=r4),                   intent(in)  :: xcen
   real(KIND=r4),                   intent(in)  :: ycen
   character(LEN=32),               intent(in)  :: storm_name

   namelist /roughness_land/ const_land_rough

   integer :: i, j, num_fill, count
   integer :: io_unit, ierr

   character(LEN=32) :: land_rough_file = 'INPUT/land_rough_storm_domain.nc'

   if ( TRIM(storm_name) == 'ideal') then
     do j=1,ny
       do i=1,nx
         land_rough(i,j) = flagreal ! over water
       enddo
     enddo
     return
   endif

   if ( first_time_land_rough ) then

     call get_unit_num (io_unit)
     open (io_unit,FILE='input.nml',STATUS='UNKNOWN')
     read (io_unit,NML=roughness_land)
     close(io_unit,STATUS='KEEP')

     if ( rank == master ) then
       write(stdout, 1000) const_land_rough
     endif
1000 format( 2x, 'land roughness, read in: land_roughness =', e13.5)

   endif

!  determine if water only
   count = 0
   do j=1,ny
     do i=1,nx
       if ( land_mask(i,j) == 0.0 ) then
         count = count + 1
       endif
     enddo
   enddo 
   if ( count == nx*ny ) then
     do j=1,ny
       do i=1,nx
         land_rough(i,j) = flagreal ! over water
       enddo
     enddo
     first_time_land_rough = .false.
     return
   endif

   if ( const_land_rough /= 0.0_r4 ) then
     do j=1,ny
       do i=1,nx
         if ( land_mask(i,j) == 0.0 ) then
           land_rough(i,j) = flagreal         ! over water
         else
           land_rough(i,j) = const_land_rough ! over land
         endif
       enddo
     enddo
     first_time_land_rough = .false.
     return
   endif
         
   if (first_time_land_rough) then

     call get_file_dims ( land_rough_file, im, 'lon', jm, 'lat' )

     allocate( land_rough_data(im,jm))
     allocate( lons(im))
     allocate( lats(jm))
     allocate( xdist(im))
     allocate( ydist(jm))

     call mpp_read_axes ( lons, lats, land_rough_file, im, jm ) 
     call mpp_read_data ( land_rough_data, 'land_rough', land_rough_file, im, jm )

     do j=1,jm
       do i=1,im
         if ( land_rough_data(i,j) == 0 ) then
           land_rough_data(i,j) = flagreal
         endif
       enddo
     enddo
     first_time_land_rough = .false.
   endif

!  create roughness grid (cartesian)
   do i=1,im
     xdist(i) = rearth * COS(ycen*deg2rad) * (lons(i) - xcen) * deg2rad 
   enddo

   do j=1,jm
     ydist(j) = rearth*(lats(j) - ycen) * deg2rad
   enddo

! interpolate land roughness onto parametric grid 

   do j=1,ny
     do i=1,nx

       ig1%lo = flagint
       ig1%up = flagint
       jg1%lo = flagint
       jg1%up = flagint

       call bilinear_interp (flagreal, flagreal, land_rough_data, land_rough(i,j), xdist, ydist, &
                                              xgrid(i), ygrid(j), im, jm, 1, 1, ig1, jg1, flagint)

       if ( land_mask(i,j) == flagreal ) then
         land_rough(i,j) = flagreal ! over water
       endif

     enddo
   enddo

!  fill in land points with no defined roughness
   num_fill = 0

   do j=1,ny
     do i=1,nx
       if( land_rough(i,j) == flagreal .AND. land_mask(i,j) == 1.0_r4 ) then
         land_rough(i,j) = land_rough_fill
         num_fill        = num_fill + 1
       endif
     enddo
   enddo

#ifdef DEBUG
   if ( num_fill /= 0 ) then
     if ( rank == master ) then
       write(stdout, 1000) nx*ny, num_fill 
     endif
   endif

1000 format( 2x, 'number of points in domain = ', i8, /, &
             2x, 'number of land fill points = ', i8 / )
#endif

end subroutine land_roughness

end module land_roughness_mod
