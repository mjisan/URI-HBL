module topography_mod

use shr_kind_mod,        only : r4 => shr_kind_r4
use shr_kind_mod,        only : r8 => shr_kind_r8
use mpp_setup_mod,       only : master, rank
use mpp_io_mod,          only : stdout, error_handler, get_unit_num
use mpp_io_utils_mod,    only : get_file_dims

use bilinear_interp_mod, only : bilinear_interp, pair
use mpp_read_mod,        only : mpp_read_axes, mpp_read_data
use grid_mod,            only : nx, ny
use constants_mod,       only : rearth, deg2rad

implicit none

private

include 'netcdf.inc'

logical :: first_time_topog=.true.

integer :: im, jm

integer, parameter :: flagint  = -9999
real(KIND=r4)      :: flagreal = -9999_r4

real(KIND=r4), allocatable,  dimension(:,:) :: topog_data

real(KIND=r4), allocatable, dimension(:) :: lons
real(KIND=r4), allocatable, dimension(:) :: lats
real(KIND=r4), allocatable, dimension(:) :: xdist
real(KIND=r4), allocatable, dimension(:) :: ydist

real(KIND=r4) :: topog_null

type(pair) :: ig1, jg1

public :: topography

contains

subroutine topography (topog, land_mask, xgrid, ygrid, xcen, ycen, &
                                                        storm_name )

   real(KIND=r4), dimension(nx,ny), intent(out) :: topog
   real(KIND=r4), dimension(nx,ny), intent(out) :: land_mask
   real(KIND=r4), dimension(nx),    intent(in)  :: xgrid
   real(KIND=r4), dimension(ny),    intent(in)  :: ygrid
   real(KIND=r4),                   intent(in)  :: xcen
   real(KIND=r4),                   intent(in)  :: ycen
   character(LEN=32),               intent(in)  :: storm_name

   integer :: i, j
   integer :: io_unit

   character(LEN=32) :: topog_file = 'INPUT/topog_storm_domain.nc'
   namelist /topog_elevation/ topog_null

   if ( TRIM(storm_name) == 'ideal') then

     do j=1,ny
       do i=1,nx
         land_mask(i,j) = 0.0_r4
         topog(i,j)     = 0.0_r4
       enddo
     enddo
     return
   endif

   if (first_time_topog) then

     call get_unit_num (io_unit)
     open (io_unit,FILE='input.nml',STATUS='UNKNOWN')
     read (io_unit,NML=topog_elevation)
     close(io_unit,STATUS='KEEP')

     if ( rank == master ) then
       write(stdout, 1000) topog_null
     endif
1000 format( 2x, 'topography elevation, read in: topography =', e13.5)

   endif

   if ( topog_null == flagreal ) then
     do j=1,ny
       do i=1,nx
         land_mask(i,j) = 0.0_r4
         topog(i,j)     = 0.0_r4
       enddo
     enddo
     first_time_topog = .false.
     return
   endif

   if (first_time_topog) then
     call get_file_dims ( topog_file, im, 'lon', jm, 'lat' )

     allocate( topog_data(im,jm))
     allocate(       lons(im))
     allocate(       lats(jm))
     allocate(      xdist(im))
     allocate(      ydist(jm))

     call mpp_read_axes ( lons, lats, topog_file, im, jm ) 
     call mpp_read_data ( topog_data, 'topo', topog_file, im, jm )

     do j=1,jm
       do i=1,im
        if ( topog_data(i,j) <= 0 ) then
          topog_data(i,j) = flagreal
        endif
       enddo
     enddo

     first_time_topog = .false.

   endif

   if ( topog_null == flagreal ) then
     do j=1,ny
       do i=1,nx
         land_mask(i,j) = 0.0_r4
         topog(i,j)     = 0.0_r4
       enddo
     enddo
     return
   endif

!  create topog grid (cartesian)
    do i=1,im
      xdist(i) = rearth * COS(ycen*deg2rad) * (lons(i) - xcen) * deg2rad 
    enddo

    do j= 1,jm
      ydist(j) = rearth*(lats(j) - ycen) * deg2rad
    enddo

! interpolate topogography and land mask onto parametric grid 

   do j=1,ny
     do i=1,nx

       ig1%lo = flagint
       ig1%up = flagint
       jg1%lo = flagint
       jg1%up = flagint

       call bilinear_interp (flagreal, flagreal, topog_data, topog(i,j), xdist, ydist, &
                                    xgrid(i), ygrid(j), im, jm, 1, 1, ig1, jg1, flagint)

       if ( topog(i,j) <= 0 ) then
         land_mask(i,j) = 0.0_r4 ! over water
         topog(i,j)     = 0.0_r4
       else
         land_mask(i,j) = 1.0_r4 ! over land
       endif 

     enddo
   enddo

end subroutine topography

end module topography_mod
