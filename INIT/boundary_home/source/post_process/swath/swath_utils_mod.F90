module swath_utils_mod

use shr_kind_mod,     only : r4 => shr_kind_r4
use mpp_setup_mod,    only : size, master, rank
use mpp_io_mod,       only : stdout, error_handler, get_unit_num
use mpp_read_mod,     only : mpp_read_data
use mpp_io_utils_mod, only : get_file_dims

implicit none

#include 'storm_domain.h'

private

  public :: swath_des_dims, swath_src_dims
  public :: swath_des_grid, swath_src_grid

 character(LEN=48) :: read_swath_file

 character(LEN=32), parameter :: x_name    = 'x' ! name of x axis
 character(LEN=32), parameter :: y_name    = 'y' ! name of y axis
 character(LEN=32), parameter :: rec_name  = 't' ! name of time axis

contains

  subroutine swath_des_dims ( nxdes, nydes )

  integer, intent(out) :: nxdes
  integer, intent(out) :: nydes
   
  nxdes = im 
  nydes = jm 

  end subroutine swath_des_dims 

  subroutine swath_src_dims ( nxsrc, nysrc, nrecs )


  integer, intent(out) :: nxsrc
  integer, intent(out) :: nysrc
  integer, intent(out) :: nrecs
 
  namelist /swath_read/ read_swath_file

! loop indexes, and error handling
  integer :: i, j, k, retval, rec, id, io, io_unit

! set swath file name
  call get_unit_num ( io_unit )

  open (io_unit,FILE='input.nml',STATUS='UNKNOWN')
  read (io_unit,NML=swath_read)
  close(io_unit,STATUS='KEEP')

  if ( rank == master ) then
    write(stdout, 1000) read_swath_file
1000 format( 2x, 'name of swath file written in: read_boundary_swath = ', a32 )
  endif

  call get_file_dims ( read_swath_file, nxsrc, x_name, nysrc, y_name, &
                                     rec_dim=nrecs, rec_name=rec_name )
 
  end subroutine swath_src_dims 

  subroutine swath_des_grid ( uvbot, rainfall, lon_des, lat_des, nxdes, nydes )

  real(KIND=r4), dimension(nxdes,nydes), intent(out) :: uvbot
  real(KIND=r4), dimension(nxdes,nydes), intent(out) :: rainfall
  real(KIND=r4), dimension(nxdes),       intent(out) :: lon_des
  real(KIND=r4), dimension(nydes),       intent(out) :: lat_des
  integer,                               intent(in)  :: nxdes
  integer,                               intent(in)  :: nydes

  integer :: i, j

! initialize destination grid
  do i=1,nxdes
    lon_des(i) =  lonstart + (i-1)*resol 
  enddo

  do j=1,nydes
    lat_des(j) =  latstart + (j-1)*resol 
  enddo

! initialize destination fields
  do j=1,nydes
    do i=1,nxdes
         uvbot(i,j) = 0.0_r4
      rainfall(i,j) = 0.0_r4
    enddo
  enddo
    
  end subroutine swath_des_grid

  subroutine swath_src_grid ( uvbot_src, rainfall_src, lon_src, lat_src, nxsrc, nysrc, recs )

  real(KIND=r4), dimension(nxsrc,nysrc), intent(out) :: uvbot_src
  real(KIND=r4), dimension(nxsrc,nysrc), intent(out) :: rainfall_src
  real(KIND=r4), dimension(nxsrc),       intent(out) :: lon_src
  real(KIND=r4), dimension(nysrc),       intent(out) :: lat_src
  integer,                               intent(in)  :: nxsrc
  integer,                               intent(in)  :: nysrc
  integer,                               intent(in)  :: recs 

  real(KIND=r4) :: x_corner_west, x_corner_east, y_corner_south, y_corner_north
  integer :: i, j

  real(KIND=r4), dimension(nxsrc,nysrc) :: ubot_src
  real(KIND=r4), dimension(nxsrc,nysrc) :: vbot_src
  real(KIND=r4) :: deltx, delty

! read source grid
  
  call  mpp_read_data ( x_corner_west,  'x_corner_west',  read_swath_file, recs )
  call  mpp_read_data ( x_corner_east,  'x_corner_east',  read_swath_file, recs )
  call  mpp_read_data ( y_corner_south, 'y_corner_south', read_swath_file, recs )
  call  mpp_read_data ( y_corner_north, 'y_corner_north', read_swath_file, recs )

  deltx = ABS((x_corner_east  - x_corner_west ))/FLOAT(nxsrc-1)
  delty = ABS((y_corner_north - y_corner_south))/FLOAT(nysrc-1)

  do i=1,nxsrc
    lon_src(i) = x_corner_west + deltx*(i-1)
  enddo

  do j=1,nysrc
    lat_src(j) = y_corner_south + delty*(j-1)
  enddo

! read source fields

#ifdef NO_WIND_FIELD
  ubot_src = 0.0
  vbot_src = 0.0
#else
  call mpp_read_data (ubot_src,     'ubot',     read_swath_file, nxsrc, nysrc, recs )
  call mpp_read_data (vbot_src,     'vbot',     read_swath_file, nxsrc, nysrc, recs )
#endif

#ifdef NO_RAINFALL
  rainfall_src = 0.0
#else
  call mpp_read_data (rainfall_src, 'rainfall', read_swath_file, nxsrc, nysrc, recs )
#endif

  do j=1,nysrc
    do i=1,nxsrc
       uvbot_src(i,j) = SQRT(ubot_src(i,j)**2 + vbot_src(i,j)**2)
      rainfall_src(i,j) = 10.0
    enddo
  enddo

end subroutine swath_src_grid 

end module swath_utils_mod
