module storm_stats_exp_mod

use shr_kind_mod,          only : r4 => shr_kind_r4
use shr_kind_mod,          only : r8 => shr_kind_r8
use mpp_setup_mod,         only : master, rank
use mpp_io_mod,            only : get_unit_num, error_handler, stdout
use bilinear_interp_mod,   only : bilinear_interp, pair
use storm_stats_ideal_mod, only : storm_stats_ideal
use storm_stats_named_mod, only : storm_stats_named
use storm_stats_utils_mod, only : storm_stats_init, date_to_day, horiz_interp, &
                                  interp_quad, expwind_expo, expwind_poly
use grid_mod,              only : nx, ny, rpts

implicit none
private

public :: storm_stats_exp  

contains

subroutine storm_stats_exp (xgrid, ygrid, time, ubot, vbot, pout, xcen, ycen, xcor, ycor, &
                                     utx, uty, rmax, wsmax, storm_name, storm_start_time)

   real(KIND=r4), dimension(nx),    intent(in)    :: xgrid
   real(KIND=r4), dimension(ny),    intent(in)    :: ygrid
   real(KIND=r4),                   intent(in)    :: time
   real(KIND=r4), dimension(nx,ny), intent(out)   :: ubot
   real(KIND=r4), dimension(nx,ny), intent(out)   :: vbot
   real(KIND=r4), dimension(nx,ny), intent(out)   :: pout 
   real(KIND=r4),                   intent(out)   :: xcen, ycen
   real(KIND=r4), dimension(2),     intent(out)   :: xcor, ycor
   real(KIND=r4),                   intent(out)   :: utx, uty
   real(KIND=r4),                   intent(out)   :: rmax
   real(KIND=r4),                   intent(out)   :: wsmax
   character(LEN=32),               intent(out)   :: storm_name
   character(LEN=32),               intent(out)   :: storm_start_time

   character(LEN=32) :: storm_track_file
   character(LEN=32) :: storm_wind_scale_ratio_file
   real(KIND=r4), dimension(rpts) :: wind_scale_ratio

   call storm_stats_init (storm_name, storm_track_file, storm_start_time, &
                             storm_wind_scale_ratio_file, wind_scale_ratio)

   if ( storm_name == 'ideal' ) then
     call storm_stats_ideal (xgrid, ygrid, time, ubot, vbot, xcen, ycen, xcor, ycor, &
                   utx, uty, rmax, storm_name, storm_track_file, storm_start_time, &
                                      storm_wind_scale_ratio_file, wind_scale_ratio)
   else
     call storm_stats_named (xgrid, ygrid, time, ubot, vbot, pout, xcen, ycen, xcor, ycor, &
           utx, uty, rmax, wsmax,  storm_name, storm_track_file, storm_start_time, &
                                            storm_wind_scale_ratio_file, wind_scale_ratio)
   endif

end subroutine storm_stats_exp

end module storm_stats_exp_mod
