module storm_stats_utils_mod

use shr_kind_mod,        only : r4 => shr_kind_r4
use shr_kind_mod,        only : r8 => shr_kind_r8
use mpp_setup_mod,       only : master, rank
use mpp_io_mod,          only : get_unit_num, error_handler, stdout
use bilinear_interp_mod, only : bilinear_interp, pair
use grid_mod,            only : nx, ny, rpts

!#define DEBUG 1

implicit none
private

logical :: first_time = .true.

public :: storm_stats_init, date_to_day, horiz_interp, interp_quad, expwind_poly, expwind_expo

contains

subroutine storm_stats_init (storm_name, storm_track_file, storm_start_time, &
                                storm_wind_scale_ratio_file, wind_scale_ratio)

   character(LEN=32), intent(out) :: storm_name 
   character(LEN=32), intent(out) :: storm_track_file
   character(LEN=32), intent(out) :: storm_start_time
   character(LEN=32), intent(out) :: storm_wind_scale_ratio_file

   real(KIND=r4), dimension(rpts), intent(out) :: wind_scale_ratio

   integer :: io_unit, ierr
   integer :: unit_track, unit_scale
   integer :: i

   namelist /storm_stats/ storm_name, storm_track_file, storm_start_time, &
                          storm_wind_scale_ratio_file
   
   call get_unit_num (io_unit)

   storm_wind_scale_ratio_file = ''

   if ( first_time ) then
     open (io_unit,FILE='input.nml',STATUS='UNKNOWN')
     read (io_unit,NML=storm_stats)
     close(io_unit,STATUS='KEEP')

     if ( rank == master ) then
       write(stdout, 1000) storm_name, storm_track_file, storm_start_time, &
                           storm_wind_scale_ratio_file
1000 format( 2x, 'name of storm read in: storm_stats_init                       = ', a32, / &
             2x, 'name of storm track file read in: storm_stats_init            = ', a32, / &
             2x, 'storm start time read in: storm_stats_init                    = ', a32, / &
             2x, 'name of storm wind_scale_ratio file read in: storm_stats_init = ', a32)
     endif

     if (storm_wind_scale_ratio_file .EQ. '') then
       do i=1,rpts
         wind_scale_ratio(i) = 1.0_r4
       enddo
     else
       call get_unit_num (unit_scale)
       open (unit_scale, file=storm_wind_scale_ratio_file, status='old')
       do i=1,rpts
         read (unit_track, 2000, iostat=ierr) wind_scale_ratio(i)
       enddo
     endif

2000 format(f14.6)

     first_time = .false.
   endif

end subroutine storm_stats_init

subroutine date_to_day ( year, julday, date )

   integer ,        intent(out) :: year
   real(KIND=r4),  intent(out) :: julday
   integer (KIND=4), intent(in) :: date

   integer , dimension(12) :: dat2day
   integer , dimension(12) :: dat2dayl
   integer  :: month, hour
   integer  ni

   real(KIND=r8) :: date_fp, date_dif, month_fp, month_dif, hour_dec
   integer :: date_int, month_int
        

   data dat2day  /31,28,31,30,31,30,31,31,30,31,30,31/
   data dat2dayl /31,29,31,30,31,30,31,31,30,31,30,31/

   year   = int(date/1000000.)
   month_fp = date/1000000.0_r8
   month_int = int(date/1000000)
   month_dif = (month_fp - month_int)*100
   month = NINT(month_dif)

   julday = 0
   if(mod(year,4)  .eq. 0) then
     do ni=1,month-1
       julday = julday+dat2dayl(ni)
     enddo
   else
     do ni=1,month-1
       julday = julday+dat2day(ni)
     enddo
   endif
    date_fp = date/10000.0_r8
    date_int = int(date/10000)
    date_dif = (date_fp - date_int)*100
    date_dif = NINT(date_dif) 
    julday = julday + date_dif

   hour = date-int(date/100)*100 
   hour_dec = hour/24.
   julday = julday+hour_dec

end subroutine date_to_day
 
subroutine horiz_interp ( nv, nvi, xv, yv, xvi, yvi,dop )

!CN: nv (indx of last step) 
!CN: nvi (doesnt get used, I think it was supposed to just represent the index of the first step (1) )
!CN: xv (Vector: with all of the time or space increments)
!CN: yv (Vector: with all of the values for the variable over space or time)
!CN: xvi (Time or location at which we want to interpolate a value for)
!CN: yvi (end value for the variable that is being interpolated for)

   integer,                      intent(in) :: nv, nvi
   real(KIND=r4), dimension(nv), intent(in) :: xv
   real(KIND=r4), dimension(nv), intent(in) :: yv
   real(KIND=r4),                intent(in) :: xvi
   real(KIND=r4),               intent(out) :: yvi
   logical, optional :: dop
   real(KIND=r4), parameter :: flagreal=-999.0_r4
   real(KIND=r4) tmpv(nv)

   integer :: n, ivi
   integer :: iv, jv
   !CN: Check to see if the time or distance is outside of the time or distance Vector supplied:
   if((xvi .gt. xv(nv) .and. xvi .gt. xv(1)) .or.      &
      (xvi .lt. xv(1) .and. xvi .lt. xv(nv))) then
    !CN: If it is greater then the last value in the vector:
     if(xvi .gt. xv(nv) .and. xvi .gt. xv(1)) then
       if(xv(nv) .gt. xv(1)) then
         yvi = yv(nv)
       else
         yvi = yv(1)
       endif
    !CN: If it is less then the last value in the vector:
     else
       if(xv(nv) .gt. xv(1)) then
         yvi = yv(1)
       else
         yvi = yv(nv)
       endif
     endif
  !CN: If it is in the range of the vector supplied:  
   else
     do jv=1,nv-1
       tmpv(jv) = (xvi-xv(jv))*(xvi-xv(jv+1))
     enddo
     do jv=1,nv-1
       if(tmpv(jv) .le. 0.0_r4) then
         ivi = jv
       endif
     enddo
     !this is the interpolation
     yvi = (yv(ivi)*abs(xv(ivi+1)-xvi)+yv(ivi+1)*       &
               abs(xvi-xv(ivi)))/abs(xv(ivi+1)-xv(ivi))
   endif
end subroutine horiz_interp
 
subroutine interp_quad ( flagreal, nv, nvi, xv, yv, xvi, yvi )

!   this subroutine determines ni values yi at the points xi
!   interpolating between the n values y at the points x
!   values .eq. al to mask are ig .or. d

   integer,                      intent(in)  :: nv,nvi
   real(KIND=r4), dimension(nv), intent(in)  :: xv, yv
   real(KIND=r4),                intent(in)  :: flagreal, xvi
   real(KIND=r4),                intent(out) :: yvi
   real(KIND=r4) tmpv(nv) 
   integer ::  ivi, iv, jv

   do iv=1,nvi
     if((xvi .gt. xv(nv) .and. xvi .gt. xv(1)) .or.    &
        (xvi .lt. xv(1) .and. xvi .lt. xv(nv))) then
        if(xvi .gt. xv(nv) .and. xvi .gt. xv(1)) then
          if(xv(nv) .gt. xv(1)) then
            yvi = yv(nv)
          else
            yvi = yv(1)
          endif
        else
          if(xv(nv) .gt. xv(1)) then
            yvi = yv(1)
          else
            yvi = yv(nv)
          endif
        endif
      else
        do jv=1,nv-1
          tmpv(jv) = (xvi-xv(jv))*(xvi-xv(jv+1))
        enddo
        do jv=1,nv-1
          if(tmpv(jv) .le. 0.0_r4) ivi = jv
        enddo
        if(yv(ivi) .eq. flagreal .or. yv(ivi+1) .eq. flagreal) then
          yvi = flagreal
        else
          yvi = (yv(ivi)*abs(xv(ivi+1)-xvi)+yv(ivi+1)*  &
                abs(xvi-xv(ivi)))/abs(xv(ivi+1)-xv(ivi))
        endif
      endif
    enddo

end subroutine interp_quad




function expwind_poly (r, rmax, rref18, rref26, wsmax, ws18, ws26, n )
 
   real(KIND=r4), intent(in) ::  r, rmax, rref18, rref26, wsmax, ws18, ws26

   real(KIND=r4) :: b, expwind_poly, r1, ws, n1, n2

   real, intent(out) :: n

   
   r1 = 0.5*(rref18+rref26)
   ws = 0.5*(ws18+ws26)

   n1 = abs(log10(ws/wsmax))

   n2 = abs(log10(r1/rmax))

   n  = n1/n2


   if(rref18 .le. 0.0_r4) then
     r1 = rref26
     ws = ws26
   endif

   if(rref26 .le. 0.0_r4) then
     r1 = rref18
     ws = ws18
   endif

   if(r .ge. rmax) then

      expwind_poly = wsmax*((rmax/r)**n)

   else

      expwind_poly = r*wsmax/rmax

   endif

end function expwind_poly


function expwind_expo (r, rmax, rref18, rref26, wsmax, ws18, ws26 )

   real(KIND=r4), intent(in) ::  r, rmax, rref18, rref26, wsmax, ws18, ws26

   real(KIND=r4) :: b, expwind_expo, r1, ws

   r1 = 0.5*(rref18+rref26)
   ws = 0.5*(ws18+ws26)

   if(rref18 .le. 0.0_r4) then
     r1 = rref26
     ws = ws26
   endif

   if(rref26 .le. 0.0_r4) then
     r1 = rref18
     ws = ws18
   endif

   if(r .ge. rmax) then
     b      = (rmax-r1)/log(ws/wsmax)
     expwind_expo = wsmax*exp((rmax-r)/b)
   else
     expwind_expo = r*wsmax/rmax
   endif

end function expwind_expo


end module storm_stats_utils_mod
