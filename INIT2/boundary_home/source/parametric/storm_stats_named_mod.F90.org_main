module storm_stats_named_mod

use shr_kind_mod,         only : r4 => shr_kind_r4
use shr_kind_mod,         only : r8 => shr_kind_r8
use mpp_setup_mod,        only : master, rank
use mpp_io_mod,           only : get_unit_num, error_handler, stdout
use bilinear_interp_mod,  only : bilinear_interp, pair
use storm_stats_utils_mod,only : storm_stats_init, date_to_day, horiz_interp, &
                                 interp_quad, expwind_poly, expwind_expo
use grid_mod,             only : nx, ny, icen, jcen, rpts
use constants_mod,        only : pie, deg2rad, rearth, rhoair, fcor

#define DEBUG 1

implicit none
private

#include 'storm_domain.h'

public :: storm_stats_named

integer, parameter :: flagint  = -999
real(KIND=r4)      :: flagreal = -999.0_r4

real(KIND=r4) :: xcen_init, ycen_init

contains

subroutine storm_stats_named (xgrid, ygrid, time, ubot, vbot, pres_sea, xcen, ycen, xcor, ycor, &
                       utx, uty, rmax, wsmax, storm_name, storm_track_file, storm_start_time, & 
                                                 storm_wind_scale_ratio_file, wind_scale_ratio)

   real(KIND=r4), dimension(nx),    intent(in)  :: xgrid
   real(KIND=r4), dimension(ny),    intent(in)  :: ygrid
   real(KIND=r4),                   intent(in)  :: time
   real(KIND=r4), dimension(nx,ny), intent(out) :: ubot
   real(KIND=r4), dimension(nx,ny), intent(out) :: vbot
   real(KIND=r4), dimension(nx,ny), intent(out) :: pres_sea
   real(KIND=r4),                   intent(out) :: xcen, ycen
   real(KIND=r4), dimension(2),     intent(out) :: xcor, ycor
   real(KIND=r4),                   intent(out) :: utx, uty
   real(KIND=r4),                   intent(out) :: rmax
   real(KIND=r4),                   intent(out) :: wsmax 
   character(LEN=32),               intent(in)  :: storm_name
   character(LEN=32),               intent(in)  :: storm_track_file
   character(LEN=32),               intent(in)  :: storm_start_time
   character(LEN=32),               intent(in)  :: storm_wind_scale_ratio_file 

   real(KIND=r4), dimension(rpts),  intent(in)  :: wind_scale_ratio


   integer         :: hour, lat, lon, mx, rmw
   integer(KIND=4) :: startdate, date
   integer         :: ni, i, j, n1, nlines, nm1, np1, ii, ierr, indx
   integer         :: day, month, year, end_flag
   integer         :: garb(5) 
   integer         :: rd1(4)  
   integer         :: rd2(4)  
   real            :: rvg, r18n, r26n, delv18, delv26 
   
   character(LEN=7 ) :: id
   character(LEN=3 ) :: name_a3
   character(LEN=4 ) :: name_a4
   character(LEN=5 ) :: name_a5
   character(LEN=1 ) :: ew, ns

   character(LEN=8 ) :: model_date

   integer, allocatable, dimension(:) :: track_time
   
   ! rmaxa  - radius of maximum wind 
   ! wspmax - maximum wind from the track file
   real(KIND=r4), allocatable, dimension(:)   :: x, y, tm,tm_ut, pres, pres0, rmaxa, wspmax, radcls 

   real(KIND=r4), allocatable, dimension(:)   :: utx_vec, uty_vec
   real(KIND=r4), allocatable, dimension(:,:) :: r18v, r26v ! dimensions are: (the number of time steps from track file, the 4 quadrants)

   real(KIND=r4), dimension(im,jm) :: wndx ! wind after scaling and translation speed is applied 
   real(KIND=r4), dimension(im,jm) :: wndy ! wind after scaling and translation speed is applied 
   real(KIND=r4), dimension(im,jm) :: pres_sea_pnt

   real(KIND=r4), dimension(im) :: lons
   real(KIND=r4), dimension(jm) :: lats
   real(KIND=r4), dimension(im) :: xdist
   real(KIND=r4), dimension(jm) :: ydist

   real(KIND=r4), dimension(im,jm) :: coriolis

   real(KIND=r4), dimension(5)  :: rref18v, rref26v ! interpolated four quadrant winds 
   real(KIND=r4), dimension(5)  :: alphv 
   real(KIND=r4), dimension(14) :: rad, ws, radm, wsm, angl
   real(KIND=r4), dimension(3)  :: r_contin, wnd_contin

   real(KIND=r4) :: timev, delt_track
   real(KIND=r4) :: rcls 
   real(KIND=r4) :: cmp, t1, t2, f0, f1, l0, l1, r, a7, b, c, e, delp, x0, y0
   real(KIND=r4) :: a_contin, b_contin 
   real(KIND=r4) :: c1, c2, c3
   real(KIND=r4) :: pres1, pres2
   real(KIND=r4) :: cl0, cnorth_e, x1, x2, y1, y2
   real(KIND=r4) :: prsmin, wndmax, ws18, ws26, uts
   real(KIND=r4) :: wtan, wrad, rangl, alphaw, wnd_scale
   real(KIND=r4) :: wx, wy, wm, wind_scale
   real(KIND=r4) :: rref18, rref26, wnd, wnd1, wnd2, utxa, utya 
   real(KIND=r4) :: deltax, deltax1, deltay, deltay1, dxdy, julday
   real(KIND=r4) :: xgridc, ygridc

   integer :: unit_track

   DATA rad  /0.,.4,.7,.8,.95,1.,1.35,2.7,4.05,5.4,6.75,8.1,10.8,13.5/
   DATA angl /0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./

   type(pair) :: ig1, jg1

   e          = exp(1.0_r4)
   timev      = time
   wind_scale = 1.0_r4

   model_date=storm_start_time(1:4)//storm_start_time(6:7)//storm_start_time(9:10)//storm_start_time(12:13)

   read(model_date,'(i10)') startdate

   call get_unit_num (unit_track)

   open (unit_track, file=storm_track_file, status='old') 

   nlines = 0
   do i=1,1000
     read (unit_track, 1000, iostat=ierr)
     if ( ierr .eq. -1 ) cycle
     nlines = nlines + 1
   enddo

   allocate (         x(nlines))
   allocate (         y(nlines))
   allocate (        tm(nlines))
   allocate (   utx_vec(nlines))
   allocate (   uty_vec(nlines))
   allocate (      pres(nlines))
   allocate (     pres0(nlines))
   allocate (     rmaxa(nlines))
   allocate (    wspmax(nlines))
   allocate (    radcls(nlines))
   allocate (track_time(nlines))
   allocate (      r18v(nlines,4))
   allocate (      r26v(nlines,4))

   rewind(unit_track)

   do i=1,nlines
     select case (storm_name)

       case ('bob')
         read (unit_track,1000) id, name_a3, date, hour, lat, ns, lon, ew, garb, mx, rmw, rd1, rd2
1000   format (a6, 2x, a3, 2x, i8, 1x, i4, 1x, i3, a1, 1x, 1x, i3, a1, 1x, i3, 1x, i3, 3i5,  &
               1x, i2, 1x, i3, 1x, i3, 1x, i3, 1x, i3, 1x, i3, 1x, i3, 1x, i3, 1x, i3, 1x, i3)

       case ('rhody')
         read (unit_track,2000) id, name_a5, date, hour, lat, ns, lon, ew, garb, mx, rmw, rd1, rd2
2000   format (a6, 2x, a5, 2x, i8, 1x, i4, 1x, i3, a1, 1x, i3, a1, 1x, i3, 1x, i3, 3i5,  &
               1x, i2, 1x, i3, 1x, i3, 1x, i3, 1x, i3, 1x, i3, 1x, i3, 1x, i3, 1x, i3, 1x, i3)
     
       case ('irma') 
         read (unit_track,3000) id, name_a4, date, hour, lat, ns, lon, ew, garb, mx, rmw, rd1, rd2
3000   format (a7, 2x, a4, 2x, i8, 1x, i4, 1x, i3, a1, 1x, i3, a1, 1x, i3, 1x, i3, 3(i4, 1x), &
               1x, i3, 1x, i3, 1x, i4, 1x, i4, 1x, i4, 1x, i4, 1x, i4, 1x, i4, 1x, i4, 1x, i4)

       case ('storm')
         read (unit_track,4000) id, name_a5, date, hour, lat, ns, lon, ew, garb, mx, rmw, rd1, rd2
4000   format (a6, 2x, a5, 2x, i8, 1x, i4, 1x, i3, a1, 1x, i3, a1, 1x, i3, 1x, i3, 3i5,  &
               1x, i2, 1x, i3, 1x, i3, 1x, i3, 1x, i3, 1x, i3, 1x, i3, 1x, i3,1x, i3, 1x, i3)

       case ('carol')
         read (unit_track,5000) id, name_a5, date, hour, lat, ns, lon, ew, garb, mx, rmw, rd1, rd2
5000   format (a3, 1x, a5, 1x, i8, 1x, i4, 1x, i3, a1, 1x, i3, a1, 1x, i3, 1x, i3, 3i5,  &
               1x, i2, 1x, i3, 1x, i3, 1x, i3, 1x, i3, 1x, i3, 1x, i3, 1x, i3,1x, i3, 1x, i3)

       case ('michael')
         read (unit_track,6000) id, name_a4, date, hour, lat, ns, lon, ew, garb, mx, rmw, rd1, rd2
6000   format (a3, 2x, a7, 2x, i8, 1x, i4, 1x, i3, a1, 1x, i3, a1, 1x, i3, 1x, i3, 3(i4, 1x), &
               1x, i3, 1x, i3, 1x, i4, 1x, i4, 1x, i4, 1x, i4, 1x, i4, 1x, i4, 1x, i4, 1x, i4)

       case ('florence')
         read (unit_track,7000) id, name_a4, date, hour, lat, ns, lon, ew, garb, mx, rmw, rd1, rd2, rvg, r18n, r26n, delv18, delv26
7000   format (a3, 2x, a8, 2x, i8, 1x, i4, 1x, i3, a1, 1x, i3, a1, 1x, i3, 1x, i3, 3(i4, 1x), &
               1x, i3, 1x, i3, 1x, i4, 1x, i4, 1x, i4, 1x, i4, 1x, i4, 1x, i4, 1x, i4, 1x, i4, 1x, f4.2, 1x, f5.2, 1x, f5.2, 1x, f4.2, 1x, f3.2)

       case ('track')
         read (unit_track,8000) id, name_a4, date, hour, lat, ns, lon, ew, garb, mx, rmw, rd1, rd2
8000   format (a3, 1x, a5, 1x, i8, 1x, i4, 1x, i3, a1, 1x, i3, a1, 1x, i3, 1x, i3, 3(i4, 1x), &
               1x, i3, 1x, i3, 1x, i4, 1x, i4, 1x, i4, 1x, i4, 1x, i4, 1x, i4, 1x, i4, 1x, i4)

       case ('dorian')
         read (unit_track,9000) id, name_a4, date, hour, lat, ns, lon, ew, garb, mx, rmw, rd1, rd2
9000   format (a3, 1x, a6, 1x, i8, 1x, i4, 1x, i3, a1, 1x, i3, a1, 1x, i3, 1x, i3, 3(i4, 1x), &
               1x, i3, 1x, i3, 1x, i4, 1x, i4, 1x, i4, 1x, i4, 1x, i4, 1x, i4, 1x, i4, 1x, i4)

     
     end select

#ifdef DEBUG
     write(stdout,*) 'id     = ', id
     write(stdout,*) 'name   = ', storm_name
     write(stdout,*) 'date   = ', date
     write(stdout,*) 'hour   = ', hour
     write(stdout,*) 'lat    = ', lat
     write(stdout,*) 'ns     = ', ns 
     write(stdout,*) 'lon    = ', lon
     write(stdout,*) 'ew     = ', ew 
     write(stdout,*) 'garb   = ', garb
     write(stdout,*) 'mx     = ', mx
     write(stdout,*) 'rmw    = ', rmw
     write(stdout,*) 'rd1    = ', rd1
     write(stdout,*) 'rd2    = ', rd2
     write(stdout,*) 'rvg    = ', rvg
     write(stdout,*) 'r18n   = ', r18n
     write(stdout,*) 'r26n   = ', r26n
     write(stdout,*) 'delv18 = ', delv18
     write(stdout,*) 'delv26 = ', delv26
#endif

     track_time(i) = hour/100

     date=date*100+hour/100
     call date_to_day( year, julday, date)
     tm(i) = julday               ! julian Day Number
     x(i)  = lon/10.0_r4          ! eastern hemsiphere
     y(i)  = lat/10.0_r4          ! northern hemisphere

     if (ew .eq. 'W') x(i)=-x(i)  ! western hemisphere
     if (ns .eq. 'S') y(i)=-y(i)  ! southern hemisphere

     pres(i)   = float(garb(3))
     pres0(i)  = float(garb(4))
     radcls(i) = float(garb(5))*1000
     wspmax(i) = float(mx)
     rmaxa(i)  = float(rmw)
     
!    reformatting the four quadrant radius of 18 and 26 m/s wind.
!    end format will be a table with rows for each TCvitals timestep and columns for the four quadrants   

     do ni=1,4
       n1         = ni+(1-mod(ni,2))*sign(2,3-ni) !n1 values in order over the loop of ni: 1, 4, 3, 2
       r18v(i,ni) = rd1(n1)
       r26v(i,ni) = rd2(n1)

       !if the max wind is less then 26 m/s or the radius of 26 m/s wind is less then the radius of max wind (RMW) set value equal to flagreal

       if (wspmax(i) .le. 26.or.r26v(i,ni) .le. rmaxa(i)) r26v(i,ni) = flagreal
       if (wspmax(i) .le. 18.or.r18v(i,ni) .le. rmaxa(i)) r18v(i,ni) = flagreal

       ! if this condition were true it would not be a hurricane any more because the faster wind is closer to the center of the hurricane

       if (r26v(i,ni) .gt. r18v(i,ni))                    r26v(i,ni) = flagreal 
     enddo 
   enddo 

   close(unit_track)

!  calculate time interval of track time

   delt_track = 0.0_r4
   do i=1,nlines
     if ( track_time(i) == 0.0_r4 ) then
       delt_track = (track_time(i+1) - track_time(i))*3600.0_r4
     endif
     if ( delt_track .ne. 0.0_r4 ) exit
   enddo
 
!  calculating starting day 

!  add the number of hours to the start date so it is correct format for date_to_year subroutine
   
   read(storm_start_time(12:13), '(i2)') hour
   startdate=startdate*100+hour
 
   call date_to_day (year, julday, startdate)
   do i=1,nlines
     tm(i) = tm(i) - julday
   enddo 

!  interpolation of storm path to determine the current position
   cmp=tm(nlines)

   if (cmp .lt. time .or. time .lt. tm(1)) return

   call horiz_interp (nlines, 1, tm, x, timev, f0)
   call horiz_interp (nlines, 1, tm, y, timev, l0)
   xcen=f0
   ycen=l0

   call horiz_interp (nlines, 1, tm, pres,   timev, pres1)
   call horiz_interp (nlines, 1, tm, pres0,  timev, pres2)
   call horiz_interp (nlines, 1, tm, radcls, timev, rcls)
   delp=(pres2-pres1)*100.0_r4
        
   if ( delp .lt. 100.0_r4 ) delp = 100.0_r4
   call horiz_interp (nlines, 1, tm, wspmax, timev, wsmax)
   prsmin = pres1
   wndmax = wsmax
   wsmax  = wsmax
   ws18   = r18n+delv18
   ws26   = r26n+delv26

   call horiz_interp (nlines, 1, tm, rmaxa, timev, rmax)
   rmax=rmax*1.e3


!  interpolate the four quadrant winds beteen the last TCvitals step and the current time step 
   do ni=1,4
     call interp_quad(flagreal, nlines, 1, tm, r18v(1,ni), timev, rref18v(ni))
     if (rref18v(ni) .ne. flagreal) rref18v(ni) = rref18v(ni)*1.e3

     call interp_quad(flagreal, nlines, 1, tm, r26v(1,ni), timev, rref26v(ni))
     if (rref26v(ni) .ne. flagreal) rref26v(ni) = rref26v(ni)*1.e3
!   every 90 degrees around the unit circle, later in the code 2pie is added in the 5th location 
     alphv(ni) = (ni-1)*pie/2
   enddo 
   
   do ni=2,6
     n1  = mod(ni-1,4)+1
     nm1 = mod(ni-2,4)+1
     np1 = mod(ni,4)+1
      
     if (rref18v(n1) .eq. flagreal) then
       if (rref18v(nm1) .ne. flagreal) then
         if (rref18v(np1) .ne. flagreal) then
           rref18v(n1) = 0.5*(rref18v(nm1)+rref18v(np1))
         else
           rref18v(n1) = rref18v(nm1)
         endif 
       else
         if (rref18v(np1) .ne. flagreal) then
           rref18v(n1) = rref18v(np1)
         else
           rref18v(n1) = flagreal
         endif 
       endif 
     endif 
     
     if (rref26v(n1) .eq. flagreal) then
       if(rref26v(nm1) .ne. flagreal) then
         if(rref26v(np1) .ne. flagreal) then
           rref26v(n1) = 0.5*(rref26v(nm1)+rref26v(np1))
         else
           rref26v(n1) = rref26v(nm1)
         endif 
       else
         if(rref26v(np1) .ne. flagreal) then
           rref26v(n1) = rref26v(np1)
         else
           rref26v(n1) =flagreal 
         endif 
       endif 
     endif 
   enddo 


!  this is completes the circle that makes the hurricane - this is to make loops for interpolation simpler

   rref18v(5) = rref18v(1)
   rref26v(5) = rref26v(1)
     alphv(5) = alphv(4)+pie/2

!  create wind grid (spherical) on parametric domain

   do i= 1,im
     lons(i) = ( lonstart + (i-1)*resol )
   enddo

   do j= 1,jm
     lats(j) = ( latstart + (j-1)*resol )
   enddo

!  THIS IS THE OLD METHOD, JUST ADDING IT IN TO DO RUNS FOR CONFERANCE

!  f0,l0 are the current (lonitude,lattitude) center coordinates of the storm
   x0  = f0
   y0  = l0
   f0  = f0*deg2rad
   l0  = l0*deg2rad
   cl0 = cnorth_e*deg2rad

!  calculating utx and uty (storm speed)
   cmp = tm(nlines)
   do i=1,nlines
     if(abs(tm(i)-time) .le. cmp) then
       cmp = abs(tm(i)-time)
       ii  = i
     endif
   enddo
   if((tm(ii)-time).le.0 .and. ii .ne. nlines) then
     t1 = tm(ii)
     t2 = tm(ii+1)
     x1 = x(ii)
     x2 = x(ii+1)
     y1 = y(ii)
     y2 = y(ii+1)
   else
     t2 = tm(ii)
     t1 = tm(ii-1)
     x2 = x(ii)
     x1 = x(ii-1)
     y2 = y(ii)
     y1 = y(ii-1)
   endif

   if(ifplane .eq. 1) then
     deltax1 = rearth*cos(cl0)*(x2-x1)*deg2rad
   else
     deltax1 = rearth*cos(l0)*(x2-x1)*deg2rad
   endif
   deltay1 = rearth*(y2-y1)*2.0_r4*pie/360.0_r4
   utx     = deltax1/((t2-t1)*24.0_r4*3600.0_r4)
   uty     = deltay1/((t2-t1)*24.0_r4*3600.0_r4)


   uts     = sqrt(utx**2 + uty**2)


  b = 1.3_r4
  a7 = rmax**b
  c = (b/e)**0.5
  delp = (wsmax/c)**2

   do i=1,14
     radm(i) = rmax*rad(i)
   enddo 

   do i=1,im
     do j=1,jm
       coriolis(i,j) = fcor
     enddo
   enddo



! variables for inside the RMW 

   r_contin(1) = rmax-1.e3
   r_contin(2) = rmax
   r_contin(3) = rmax+1.e3

   do r=1,3
     do i=1,im
       do j=1,jm
!         wnd_contin = sqrt( a7*b*delp*exp(-a7/r_contin(r)**b)/(rhoair*r_contin(r)**b) )

          wnd_contin = sqrt(a7*b*delp*exp(-a7/r_contin(r)**b)/(r_contin(r)**b)+r_contin(r)**2* &
                                               coriolis(i,j)**2/4)-r_contin(r)*coriolis(i,j)/2.0
       enddo
     enddo
   enddo

   a_contin=wnd_contin(2)
   b_contin=(wnd_contin(3)+wnd_contin(1)-2*wnd_contin(2))/(1.e3)**2
   c3=(a_contin+b_contin/2*rmax**2)/rmax**3
   c2=b_contin/2-3*rmax*c3
   c1=-b_contin*rmax+3*c3*rmax**2

!  start looping through all of the points
   do i=1,im
     do j=1,jm
      ! calculate the distance to the center of the hurricane and alphaw-(the angle of that distance?)
       f1 = lons(i)*deg2rad
       l1 = lats(j)*deg2rad
       if(ifplane .eq. 1) then
         deltax = rearth*cos(cl0)*(f1-f0)
       else
         deltax=  rearth*cos(l0)*(f1-f0)
       endif 
       if(abs(deltax) .lt. 1.e-8) deltax=sign(1.e-8,deltax)
       deltay = rearth*(l1-l0)
       if(abs(deltay) .lt. 1.e-8) deltay=sign(1.e-8,deltay)
       dxdy   = deltax*deltay
       r      = sqrt(deltax**2+deltay**2)
       alphaw = atan(abs(deltay/deltax))*sign(1.,dxdy) +      &
                (1-sign(1.,dxdy))*pie/2+(1-sign(1.,deltay))*pie/2
       if(alphaw .ge. pie/4) then
         alphaw = alphaw-pie/4
       else
         alphaw = alphaw-pie/4+2*pie
       endif 

!      call horiz_interp(5, 1, alphv, rref18v, alphaw, rref18) !CN Commented out to make axisymmetric exponential wind
!      call horiz_interp(5, 1, alphv, rref26v, alphaw, rref26) !CN Commented out to make axisymmetric exponential wind
       call horiz_interp(14, 1, radm, angl, r, rangl) ! this is the angle of deflection 

!      take the average the 4 quadrants !!This bit is just for makeing axisymmetric exponential wind -CN

       rref18 =( rref18v(1) + rref18v(2) + rref18v(3) + rref18v(4) )/4!(didn't use the 5th value because it is the same as the first)
       rref26 =( rref26v(1) + rref26v(2) + rref26v(3) + rref26v(4) )/4


!      calclate wind speed
       wnd  = 0.0_r4
       utxa = 0.0_r4
       utya = 0.0_r4
     
!      calculate sea-level pressure 
       pres_sea_pnt(i,j)=prsmin*100+delp*exp(-a7/(r**b)) 

       if(rref18 .le. 0 .and. rref26 .le. 0) then

         ! for inside the RMW !kepert & wang 2001 (from Kun's matlab code (maybe there is something about this in his papers...))
         if(r.le.rmax) then
           wnd = c1*r + c2*r**2 + c3*r**3
         else
         ! for Outside !Holland Wind Profile
         !  wnd = sqrt(a7*b*delp*exp(-a7/r**b)/(rhoair*r**b)+r**2* &
         !             coriolis(i,j)**2/4.)-r*coriolis(i,j)/2.0_r4
           wnd = sqrt(a7*b*delp*exp(-a7/r**b)/(r**b)+r**2* &        !removed rhoair to match with Kun's parametric code.
                      coriolis(i,j)**2/4.)-r*coriolis(i,j)/2.0_r4

         endif

         ! translation speed

         utxa = utx !/2.0_r4
         utya = uty !/2.0_r4

      else
         
         ! exponential wind profile (used when there is four quadrent information)

         wnd1=expwind_poly(r, rmax, rref18, rref26, wsmax, ws18, ws26)
         wnd2=expwind_expo(r, rmax, rref18, rref26, wsmax, ws18, ws26)

         if (wnd2 .lt. wnd1) then
            wnd = wnd2
         else
            wnd = wnd1
         end if   
         ! translation speed
         
            
         utxa = 0.0_r4
         utya = 0.0_r4
       endif 
       ! elect to make wind in f-plane runs perfectly axisymmetric
       if(ifplane .eq. 1) then
         utxa = 0.0_r4
         utya = 0.0_r4
       endif 

       rangl      = rangl*deg2rad
       rangl      = 0.0
       wtan       = wnd*cos(rangl)
       wrad       =-wnd*sin(rangl)
       wx         = wrad*(deltax)/r-wtan*(deltay)/r
       wy         = wtan*(deltax)/r+wrad*(deltay)/r
       wm         = SQRT(wx*wx + wy*wy)

       wind_scale = 1.0
       wndx(i,j)= wx + utx
       wndy(i,j)= wy + uty

!       wndx(i,j)= wx + 0.0
!       wndy(i,j)= wy + 0.0
       
     enddo 
   enddo 

   ! create wind grid (cartesian) on parametric domain
   do i=1,im
     xdist(i) = rearth*COS(ycen*deg2rad)*(lons(i) - xcen)*deg2rad  
   enddo

   do j= 1,jm
     ydist(j) = rearth*(lats(j) - ycen)*deg2rad 
   enddo

   xgridc = (xgrid(nx) - xgrid(1))*0.5_r4
   ygridc = (ygrid(ny) - ygrid(1))*0.5_r4

   xcor(1) = xcen - acos(xgridc/rearth)/(deg2rad*ycen) 
   ycor(1) = ycen - ygridc/(rearth*deg2rad) 

   xcor(2) = xcor(1) + 2 * (xcen - xcor(1))
   ycor(2) = ycor(1) + 2 * (ycen - ycor(1))

!  interpolate winds onto parametric domain
   do j=1,ny
     do i= 1,nx
       ig1%lo = flagint
       ig1%up = flagint
       jg1%lo = flagint
       jg1%up = flagint
       call bilinear_interp(flagreal, 0.0, wndx, ubot(i,j), xdist, ydist, &
                      xgrid(i), ygrid(j), im, jm, 1, 1, ig1, jg1, flagint)

       call bilinear_interp(flagreal, 0.0, wndy, vbot(i,j), xdist, ydist, &
                      xgrid(i), ygrid(j), im, jm, 1, 1, ig1, jg1, flagint)

       call bilinear_interp(flagreal, 0.0, pres_sea_pnt, pres_sea(i,j), xdist, ydist, &
                      xgrid(i), ygrid(j), im, jm, 1, 1, ig1, jg1, flagint)
     enddo
   enddo
   ! set center of the hurricane to zero
   ubot(icen,jcen) = 0.0_r4
   vbot(icen,jcen) = 0.0_r4

end subroutine storm_stats_named

end module storm_stats_named_mod
