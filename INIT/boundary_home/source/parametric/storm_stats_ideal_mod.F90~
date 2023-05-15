module storm_stats_ideal_mod

use shr_kind_mod,         only : r4 => shr_kind_r4
use shr_kind_mod,         only : r8 => shr_kind_r8
use mpp_setup_mod,        only : master, rank
use mpp_io_mod,           only : get_unit_num, error_handler, stdout
use bilinear_interp_mod,  only : bilinear_interp, pair
use storm_stats_utils_mod,only : storm_stats_init, date_to_day, horiz_interp

use grid_mod,             only : nx, ny, icen, jcen, rpts
use constants_mod,        only : fcor, rhoair, pie, rearth, deg2rad

!#define DEBUG 1

implicit none
private

#include 'storm_domain.h'           !Previous version, it's domain_wind_grid.h. Need to rename the domain_wind_grid.h as storm_domain.h in the parameteric grid directory. 

public :: storm_stats_ideal

integer, parameter :: flagint  = -999
real(KIND=r4)      :: flagreal = -999.0_r4

contains

subroutine storm_stats_ideal (xgrid, ygrid, time, ubot, vbot, xcen, ycen, xcor, ycor, &
                    utx, uty, rmax, storm_name, storm_track_file, storm_start_time, &
                                       storm_wind_scale_ratio_file, wind_scale_ratio)

   real(KIND=r4), dimension(nx),    intent(in)  :: xgrid
   real(KIND=r4), dimension(ny),    intent(in)  :: ygrid
   real(KIND=r4),                   intent(in)  :: time
   real(KIND=r4), dimension(nx,ny), intent(out) :: ubot
   real(KIND=r4), dimension(nx,ny), intent(out) :: vbot
   real(KIND=r4),                   intent(out) :: xcen, ycen
   real(KIND=r4), dimension(2),     intent(out) :: xcor, ycor
   real(KIND=r4),                   intent(out) :: utx, uty
   real(KIND=r4),                   intent(out) :: rmax
   character(LEN=32),               intent(in)  :: storm_name
   character(LEN=32),               intent(in)  :: storm_track_file
   character(LEN=32),               intent(in)  :: storm_start_time
   character(LEN=32),               intent(in)  :: storm_wind_scale_ratio_file 
   real(KIND=r4), dimension(rpts),  intent(in)  :: wind_scale_ratio

   integer            :: hour, lat, lon, mx, rmw
   integer(KIND=4)    :: startdate, date
   integer            :: ni, i, j, n1, nlines, nm1, np1, ii, ierr, indx
   integer            :: day, month, year, end_flag
   integer            :: garb(5), rd1(4), rd2(4) !garb(5) - extra columns in the track file we don't need !rd1 radius of 18m/s wind in the four quadrants !rd2 radius of 26m/s wind in the four quadrants 

   character(LEN=6 ) :: id
   character(LEN=3 ) :: name !when switching from running bob this will need changed
   character(LEN=1 ) :: ew, ns

   character(LEN=8 ) :: model_date
   
   real(KIND=r4), allocatable, dimension(:)   :: x, y, tm, pres, pres0, rmaxa, wspmax, radcls !rmaxa - radius of maximum wind !wspmax - maximum wind from the track file
   real(KIND=r4), allocatable, dimension(:,:) :: r18v, r26v !Demensions are: (the number of time steps from track file, the 4 quadrants)

   real(KIND=r4), dimension(im,jm) :: wndx !Wind after scalling and translation speed is applied 
   real(KIND=r4), dimension(im,jm) :: wndy !Wind after scalling and translation speed is applied 

   real(KIND=r4), dimension(im) :: lons
   real(KIND=r4), dimension(jm) :: lats
   real(KIND=r4), dimension(im) :: xdist
   real(KIND=r4), dimension(jm) :: ydist

   real(KIND=r4), dimension(im,jm) :: coriolis

   real(KIND=r4), dimension(5)  :: rref18v, rref26v, alphv !rref18v rref26v - Interpolated four quadrant winds 
   real(KIND=r4), dimension(14) :: rad, ws, radm, wsm, angl
   real(KIND=r4), dimension(3)  :: r_contin, wnd_contin

   real(KIND=r4) :: timev
   real(KIND=r4) :: rcls 
   real(KIND=r4) :: cmp, t1, t2, f0, f1, l0, l1, r, a7, b, c, e, delp, x0, y0 
   real(KIND=r4) :: a_contin, b_contin
   real(KIND=r4) :: pres1, pres2, wsmax
   real(KIND=r4) :: cl0, cnorth_e, x1, x2, y1, y2
   real(KIND=r4) :: prsmin, wndmax, ws18, ws26
   real(KIND=r4) :: wtan, wrad, rangl, alphaw
   real(KIND=r4) :: wx, wy, wm, wind_scale
   real(KIND=r4) :: rref18, rref26, wnd, utxa, utya !utxa and utya both have to do with translation speed
   real(KIND=r4) :: c1, c2, c3
   real(KIND=r4) :: deltax, deltax1, deltay, deltay1, dxdy, julday
   real(KIND=r4) :: xgridc, ygridc

   integer :: unit_track

   DATA rad  /0.,.4,.7,.8,.95,1.,1.35,2.7,4.05,5.4,6.75,8.1,10.8,13.5/
  
   !Original Code:
  ! DATA angl /0.,2.,4.,6.,7.,7.,14.,23.,24.,22.,21.,21.,21.,21./

   !For ideal case
   DATA angl /0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./
   
   type(pair) :: ig1, jg1

   rmax       = 40.e3
   e          = exp(1.0_r4)
   timev      = time
   wind_scale = 1.0_r4

   model_date=storm_start_time(1:4)//storm_start_time(6:7)//storm_start_time(9:10)//storm_start_time(12:13)

!   print *, storm_start_time(1:4)
!   print *, storm_start_time(6:7)
!   print *, storm_start_time(9:10)
!   print *, storm_start_time(12:13)
   read(model_date,'(i10)') startdate

   call get_unit_num (unit_track)

!  reading message file
   !a -> character
   !i -> integer
   !x -> skip

   open (unit_track, file=storm_track_file, status='old') 

  nlines = 1 !for ideal

   allocate (    x(nlines))
   allocate (    y(nlines))
   allocate (   tm(nlines))
   allocate ( pres(nlines))
   allocate ( pres0(nlines))
   allocate ( rmaxa(nlines))
   allocate (wspmax(nlines))
   allocate (radcls(nlines))
   allocate (  r18v(nlines,4))
   allocate (  r26v(nlines,4))


   do i=1,nlines
         read (unit_track,1000) id, name, date, hour, lat, ns, lon, ew, garb, mx, rmw, rd1, rd2
1000   format (a6, 2x, a3, 2x, i8, 1x, i4, 1x, i3, a1, 1x, 1x, i3, a1, 1x, i3, 1x, i3, 3i5,  &
               1x, i2, 1x, i3, 1x, i3, 1x, i3, 1x, i3, 1x, i3, 1x, i3, 1x, i3, 1x, i3, 1x, i3)

#ifdef DEBUG
     write(stdout,*) 'id     = ', id
     write(stdout,*) 'name   = ', name
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
#endif
          
     date=date*100+hour/100
     call date_to_day( year, julday, date)
     tm(i) = julday !Julian Day Number
     x(i)  = lon/10.0_r4                    ! eastern hemsiphere
     y(i)  = lat/10.0_r4                    ! northern hemisphere

     if (ew .eq. 'W') x(i)=-x(i)    ! western hemisphere
     if (ns .eq. 'S') y(i)=-y(i)    ! southern hemisphere

     pres(i)   = float(garb(3))
     pres0(i)  = float(garb(4))
     radcls(i) = float(garb(5))*1000
     wspmax(i) = float(mx)
     rmaxa(i)  = float(rmw)

    enddo 

   close(unit_track)
 
!  calculating starting day 

! Add the number of hours to the start date so it is correct format for
! date_to_year subroutine
   
   read(storm_start_time(12:13), '(i2)') hour
   startdate=startdate*100+hour/100 
   
   call date_to_day (year, julday, startdate)
   do i=1,nlines
     tm(i) = tm(i) - julday
   enddo 

   xcen= x(1)
   ycen= y(1)

   f0 = xcen
   l0 = ycen

   rcls = radcls(1)

   wsmax = wspmax(1)
   wndmax = wsmax
   wsmax  = wsmax*wind_scale

   rmax = rmaxa(1)
   rmax=rmax*1.e3

!  create wind grid (spherical)) (on parametric domain)
   do i= 1,im
     lons(i) = ( lonstart + (i-1)* resol )
   enddo

   do j= 1,jm
     lats(j) = ( latstart + (j-1)* resol )
   enddo

!  f0,l0 are the current (lonitude,lattitude) center coordinates of the storm
   x0  = f0
   y0  = l0
   f0  = f0*deg2rad
   l0  = l0*deg2rad
   cl0 = cnorth_e*deg2rad

   utx     = 0.0_r4 !for ideal
   uty     = 0.0_r4 !for ideal

!!!! Setting up to calculate the wind !!!!!!!

   do i=1,14
     radm(i) = rmax*rad(i)
   enddo 

   do i=1,im
     do j=1,jm
       coriolis(i,j) = fcor
     enddo
   enddo

  !!! Variables for outside the RMW !!!
  b = 1.3_r4
  a7 = rmax**b
  c = (b/e)**0.5
  delp = (wsmax/c)**2
  print *, 'b: ', b, 'a7: ', a7, 'c: ', c, 'delp: ',delp

  !!! Variables for inside the RMW !!!
  r_contin(1) = rmax-1.e3
  r_contin(2) = rmax
  r_contin(3) = rmax+1.e3

 do r = 1,3
     do i=1,im
     do j=1,jm
!  wnd_contin = sqrt( a7*b*delp*exp(-a7/r_contin(r)**b)/(rhoair*r_contin(r)**b) ) !Equation was partially written. Corrected it based on Holland(1980). rhoair was not used to match with Kun's code.
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


   !Start looping through all of the points in the grid
   do i=1,im
     do j=1,jm
      !Calculate the distance to the center of the hurricane and alphaw-(the angle of that distance?)
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
       alphaw = atan(abs(deltay/deltax))*sign(1.,dxdy) + &
                (1-sign(1.,dxdy))*pie/2+(1-sign(1.,deltay))*pie/2
       if(alphaw .ge. pie/4) then
         alphaw = alphaw-pie/4
       else
         alphaw = alphaw-pie/4+2*pie
       endif 

       call horiz_interp(14, 1, radm, angl, r, rangl) !This is the angle of deflection 

! calculate wind speed
       wnd  = 0.0_r4
       utxa = 0.0_r4
       utya = 0.0_r4

! inside the RMW ! kepert & wang 2001 (from Kun's matlab code (maybe there is something about this in his papers...))
        if(r.le.rmax) then
         wnd = c1*r + c2*r**2 + c3*r**3
        else
! outside !Holland Wind Profile
         wnd = sqrt(a7*b*delp*exp(-a7/r**b)/(r**b)+r**2* &        !removed rhoair to match with Kun's parametric code. 
               coriolis(i,j)**2/4.)-r*coriolis(i,j)/2.0_r4
       endif

         !translation speed
         utxa = utx !/2.0_r4
         utya = uty !/2.0_r4

         
!      scale wind ratio
!       do indx=1,rpts
!        if(r.ge.((indx-1)*stepx) .and. r.lt.(indx*stepx) ) then
!         wnd = wnd*wind_scale_ratio(indx)
!        endif
!       enddo

       rangl      = rangl*deg2rad
       wtan       = wnd*cos(rangl)
       wrad       =-wnd*sin(rangl)
       wx         = wrad*(deltax)/r-wtan*(deltay)/r
       wy         = wtan*(deltax)/r+wrad*(deltay)/r
       wm         = SQRT(wx*wx + wy*wy)

      ! method one
      ! wind_scale = MAX(1.-MAX(wm-15.,0.)/25.*0.2,0.75)
      ! wndx(i,j)= wx*wind_scale + MIN(1.0, (rmax/r))*utx/2
      ! wndy(i,j)= wy*wind_scale + MIN(1.0, (rmax/r))*uty/2

      ! method two
      ! wnd_scale = MAX(1.-MAX(wm-15.,0.)/25.*0.2,0.75)
      ! wndx(i,j)= wx*wind_scale + MIN(1.0, (rmax/r))*utx
      ! wndy(i,j)= wy*wind_scale + MIN(1.0, (rmax/r))*uty

      ! method three
       wind_scale = 1.0
       wndx(i,j)= wx + utxa
       wndy(i,j)= wy + utya
     enddo 
   enddo 

!  create wind grid (cartesian) (on parametric domain)
   do i=1,im
     xdist(i) = rearth * COS(ycen*deg2rad) * (lons(i) - xcen) * deg2rad  !
   enddo

   do j= 1,jm
     ydist(j) = rearth*(lats(j) - ycen) * deg2rad
   enddo

   xgridc = (xgrid(nx) - xgrid(1))*0.5_r4
   ygridc = (ygrid(ny) - ygrid(1))*0.5_r4

   xcor(1) = xcen - acos(xgridc/rearth)/(deg2rad*ycen)
   ycor(1) = ycen - ygridc/(rearth*deg2rad)
   xcor(2) = xcor(1) + 2 * (xcen - xcor(1))
   ycor(2) = ycor(1) + 2 * (ycen - ycor(1))

! interpolate winds onto parametric domain
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
     enddo
   enddo

   !Set center of the hurricane to zero
   ubot(icen,jcen) = 0.0_r4
   vbot(icen,jcen) = 0.0_r4

end subroutine storm_stats_ideal

end module storm_stats_ideal_mod
