module top_wind_mod

use shr_kind_mod,      only : r4 => shr_kind_r4
use shr_kind_mod,      only : r8 => shr_kind_r8
use mpp_io_mod,        only : get_unit_num
use linear_interp_mod, only : linear_interp
use grid_mod,          only : nx, ny, rpts

implicit none
private

logical :: first_time=.true.

real(KIND=r4), dimension(rpts) :: normrad
real(KIND=r4), dimension(rpts) :: ratio

logical, parameter :: add_storm_trans = .true.
real(KIND=r4)      :: flagreal = -999.0_r4

public :: norm_dist_ratio

contains

subroutine norm_dist_ratio (xgrid, ygrid, dist, ubot, vbot, utx, uty, rmax, utop, vtop)

   real(KIND=r4), dimension(nx),    intent(in)  :: xgrid(nx)
   real(KIND=r4), dimension(ny),    intent(in)  :: ygrid(ny)
   real(KIND=r4), dimension(nx,ny), intent(in)  :: dist(nx, ny)
   real(KIND=r4), dimension(nx,ny), intent(in)  :: ubot
   real(KIND=r4), dimension(nx,ny), intent(in)  :: vbot
   real(KIND=r4),                   intent(in)  :: utx 
   real(KIND=r4),                   intent(in)  :: uty 
   real(KIND=r4),                   intent(in)  :: rmax
   real(KIND=r4), dimension(nx,ny), intent(out) :: utop
   real(KIND=r4), dimension(nx,ny), intent(out) :: vtop

   real(KIND=r4) :: uvbot(nx,ny)
   real(KIND=r4) :: dmax(1), scale(1)
   real(KIND=r4) :: trfact
   real(KIND=r4) :: vmax

   integer :: i, j, imax, jmax

   call norm_dist_ratio_init

 ! Setup initial parameters
   imax = 0
   jmax = 0
   vmax = flagreal

!  Calculate uvbot and also track what the max value of uvbot is using vmax
   do j=1,ny
     do i=1,nx  
       uvbot(i,j) = SQRT(ubot(i,j) * ubot(i,j) + vbot(i,j) * vbot(i,j))
       if (uvbot(i,j) > vmax) then
         vmax = uvbot(i,j)
         imax = i
         jmax = j
       endif
     enddo
   enddo

   do j=1,ny
     do i=1,nx
       dmax(1) = dist(i,j)/dist(imax,jmax)

       if (dmax(1) >= normrad(1) .AND. dmax(1)<= normrad(rpts)) then

           call linear_interp ( flagreal, flagreal, ratio, scale, normrad, dmax, rpts, 1 )
       else
           scale(1) = 1.0
       endif
       utop(i,j) = ubot(i,j) * 1  !scale(1) gave overestimated value for wind speed. So, replaced it with 1. 
       vtop(i,j) = vbot(i,j) * 1
     enddo
   enddo

!  apply storm translation speed

   if (add_storm_trans) then

     do j=1,ny
       do i=1,nx
         if ( dist(i,j) <= rmax ) then
           trfact = 1.0
         else
           trfact = MIN(1.0, rmax/dist(i,j))
           trfact = rmax/dist(i,j)
         endif
         utop(i,j) = utop(i,j) + 0        !utx is zero for static ideal case. -5 m/s for ideal moving case
         vtop(i,j) = vtop(i,j) + 0        !uty is zero for both static and moving ideal case
       enddo
     enddo
   endif

end subroutine norm_dist_ratio

subroutine norm_dist_ratio_init 

   integer :: n, unit_topwind

   if (first_time) then
     call get_unit_num (unit_topwind)
     open (unit_topwind, file='INPUT/topwindrationorm.dat', status='OLD', form='FORMATTED')

     do n=1,rpts
       read (unit_topwind, '(f12.4, f12.4)') normrad(n), ratio(n)
     enddo

     close (unit_topwind)

     first_time = .false.
   endif

end subroutine norm_dist_ratio_init 

end module top_wind_mod

