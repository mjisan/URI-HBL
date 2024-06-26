module model_step_mod

use shr_kind_mod,                    only : r8 => shr_kind_r8
use mpp_setup_mod,                   only : size, master, rank
use mpp_update_domain_mod,           only : update_domain
use mpp_io_boundary_parametric_mod,  only : read_boundary_parametric
use mpp_io_boundary_model_mod,       only : write_boundary_model
use mpp_io_mod,                      only : get_unit_num, error_handler, stdout, to_lower_case
use mpp_utils_mod,                   only : mpp_max, mpp_min, mpp_bcast
use mpp_setup_mod,                   only : mpp_barrier

use constants_mod,                   only : fcor
use roughness_mod,                   only : roughness_2
use grid_mod,                        only : nx, ny, nz, ngx, ngy, npx, npy 
use grid_mod,                        only : stepx, stepy, stepz
use time_mod,                        only : delt, nrecs, spinup_nrecs, spinup_nsteps, spinup_recs_int
use time_mod,                        only : forecast_recs_int,forecast_nrecs, spinup_nrecs

implicit none
private

#include 'grid_params_r8.h'

public model_step 

logical           :: model_initial = .true.

character(LEN=32) :: storm_name, storm_name_lower
logical           :: bot_top = .false.

real(KIND=r8), dimension(nx,ny,forecast_nrecs+1) :: ubot
real(KIND=r8), dimension(nx,ny,forecast_nrecs+1) :: vbot
real(KIND=r8), dimension(nx,ny,forecast_nrecs+1) :: ug
real(KIND=r8), dimension(nx,ny,forecast_nrecs+1) :: vg
real(KIND=r8), dimension(nx,ny,forecast_nrecs+1) :: pgfx
real(KIND=r8), dimension(nx,ny,forecast_nrecs+1) :: pgfy
real(KIND=r8), dimension(nx,ny,forecast_nrecs+1) :: topog
real(KIND=r8), dimension(nx,ny,forecast_nrecs+1) :: mask
real(KIND=r8), dimension(nx,ny,forecast_nrecs+1) :: landzo

real(KIND=r8), dimension(nx,ny) :: rain_acc

real(KIND=r8), dimension(forecast_nrecs+1)       :: xcen    
real(KIND=r8), dimension(forecast_nrecs+1)       :: ycen    

real(KIND=r8), dimension(forecast_nrecs+1,2)     :: xcor    
real(KIND=r8), dimension(forecast_nrecs+1,2)     :: ycor    

real(KIND=r8) :: x_center, y_center
real(KIND=r8), dimension(2) :: x_corner
real(KIND=r8), dimension(2) :: y_corner

real(KIND=r8), dimension(nrecs) :: mask_time

real(KIND=r8) :: pmix_len  = 8.0e+01
real(KIND=r8) :: turh      = 5.0e+04

real(KIND=r8) :: check_bound_time_beg  = -1.0_r8
real(KIND=r8) :: check_bound_time_end  = -1.0_r8
real(KIND=r8) :: flagcheck = 0.0_r8

real(KIND=r8) :: flagreal = -9999_r8

contains

subroutine model_step ( model_time) 

   real(KIND=r8), intent(in) :: model_time

   real(KIND=r8), dimension(nx,ny,nz) :: um
   real(KIND=r8), dimension(nx,ny,nz) :: vm
   real(KIND=r8), dimension(nx,ny,nz) :: wm
   real(KIND=r8), dimension(nx,ny,nz) :: tur
   real(KIND=r8), dimension(nx,ny,nz) :: umvm
 
   real(KIND=r8), dimension(nx,ny,2) :: tur_sl
 
   real(KIND=r8), dimension(nx,ny) :: pgfx1
   real(KIND=r8), dimension(nx,ny) :: pgfy1
   real(KIND=r8), dimension(nx,ny) :: mask1
   real(KIND=r8), dimension(nx,ny) :: ug1
   real(KIND=r8), dimension(nx,ny) :: vg1
   real(KIND=r8), dimension(nx,ny) :: cd
   real(KIND=r8), dimension(nx,ny) :: fric_x
   real(KIND=r8), dimension(nx,ny) :: fric_y
   real(KIND=r8), dimension(nx,ny) :: znot2d
   real(KIND=r8), dimension(nx,ny) :: lznot
 
   real(KIND=r8), dimension(nx) :: htur_u
   real(KIND=r8), dimension(nx) :: htur_v
   real(KIND=r8), dimension(nx) :: vtur_u
   real(KIND=r8), dimension(nx) :: vtur_v
   real(KIND=r8), dimension(nx) :: adv_u
   real(KIND=r8), dimension(nx) :: adv_v
   real(KIND=r8), dimension(nx) :: difu
   real(KIND=r8), dimension(nx) :: difv
   real(KIND=r8), dimension(nx) :: difw  
 
 
   real(KIND=r8) :: vref, znot, fric, ustar, uvbot, zkap, pmix, pmixsq
   real(KIND=r8) :: deter, conv
   real(KIND=r8) :: wsmax
 
   integer :: i, j, k, recs, rloc

!  if model is being in itialized then read in the file
   if ( model_initial ) then
     call model_init ( um, vm, wm, pgfx1, pgfy1, mask1 )
     model_initial = .false.
   endif

!  update upper and lower boundary information from input file which gets
!  updated every 15 minutes

   if ( TRIM(storm_name) /= 'ideal' ) then

     if(model_time > mask_time(spinup_nrecs)) then 

       do recs=spinup_nrecs,(nrecs-1)
         if (model_time >= mask_time(recs) .AND. model_time < mask_time(recs + 1)) then

           do j=1,ny
             do i=1,nx
               if (bot_top) then
                  um(i,j,nz) = ubot(i,j,recs-(spinup_nrecs-1))
                  vm(i,j,nz) = vbot(i,j,recs-(spinup_nrecs-1))
                else
                  um(i,j,nz) = ug(i,j,recs-(spinup_nrecs-1))
                  vm(i,j,nz) = vg(i,j,recs-(spinup_nrecs-1))
                endif
                  
                pgfx1(i,j)   =   pgfx(i,j,recs-(spinup_nrecs-1))
                pgfy1(i,j)   =   pgfy(i,j,recs-(spinup_nrecs-1))
                mask1(i,j)   =   mask(i,j,recs-(spinup_nrecs-1))
                lznot(i,j)   = landzo(i,j,recs-(spinup_nrecs-1))
                x_center     =   xcen(recs-(spinup_nrecs-1))
                y_center     =   ycen(recs-(spinup_nrecs-1))
                x_corner(1)  =   xcor(recs-(spinup_nrecs-1),1)
                x_corner(2)  =   xcor(recs-(spinup_nrecs-1),2)
                y_corner(1)  =   ycor(recs-(spinup_nrecs-1),1)
                y_corner(2)  =   ycor(recs-(spinup_nrecs-1),2)
             enddo
           enddo
         endif
       enddo

       if (model_time == mask_time(nrecs)) then
         do j=1,ny
           do i=1,nx
             if (bot_top) then
               um(i,j,nz) = ubot(i,j,nrecs-spinup_nrecs+1)
               vm(i,j,nz) = vbot(i,j,nrecs-spinup_nrecs+1)
             else
                um(i,j,nz) = ug(i,j,nrecs-spinup_nrecs+1)
                vm(i,j,nz) = vg(i,j,nrecs-spinup_nrecs+1)
             endif

             pgfx1(i,j)    =   pgfx(i,j,nrecs-spinup_nrecs+1)
             pgfy1(i,j)    =   pgfy(i,j,nrecs-spinup_nrecs+1)
             mask1(i,j)    =   mask(i,j,nrecs-spinup_nrecs+1)
             lznot(i,j)    = landzo(i,j,nrecs-spinup_nrecs+1)
             x_center      =   xcen(nrecs-spinup_nrecs+1)
             y_center      =   ycen(nrecs-spinup_nrecs+1)
             x_corner(1)   =   xcor(nrecs-spinup_nrecs+1,1)
             x_corner(2)   =   xcor(nrecs-spinup_nrecs+1,2)
             y_corner(1)   =   ycor(nrecs-spinup_nrecs+1,1)
             y_corner(2)   =   ycor(nrecs-spinup_nrecs+1,2)
           enddo
         enddo
       endif
     endif 

   endif 
 
!  compute znot
   do j=2,ny-1
     do i=2,nx-1
       vref = (um(i,j,2)**2+vm(i,j,2)**2)**0.5
       if (vref .le. 0.01_r8) vref = 0.01_r8

       ! check if over land or water
       if (mask1(i,j) /= 0.0_r8) then
         ! storm is over land
         znot = lznot(i,j)
       else
         ! storm is over water
         call roughness_2( vref, znot )
       endif
 
       cd(i,j)       = (0.4_r8/log(stepz/znot))**2
       fric          = cd(i,j)*vref**2
       ustar         = sqrt(fric)
       uvbot         = ustar/0.4*log(10./znot)
       fric_y(i,j)   = fric*vm(i,j,2)/vref
       fric_x(i,j)   = fric*um(i,j,2)/vref
       vm(i,j,1)     = uvbot*vm(i,j,2)/vref
       um(i,j,1)     = uvbot*um(i,j,2)/vref
       tur_sl(i,j,1) = 0.4_r8*10.0_r8*ustar
       tur_sl(i,j,2) = 0.4*stepz*ustar
       znot2d(i,j)   = znot
     enddo
   enddo

!  k=l^2*s*(1-ri)^0.5   -->  Km
   do k=2,nz-1
     zkap   = 0.4_r8*stepz*(k-1)
     pmix   = zkap/(1.0_r8+(zkap/pmix_len))
     pmixsq = pmix*pmix

     do j=2,ny-1
       do i=2,nx-1
         deter = ((vm(i,j,k+1)-vm(i,j,k-1))*stepz2inv)**2  &
                +((um(i,j,k+1)-um(i,j,k-1))*stepz2inv)**2
         if(deter < 1e-10) deter = 1e-10
         tur(i,j,k) = pmixsq*sqrt(deter)
       enddo
       do i=2,nx-1
         if (tur(i,j,k) >= 100.0_r8) tur(i,j,k) = 100.0_r8
       enddo
     enddo
   enddo

   do j=1,ny
     do i=1,nx
      tur(i,j,1)  = tur_sl(i,j,1)
      tur(i,j,nz) = tur(i,j,nz-1)
     enddo
   enddo

!  compute vm 

   do k=2,nz-1
     do j=2,ny-1
       do i=2,nx-1
         if ( um(i,j,k) >= 0.0_r8 ) then
           difu(i) = vm(i,j,k)-vm(i-1,j,k)
         else
           difu(i) = vm(i+1,j,k)-vm(i,j,k)
         endif
       enddo
       do i=2,nx-1
         if(vm(i,j,k) >= 0.0_r8) then
           difv(i) = vm(i,j,k)-vm(i,j-1,k)
         else
           difv(i) = vm(i,j+1,k)-vm(i,j,k)
         endif
       enddo

!      difu = 0.5*(vm(i+1,j,k)-vm(i-1,j,k))
!      difv = 0.5*(vm(i,j+1,k)-vm(i,j-1,k))

       do i=2,nx-1
         if(wm(i,j,k) >= 0.0_r8) then
           difw(i) = vm(i,j,k)-vm(i,j,k-1)
         else
           difw(i) = vm(i,j,k+1)-vm(i,j,k)
         endif
       enddo

       do i=2,nx-1
         if(k == 2) then
           vtur_v(i) = ((tur(i,j,k+1)+tur(i,j,k))*0.5_r8  &
           *(vm(i,j,k+1)-vm(i,j,k))*stepzinv-fric_y(i,j)) &
           *stepzinv
         else
           vtur_v(i) = ((tur(i,j,k+1)+tur(i,j,k))*0.5_r8 &
               *(vm(i,j,k+1)-vm(i,j,k))*stepzinv         &
               -(tur(i,j,k)+tur(i,j,k-1))*0.5_r8         &
               *(vm(i,j,k)-vm(i,j,k-1))*stepzinv)        &
               *stepzinv
         endif
       enddo

       do i=2,nx-1
           htur_v(i) = turh*((vm(i+1,j,k)-2.0_r8*vm(i,j,k)  &
                       +vm(i-1,j,k))*stepxsqinv             &
         +(vm(i,j+1,k)-2.0_r8*vm(i,j,k)+vm(i,j-1,k))*stepysqinv)

!        intergrate vm

         adv_v(i) = -um(i,j,k)*difu(i)*stepxinv       &
                    -vm(i,j,k)*difv(i)*stepyinv       &
                    -wm(i,j,k)*difw(i)*stepzinv

         vm(i,j,k) =  vm(i,j,k)+delt*( adv_v(i)        &
                      +pgfy1(i,j)-fcor*um(i,j,k)       &
                      +vtur_v(i)+htur_v(i) )
       enddo

     enddo
   enddo

   call update_domain ( vm, nx, ny, nz, npx, npy )

!  compute um 

   do k=2,nz-1
     do j=2,ny-1
       do i=2,nx-1
         if(um(i,j,k) >= 0.0_r8) then
           difu(i) = um(i,j,k)-um(i-1,j,k)
         else
           difu(i) = um(i+1,j,k)-um(i,j,k)
         endif
       enddo
       do i=2,nx-1
         if(vm(i,j,k) >= 0.0_r8) then
           difv(i) = um(i,j,k)-um(i,j-1,k)
         else
           difv(i) = um(i,j+1,k)-um(i,j,k)
         endif
       enddo

!      difu = 0.5*(um(i+1,j,k)-um(i-1,j,k))
!      difv = 0.5*(um(i,j+1,k)-um(i,j-1,k))

       do i=2,nx-1
         if(wm(i,j,k) >= 0.0_r8) then
           difw(i) = um(i,j,k)-um(i,j,k-1)
         else
           difw(i) = um(i,j,k+1)-um(i,j,k)
         endif
       enddo

       do i=2,nx-1
         if(k == 2) then
           vtur_u(i) = ((tur(i,j,k+1)+tur(i,j,k))*0.5_r8   &
           *(um(i,j,k+1)-um(i,j,k))*stepzinv -fric_x(i,j)) &
           *stepzinv
         else
           vtur_u(i) = ((tur(i,j,k+1)+tur(i,j,k))*0.5_r8   &
                        *(um(i,j,k+1)-um(i,j,k))*stepzinv  &
                         -(tur(i,j,k)+tur(i,j,k-1))*0.5_r8 &
                        *(um(i,j,k)-um(i,j,k-1))*stepzinv) &
                        *stepzinv
         endif
       enddo
       do i=2,nx-1
         htur_u(i) = turh*((um(i+1,j,k)-2.0_r8 &
          *um(i,j,k)+um(i-1,j,k))*stepxsqinv   &
          +(um(i,j+1,k)-2.*um(i,j,k)+um(i,j-1,k))*stepysqinv)

!        intergrate um

         adv_u(i) = -um(i,j,k)*difu(i)*stepxinv  &
                    -vm(i,j,k)*difv(i)*stepyinv  &
                    -wm(i,j,k)*difw(i)*stepzinv

         um(i,j,k) =  um(i,j,k)+delt*( adv_u(i)  &
                      +pgfx1(i,j)+fcor*vm(i,j,k) &
                      +vtur_u(i)+htur_u(i)  )
       enddo
     enddo
   enddo

   
   call update_domain ( um, nx, ny, nz, npx, npy )

!  compute wm

   do j=2,ny-1
     do i=2,nx-1
        conv = -(um(i+1,j,1)-um(i-1,j,1) &
                +vm(i,j+1,1)-vm(i,j-1,1) &
                +um(i+1,j,2)-um(i-1,j,2) &
                +vm(i,j+1,2)-vm(i,j-1,2))*0.5_r8*stepx2inv
        wm(i,j,2) = stepz*conv
     enddo
   enddo

   do k=3,nz
     do j=2,ny-1
       do i=2,nx-1
               conv = -(um(i+1,j,k-1)-um(i-1,j,k-1)  &
                       +vm(i,j+1,k-1)-vm(i,j-1,k-1))*stepx2inv
            wm(i,j,k) = wm(i,j,k-2)+2*stepz*conv
       enddo
     enddo
   enddo

!  compute accumulated rainfall


!  call rainfall ( vmax, rmax, rad, rain_acc )
 
!  check bounds of um and vm

   call model_check_bounds ( um, vm, model_time, flagcheck) 

!  write data every 15 minutes starting at zero and include last step

   do recs=1,nrecs
     if ( model_time == mask_time(recs) .OR. &
          flagcheck == 1.0 ) then

       if (recs < spinup_nrecs+1) then
         rloc = 1
       else
         rloc = recs-(spinup_nrecs-1)
       endif

!      compute maximum velocity in global domain

       do k=1,nz
         do j=1,ny
           do i=1,nx
             umvm(i,j,k) = SQRT(um(i,j,k)**2 + vm(i,j,k)**2)
            enddo
          enddo
        enddo

       call mpp_max ( umvm, wsmax )
          
       call update_domain (  wm, nx, ny, nz, npx, npy )
       call update_domain ( tur, nx, ny, nz, npx, npy )

       if ( rank == master ) then
         write(stdout,1000) model_time, rloc
1000 format( 2x, ' MODEL DIAGNOSTICS WRITTEN AT TIME = ', e13.6, /, &
             2x, ' LOCATION OF RECORD WRITTEN        = ', i4 )
       endif

       call write_boundary_model ( um, vm, pgfx(:,:,rloc), pgfy(:,:,rloc), wm, tur, &
            znot2d, mask1, rain_acc, wsmax, x_center, y_center, x_corner, y_corner, &
                                               recs, nx, ny, nz, ngx, ngy, npx, npy )

       if ( flagcheck == 1.0 ) then
         call error_handler ( 'ERROR', 'model_check' )
       endif

     endif
   enddo

end subroutine model_step 

subroutine model_check_bounds ( um, vm, model_time, flagcheck )

   real(KIND=r8), dimension(nx,ny,nz), intent(in) :: um
   real(KIND=r8), dimension(nx,ny,nz), intent(in) :: vm
   real(KIND=r8), intent(in)  :: model_time
   real(KIND=r8), intent(out) :: flagcheck

   integer :: i, j, k
   integer :: istop, jstop, kstop
   integer :: stopflag = 0, rloc
   real(KIND=r8), parameter :: uvbound = 199.0_r8

   real(KIND=r8) :: um_max, um_min, vm_max, vm_min
   real(KIND=r8) :: stopflags

   if ( model_time .GE. check_bound_time_beg .AND. &
        model_time .LE. check_bound_time_end  ) then
     call mpp_max ( um, um_max )
     call mpp_min ( um, um_min )
     call mpp_max ( vm, vm_max )
     call mpp_min ( vm, vm_min )

     if ( rank == master ) then
       write (stdout,1000) model_time, um_max, um_min
       write (stdout,2000) model_time, vm_max, vm_min
     endif

   endif

1000 format( 2x, ' model time = ', e12.5,           &
             2x, ' maximum valve for um = ', e18.9, &
             2x, ' minimum valve for um = ', e18.9 )

2000 format( 2x, ' model time = ', e12.5,           &
             2x, ' maximum valve for vm = ', e18.9, &
             2x, ' minimum valve for vm = ', e18.9  )

   stopflag = 0
   do k=2,nz-1
     do j=2,ny-1
       do i=2,nx-1
         if ( ABS(um(i,j,k)) > uvbound .OR. ABS(vm(i,j,k)) > uvbound ) then
           istop    = i
           jstop    = j
           kstop    = k
           stopflag = 1
           goto 3000
         endif
       enddo
     enddo
   enddo

3000 continue
   stopflags = stopflag
   call mpp_barrier

   call mpp_max ( stopflags, flagcheck )
   call mpp_bcast ( flagcheck )

   if ( stopflag == 1 ) then
     write (stdout, 4000) model_time
     write (stdout, 5000) rank, istop, jstop, kstop
   endif

4000 format( 2x, ' MODEL RUN TERMINATED, AT TIME  = ', e12.5 )
5000 format( 2x, ' ON MPI-RANK = ', i3, ' (i,j,k) = ', 3i3)

end subroutine model_check_bounds

subroutine model_init ( um, vm, wm, pgfx1, pgfy1, mask1 )

   real(KIND=r8), dimension(nx,ny,nz), intent(out) :: um
   real(KIND=r8), dimension(nx,ny,nz), intent(out) :: vm
   real(KIND=r8), dimension(nx,ny,nz), intent(out) :: wm
 
   real(KIND=r8), dimension(nx,ny), intent(out) :: pgfx1
   real(KIND=r8), dimension(nx,ny), intent(out) :: pgfy1
   real(KIND=r8), dimension(nx,ny), intent(out) :: mask1

   integer :: i, j, k, recs, io_unit
 
   namelist /model/ pmix_len, turh, storm_name, bot_top, &
                    check_bound_time_beg, check_bound_time_end

   call get_unit_num (io_unit)

   open (io_unit,FILE='input.nml',STATUS='UNKNOWN')
   read (io_unit,NML=model)
   close(io_unit,STATUS='KEEP')

   storm_name_lower = to_lower_case(storm_name); storm_name = storm_name_lower

   if ( rank == master ) then
     write(stdout, 1000) pmix_len, turh, storm_name, bot_top, &
                         check_bound_time_beg, check_bound_time_end
   endif
1000 format( 2x, 'mixing length, read in: model_init         =', e13.5, /, &
             2x, 'turbulence coefficient, read in model_init =', e13.5, /, &
             2x, 'storm name, read in: model_init            =', a32,   /, &
             2x, 'logical for top wind, read in: model_init  =', l6,    /, &
             2x, 'check_bound_time_beg, read in model_init   =', e13.5  /, &
             2x, 'check_bound_time_end, read in model_init   =', e13.5)


   call read_boundary_parametric( ubot, vbot, ug, vg, pgfx, pgfy, topog, mask, landzo, &
                                                               xcen, ycen, xcor, ycor, &
                                      nx, ny, nz, ngx, ngy, npx, npy, forecast_nrecs+1 )

   do recs=1,spinup_nrecs
    mask_time(recs) = (recs-1)*spinup_recs_int
   enddo

   do recs=(spinup_nrecs+1),nrecs
    mask_time(recs) = (recs-spinup_nrecs)*forecast_recs_int + spinup_nsteps
   enddo
   
!  check to see what we are using for top wind 
   if (bot_top) then 
 
     do k=2,nz
       do j=1,ny
         do i=1,nx
           um(i,j,k) = ubot(i,j,1)
           vm(i,j,k) = vbot(i,j,1)
         enddo
       enddo
     enddo

   else 

     do k=2,nz
       do j=1,ny
         do i=1,nx
           um(i,j,k) = ug(i,j,1)
           vm(i,j,k) = vg(i,j,1)
         enddo
       enddo
     enddo
 
   endif

!  using PG term during spin up
   do j=1,ny
     do i=1,nx
       pgfx1(i,j)=pgfx(i,j,1)  
       pgfy1(i,j)=pgfy(i,j,1)
     enddo
   enddo
 
   do k=1,nz
     do j=1,ny
       do i=1,nx
         wm(i,j,k) = 0.0_r8
       enddo
     enddo
   enddo
 
!  using ocean only during spin up
   do j=1,ny
     do i=1,nx
       mask1(i,j) = 0.0_r8
     enddo
   enddo

   x_center = xcen(1)
   y_center = ycen(1)

   x_corner(1) = xcor(1,1)
   x_corner(2) = xcor(1,2)
   y_corner(1) = ycor(1,1)
   y_corner(2) = ycor(1,2)

   return
 
end subroutine model_init

end module model_step_mod
