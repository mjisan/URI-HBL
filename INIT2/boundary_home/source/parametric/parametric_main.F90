program parametric_main

   use shr_kind_mod,                   only : r4 => shr_kind_r4
   use shr_kind_mod,                   only : r8 => shr_kind_r8
   use mpp_setup_mod,                  only : mpp_init, mpp_end, mpp_barrier
   use mpp_setup_mod,                  only : size, master, rank
   use mpp_io_mod,                     only : error_handler, stdout
   use mpp_io_boundary_parametric_mod, only : write_boundary_parametric
   use tprof_mod,                      only : tstart, tstop, tnull, tprnt
   
   use storm_stats_exp_mod,            only : storm_stats_exp
   use topography_mod,                 only : topography
   use land_roughness_mod,             only : land_roughness
   use top_wind_mod,                   only : norm_dist_ratio
   use press_grad_mod,                 only : press_grad
   use bilinear_interp_mod,            only : bilinear_interp, pair
   use grid_mod,                       only : nx, ny, npx, npy, icen, jcen
   use grid_mod,                       only : stepx, stepy
   use time_mod,                       only : nrecs, delt

   implicit none

   real(KIND=r4), dimension(nx,ny) :: dist
   real(KIND=r4), dimension(nx,ny) :: pres_sea
   real(KIND=r4), dimension(nx,ny) :: ubot
   real(KIND=r4), dimension(nx,ny) :: vbot
   real(KIND=r4), dimension(nx,ny) :: utop
   real(KIND=r4), dimension(nx,ny) :: vtop
   real(KIND=r4), dimension(nx,ny) :: pgfx
   real(KIND=r4), dimension(nx,ny) :: pgfy
   real(KIND=r4), dimension(nx,ny) :: topog 
   real(KIND=r4), dimension(nx,ny) :: land_mask
   real(KIND=r4), dimension(nx,ny) :: land_rough
   real(KIND=r4), dimension(nx)    :: xgrid
   real(KIND=r4), dimension(ny)    :: ygrid
   character(LEN=32)               :: storm_name
   character(LEN=32)               :: storm_start_time

   real(KIND=r4) :: xcen, ycen
   real(KIND=r4), dimension(2) :: xcor, ycor
   real(KIND=r4) :: utx, uty, rmax, wsmax
   real(KIND=r4) :: time_sec, time_hrs, time_day

   integer :: i, j, n, recs

   type(pair) :: ig1, jg1

   call mpp_init( nx, ny)

   call create_parametric_grid (xgrid, ygrid, dist)

   call tstart('MAIN_PARAMETRIC_LOOP')

   do recs=1,nrecs 

     time_sec = (recs-1)*delt
     time_hrs = time_sec/3600.0
     time_day = time_hrs/24.0

     write(stdout, 1000) recs, time_hrs
1000 format(2x, 'recs step = ', i4, 5x, 'time in hours = ', f6.2)

! generate wind on parametric domain
     call storm_stats_exp (xgrid, ygrid, time_day, ubot, vbot, pres_sea, &
                  xcen, ycen, xcor, ycor, utx, uty, rmax, wsmax, & 
                                           storm_name, storm_start_time)

! generate topography on parametric domain
     call topography (topog, land_mask, xgrid, ygrid, xcen, ycen, &
                                                        storm_name)

! generate land roughness on parametric domain
     call land_roughness (land_rough, land_mask, xgrid, ygrid, xcen, ycen, &
                                                                 storm_name)

! calculating parametric top winds on parametric domain
     call norm_dist_ratio (xgrid, ygrid, dist, ubot, vbot, utx, uty, &
                                                   rmax, utop, vtop)

! calculating pressure gradient terms on parametric domain
     call press_grad (utop, vtop, pgfx, pgfy)

! write NetCDF parametric file
     call write_boundary_parametric (ubot, vbot, utop, vtop, pgfx, pgfy, &
                               topog, land_mask, land_rough, pres_sea, &
                        xcen, ycen, xcor, ycor, utx, uty, rmax, wsmax, &
                               time_day, storm_start_time, recs, nx, ny)

  enddo

  call tstop ('MAIN_PARAMETRIC_LOOP')

  call tprnt

  call mpp_end

contains

subroutine create_parametric_grid ( xgrid, ygrid, dist )

   real(KIND=r4), dimension(nx),    intent(out) ::  xgrid
   real(KIND=r4), dimension(ny),    intent(out) ::  ygrid
   real(KIND=r4), dimension(nx,ny), intent(out) ::   dist

   integer :: i, j

   do i=1+icen,nx
     xgrid(i) =  (i-icen) * stepx
   enddo

   do i=1,icen
     xgrid(i) = - xgrid(nx) + (i-1) * stepx
   enddo

   do j= 1+jcen,ny
     ygrid(j) =  (j-jcen) * stepy
   enddo

   do j=1,jcen
     ygrid(j) = - ygrid(ny) + (j-1) * stepy
   enddo

   do j=1,ny
     do i = 1, nx
      dist(i,j) = SQRT(xgrid(i) * xgrid(i) + ygrid(j) * ygrid(j))
     enddo
   enddo

end subroutine create_parametric_grid

end program parametric_main
