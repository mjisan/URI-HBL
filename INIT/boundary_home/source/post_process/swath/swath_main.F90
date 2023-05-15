program swath_main

   use shr_kind_mod,  only : r4 => shr_kind_r4

   use mpp_setup_mod,             only : mpp_init, mpp_end, mpp_barrier
   use mpp_setup_mod,             only : size, master, rank
   use mpp_io_mod,                only : error_handler, stdout
   use tprof_mod,                 only : tstart, tstop, tprnt

   use mpp_io_boundary_swath_mod, only : write_boundary_swath   
   use swath_mod,                 only : create_swath
   use swath_utils_mod,           only : swath_des_dims, swath_src_dims
   use swath_utils_mod,           only : swath_des_grid, swath_src_grid

   implicit none

   integer, parameter :: npx = 1
   integer, parameter :: npy = 1

   integer :: nxdes, nydes
   integer :: nxsrc, nysrc
   integer :: nrecs 

   real(KIND=r4), allocatable, dimension(:,:) :: rainfall_des
   real(KIND=r4), allocatable, dimension(:,:) :: rainfall_src

   real(KIND=r4), allocatable, dimension(:,:) :: uvbot_des
   real(KIND=r4), allocatable, dimension(:,:) :: uvbot_src

   real(KIND=r4), allocatable, dimension(:,:) :: uvbot
   real(KIND=r4), allocatable, dimension(:,:) :: rainfall

   real(KIND=r4), allocatable, dimension(:)   :: lon_des
   real(KIND=r4), allocatable, dimension(:)   :: lat_des

   real(KIND=r4), allocatable, dimension(:)   :: lon_src
   real(KIND=r4), allocatable, dimension(:)   :: lat_src

   integer :: i, j, recs

   call mpp_init( npx, npy )

!  set size of destination grid
   call swath_des_dims ( nxdes, nydes )

!  set size of source grid  
   call swath_src_dims ( nxsrc, nysrc, nrecs )

!  allocate destination and source grids

   allocate ( rainfall_des(nxdes,nydes) )
   allocate ( rainfall_src(nxsrc,nysrc) )

   allocate ( uvbot_des(nxdes,nydes) )
   allocate ( uvbot_src(nxsrc,nysrc) )

   allocate (    uvbot(nxdes,nydes) )
   allocate ( rainfall(nxdes,nydes) )

   allocate ( lon_des(nxdes) )
   allocate ( lat_des(nydes) )

   allocate ( lon_src(nxsrc) )
   allocate ( lat_src(nysrc) )
  
!  set destination grid 
   call swath_des_grid ( uvbot, rainfall, lon_des, lat_des, nxdes, nydes )

!  initialize swath fields

   call tstart('MAIN_SWATH_LOOP')

   do recs=1,nrecs 

!    set source grid 
     call swath_src_grid ( uvbot_src, rainfall_src, lon_src, lat_src, nxsrc, nysrc, recs )

     write(stdout, 1000) recs
1000 format(2x, 'recs step = ', i4)

!    create interpolated fields

     call create_swath ( uvbot_src, uvbot_des, rainfall_src, rainfall_des, &
                lon_src, lat_src, lon_des, lat_des, nxsrc, nysrc, nxdes, nydes )

!    combine interpolated fileds to create swath fields

     do j=1,nydes
       do i=1,nxdes
          uvbot(i,j)   = MAX(uvbot_des(i,j),uvbot(i,j))
         rainfall(i,j) = rainfall(i,j) + rainfall_des(i,j)
       enddo
     enddo

   enddo

  call write_boundary_swath ( uvbot, rainfall, lon_des, lat_des, nxdes, nydes )

  call tstop ('MAIN_SWATH_LOOP')
  call tprnt

  call mpp_end

contains

end program swath_main
