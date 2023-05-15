module mpp_io_rolls_mod

use shr_kind_mod,  only : r8 => shr_kind_r8
use shr_kind_mod,  only : r4 => shr_kind_r4

use mpp_setup_mod, only : rank, master
use mpp_setup_mod, only : npx, npy
use mpp_io_mod,    only : error_handler, stdout, add_rank_to_file_name
use mpp_io_mod,    only : read_diag_table

   implicit none
   private

   include 'netcdf.inc'

   public :: write_rolls

   logical :: first_time = .true.

   character(LEN=16), allocatable, dimension(:) :: diag_table_list
   integer :: num_diag_table

   character(LEN=48)  :: rolls_file='OUTPUT/rolls.nc' ! name of rolls file
   character(LEN=48)  :: file_name                    ! name of rolls file and mpi-rank

   integer, parameter :: ndims = 5  ! number of dimensions in netCDF file
   integer            :: start(ndims), count(ndims)

!  when we create netCDF files, variables and dimensions, we get back an ID for
!  each

   integer, dimension(ndims) :: dimids


   character(LEN=32), parameter :: units = 'units'

   character(LEN=32), parameter :: x_name   = 'x' ! name of x axis
   character(LEN=32), parameter :: z_name   = 'z' ! name of z axis
   character(LEN=32), parameter :: i_name   = 'i' ! name of i axis
   character(LEN=32), parameter :: j_name   = 'j' ! name of j axis
   character(LEN=32), parameter :: rec_name = 't' ! name of time axis

   integer :: ncid, x_dimid, z_dimid, i_dimid, j_dimid, rec_dimid
   integer :: x_varid, z_varid, i_varid, j_varid

   integer :: u_varid
   integer :: v_varid
   integer :: w_varid
   integer :: r_varid
   integer :: p_varid
   integer :: t_varid
   integer :: b_varid

contains

subroutine write_rolls (u_roll, v_roll, w_roll, rotor_roll, psi_roll, teta_roll, blag_roll, &
                                                                     nr, nx, nz, list, recs )

   real(KIND=r8), dimension(nr,nx,nz), intent(in) ::     u_roll
   real(KIND=r8), dimension(nr,nx,nz), intent(in) ::     v_roll
   real(KIND=r8), dimension(nr,nx,nz), intent(in) ::     w_roll
   real(KIND=r8), dimension(nr,nx,nz), intent(in) :: rotor_roll
   real(KIND=r8), dimension(nr,nx,nz), intent(in) ::   psi_roll
   real(KIND=r8), dimension(nr,nx,nz), intent(in) ::  teta_roll
   real(KIND=r8), dimension(nr,nx,nz), intent(in) ::  blag_roll

   integer,                    intent(in) :: nr
   integer,                    intent(in) :: nx
   integer,                    intent(in) :: nz
   integer, dimension(nr,2,2), intent(in) :: list
   integer,                    intent(in) :: recs

   real(KIND=r4), dimension(nr,nx,nz) :: var_r4
   real(kind=r8), dimension(nx)       :: xaxis
   real(kind=r8), dimension(nz)       :: zaxis

   integer, dimension(nr) :: iaxis
   integer, dimension(nr) :: jaxis

!  loop indexes, and error handling

   integer :: i, k, m, n, retval, rec, id

   if ( first_time == .true. ) then

     if ( rank == master ) then
       write(stdout, 1000) rolls_file
1000 format( 2x, 'name of rolls file written in: write_rolls = ', a48 )
     endif

     if ( npx*npy == 1 ) then
       file_name = rolls_file
     else
       file_name = rolls_file
       call add_rank_to_file_name(file_name)
     endif

!    read diag_table
     call read_diag_table ( diag_table_list, num_diag_table )

!    define axis 
     do i=1,nx
       xaxis(i) = i
     end do
     do k=1,nz
       zaxis(k) = k
     end do

     do n=1,nr
       iaxis(n) = list(n,1,1)
       jaxis(n) = list(n,2,1)
     end do

!    create netCDF file, overwrite if already exists

     call write_rolls_state(1)

!    define the dimensions for each axis and get back ID 

     retval = nf_def_dim(ncid, x_name  , nx, x_dimid  )
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_rolls/nf_def_dim/x_dimid' )

     retval = nf_def_dim(ncid, z_name  , nz, z_dimid  )
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_rolls/nf_def_dim/z_dimid' )

     retval = nf_def_dim(ncid, i_name  , nr, i_dimid  )
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_rolls/nf_def_dim/i_dimid' )

     retval = nf_def_dim(ncid, j_name  , nr, j_dimid  )
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_rolls/nf_def_dim/j_dimid' )

     retval = nf_def_dim(ncid, rec_name, NF_UNLIMITED, rec_dimid)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_rolls/nf_def_dim/rec_dimid' )

!    define the coordinate variables

     retval = nf_def_var(ncid, i_name, NF_DOUBLE, 1, i_dimid, i_varid)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_rolls/nf_def_var/i_varid' )

     retval = nf_def_var(ncid, j_name, NF_DOUBLE, 1, j_dimid, j_varid)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_rolls/nf_def_var/j_varid' )

     retval = nf_def_var(ncid, x_name, NF_DOUBLE, 1, x_dimid, x_varid)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_rolls/nf_def_var/x_varid' )

     retval = nf_def_var(ncid, z_name, NF_DOUBLE, 1, z_dimid, z_varid)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_rolls/nf_def_var/z_varid' )

!    the dimids array is used to pass the ID of the dimensions of the variables. 
!    note that in fortran arrays are stored in column-major

     dimids(1) = x_dimid
     dimids(2) = z_dimid
     dimids(3) = i_dimid
     dimids(4) = rec_dimid

!    define the attributes for variables

     do i=1,num_diag_table

       if ( TRIM(diag_table_list(i)) == 'u_roll') then
         retval = nf_def_var(ncid, 'u_roll',     NF_REAL, ndims, dimids, u_varid )
         if (retval .NE. 0) call error_handler ( 'ERROR', 'write_rolls/nf_def_var/u_varid' )
       endif

       if ( TRIM(diag_table_list(i)) == 'v_roll') then
         retval = nf_def_var(ncid, 'v_roll',     NF_REAL, ndims, dimids, v_varid )
         if (retval .NE. 0) call error_handler ( 'ERROR', 'write_rolls/nf_def_var/v_varid' )
       endif

       if ( TRIM(diag_table_list(i)) == 'w_roll') then
         retval = nf_def_var(ncid, 'w_roll',     NF_REAL, ndims, dimids, w_varid )
         if (retval .NE. 0) call error_handler ( 'ERROR', 'write_rolls/nf_def_var/w_varid' )
       endif

       if ( TRIM(diag_table_list(i)) == 'rotor_roll') then
         retval = nf_def_var(ncid, 'rotor_roll', NF_REAL, ndims, dimids, r_varid )
         if (retval .NE. 0) call error_handler ( 'ERROR', 'write_rolls/nf_def_var/r_varid' )
       endif

       if ( TRIM(diag_table_list(i)) == 'psi_roll') then
         retval = nf_def_var(ncid, 'psi_roll',   NF_REAL, ndims, dimids, p_varid )
         if (retval .NE. 0) call error_handler ( 'ERROR', 'write_rolls/nf_def_var/p_varid' )
       endif

       if ( TRIM(diag_table_list(i)) == 'teta_roll') then
         retval = nf_def_var(ncid, 'teta_roll',  NF_REAL, ndims, dimids, t_varid )
         if (retval .NE. 0) call error_handler ( 'ERROR', 'write_rolls/nf_def_var/t_varid' )
       endif

       if ( TRIM(diag_table_list(i)) == 'blag_roll') then
         retval = nf_def_var(ncid, 'blag_roll',  NF_REAL, ndims, dimids, b_varid )
         if (retval .NE. 0) call error_handler ( 'ERROR', 'write_rolls/nf_def_var/b_varid' )
       endif

     enddo

!    end define mode; this tells netCDF we are done defining metadata

     retval = nf_enddef(ncid)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_rolls/nf_enddef' )

!    put axis data

     retval = nf_put_var_DOUBLE(ncid, x_varid, xaxis)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_rolls/nf_put_var_DOUBLE/xaxis' )

     retval = nf_put_var_DOUBLE(ncid, z_varid, zaxis)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_rolls/nf_put_var_DOUBLE/zaxis' )

     retval = nf_put_var_DOUBLE(ncid, i_varid, iaxis)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_rolls/nf_put_var_DOUBLE/iaxis' )

     retval = nf_put_var_DOUBLE(ncid, j_varid, jaxis)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_rolls/nf_put_var_DOUBLE/jaxis' )

     first_time = .false.
   else
!  open netCDF file
     call write_rolls_state(2)
   endif

   count(1) = nx
   count(2) = nz
   count(3) = nr
   count(4) = 1 
   start(1) = 1 
   start(2) = 1 
   start(3) = 1 
   start(4) = recs 

!  put variable data 

   do i=1,num_diag_table

     if ( TRIM(diag_table_list(i)) == 'u_roll') then
       call transpose_r8_r4 (u_roll, var_r4, nr, nx, nz)
       retval = nf_put_vara_REAL(ncid, u_varid, start, count, var_r4)
       if (retval .NE. 0) call error_handler ( 'ERROR', 'write_rolls/nf_put_var_REAL/u_roll' )
     endif

     if ( TRIM(diag_table_list(i)) == 'v_roll') then
       call transpose_r8_r4 (v_roll, var_r4, nr, nx, nz)
       retval = nf_put_vara_REAL(ncid, v_varid, start, count, var_r4)
       if (retval .NE. 0) call error_handler ( 'ERROR', 'write_rolls/nf_put_var_REAL/v_roll' )
     endif

     if ( TRIM(diag_table_list(i)) == 'w_roll') then
       call transpose_r8_r4 (w_roll, var_r4, nr, nx, nz)
       retval = nf_put_vara_REAL(ncid, w_varid, start, count, var_r4)
       if (retval .NE. 0) call error_handler ( 'ERROR', 'write_rolls/nf_put_var_REAL/w_roll' )
     endif

     if ( TRIM(diag_table_list(i)) == 'rotor_roll') then
       call transpose_r8_r4 (rotor_roll, var_r4, nr, nx, nz)
       retval = nf_put_vara_REAL(ncid, r_varid, start, count, var_r4)
       if (retval .NE. 0) call error_handler ( 'ERROR', 'write_rolls/nf_put_var_REAL/rotor_roll' )
     endif

     if ( TRIM(diag_table_list(i)) == 'psi_roll') then
       call transpose_r8_r4 (psi_roll, var_r4, nr, nx, nz)
       retval = nf_put_vara_REAL(ncid, p_varid, start, count, var_r4)
       if (retval .NE. 0) call error_handler ( 'ERROR', 'write_rolls/nf_put_var_REAL/psi_roll' )
     endif

     if ( TRIM(diag_table_list(i)) == 'teta_roll') then
       call transpose_r8_r4 (teta_roll, var_r4, nr, nx, nz)
       retval = nf_put_vara_REAL(ncid, t_varid, start, count, var_r4)
       if (retval .NE. 0) call error_handler ( 'ERROR', 'write_rolls/nf_put_var_REAL/teta_roll' )
     endif

     if ( TRIM(diag_table_list(i)) == 'blag_roll') then
       call transpose_r8_r4 (blag_roll, var_r4, nr, nx, nz)
       retval = nf_put_vara_REAL(ncid, b_varid, start, count, var_r4)
       if (retval .NE. 0) call error_handler ( 'ERROR', 'write_rolls/nf_put_var_REAL/blag_roll' )
     endif

   enddo

   call write_rolls_state(3)

end subroutine write_rolls

subroutine write_rolls_state(case)
   integer, intent(in) :: case
   integer :: retval

   select CASE(case)
     CASE(1)
!    create netCDF file, overwrite if already exists
     retval = nf_create(file_name, nf_clobber, ncid)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_rolls_state/nf_create' )

     CASE(2)
!    open netCDF file
     retval = nf_open(file_name, nf_write, ncid)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_rolls_state/nf_open' )

     CASE(3)
!    close netCDF file
     retval = nf_close(ncid)
     if (retval .NE. 0) call error_handler ( 'ERROR', 'write_rolls_state/nf_close' )

    end select

end subroutine write_rolls_state

subroutine transpose_r8_r4 (var_r8, var_r4, nr, nx, nz)

   real(KIND=r8), dimension(nr,nx,nz), intent(IN)  :: var_r8
   real(KIND=r4), dimension(nx,nz,nr), intent(OUT) :: var_r4
   integer,                            intent(IN)  :: nr, nx, nz

   integer :: i, k, n

   do n=1,nr
     do k=1,nz
       do i=1,nx
         var_r4(i,k,n) = var_r8(n,i,k)
       enddo
     enddo
   enddo

end subroutine transpose_r8_r4

end module mpp_io_rolls_mod
