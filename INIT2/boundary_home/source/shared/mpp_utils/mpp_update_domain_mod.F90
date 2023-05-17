module mpp_update_domain_mod

   use shr_kind_mod,  only: r8 => shr_kind_r8
   use mpp_setup_mod, only: size, master, rank

   implicit none

   include 'mpif.h'

   public :: update_domain

   interface update_domain
     module procedure update_domain_3d
     module procedure update_domain_2d
   end interface

   private


   integer, parameter :: S=1, E=2, N=3, W=4
   integer, dimension(4) :: neighBor

   integer :: ierror, size_of_real
   integer :: comm2d

   logical :: first_time = .true.

contains

subroutine update_domain_3d ( x, nx, ny, nz, npx, npy )

   integer, intent(in)    :: nx, ny, nz
   integer, intent(in)    :: npx, npy

   real(KIND=r8), dimension(1:nx,1:ny,1:nz), intent(inout) :: x

   real(KIND=r8), dimension(nx*nz) :: xsend_n, xrecv_n
   real(KIND=r8), dimension(nx*nz) :: xsend_s, xrecv_s
   real(KIND=r8), dimension(ny*nz) :: xsend_w, xrecv_w
   real(KIND=r8), dimension(nz*ny) :: xsend_e, xrecv_e

   integer            :: flag, counts
   integer, dimension(mpi_status_size) :: status
   integer :: i, j, k, inc

   !--------------------------------
   ! update bounds of 3d-domain
   !--------------------------------

   ! create cartesian grid
   if (first_time ) then
     call cartesian_grid ( npx, npy )
     first_time = .false.
   endif

   !Send my boundary to North
   flag   = 1 
   counts = nx*nz
   inc    = 1
   do k=1,nz
     do i=1,nx
       xsend_n(inc) = x(i,ny-1,k); inc = inc + 1
     enddo
   enddo
   call MPI_SEND(xsend_n, counts, MPI_DOUBLE_PRECISION, neighBor(N), flag, comm2d, status, ierror)

   !Recv boundary from North
   call MPI_RECV(xrecv_n, counts, MPI_DOUBLE_PRECISION, neighBor(S), flag, comm2d, status, ierror)
   if ( neighBor(S) /= -1 ) then
     inc = 1
     do k=1,nz
       do i=1,nx
         x(i,1,k) = xrecv_n(inc); inc = inc + 1
       enddo
     enddo
   endif

   !Send my boundary to South
   flag   = 2 
   counts = nx*nz
   inc    = 1
   do k=1,nz
     do i=1,nx
       xsend_s(inc) = x(i,2,k); inc = inc + 1
     enddo
   enddo
   call MPI_SEND(xsend_s, counts, MPI_DOUBLE_PRECISION, neighBor(S), flag, comm2d, status, ierror)

   !Recv boundary from South
   call MPI_RECV(xrecv_s, counts, MPI_DOUBLE_PRECISION, neighBor(N), flag, comm2d, status, ierror)
   if ( neighBor(N) /= -1 ) then
     inc = 1
     do k=1,nz
       do i=1,nx
         x(i,ny,k) = xrecv_s(inc); inc = inc + 1
       enddo
     enddo
   endif

   !Send my boundary to West 
   flag   = 3 
   counts = ny*nz
   inc    = 1
   do k=1,nz
     do j=1,ny
       xsend_w(inc) = x(2,j,k); inc = inc + 1
     enddo
   enddo
   call MPI_SEND(xsend_w, counts, MPI_DOUBLE_PRECISION, neighBor(W), flag, comm2d, status, ierror)

   !Recv boundary from West 
   call MPI_RECV(xrecv_w, counts, MPI_DOUBLE_PRECISION, neighBor(E), flag, comm2d, status, ierror)
   if ( neighBor(E) /= -1 ) then
     inc = 1
     do k=1,nz
       do j=1,ny
         x(nx,j,k) = xrecv_w(inc); inc = inc + 1
       enddo
     enddo
   endif

   !Send my boundary to East 
   flag = 4 
   inc  = 1
   do k=1,nz
     do j=1,ny
       xsend_e(inc) = x(nx-1,j,k); inc = inc + 1
     enddo
   enddo
   call MPI_SEND(xsend_e, counts, MPI_DOUBLE_PRECISION, neighBor(E), flag, comm2d, status, ierror)

   !Recv boundary from East 
   call MPI_RECV(xrecv_e, counts, MPI_DOUBLE_PRECISION, neighBor(W), flag, comm2d, status, ierror)
   if ( neighBor(W) /= -1 ) then
     inc = 1
     do k=1,nz
       do j=1,ny
         x(1,j,k) = xrecv_e(inc); inc = inc + 1
       enddo
     enddo
   endif

end subroutine update_domain_3d

subroutine update_domain_2d ( x, nx, ny, npx, npy )

   integer, intent(in)    :: nx, ny
   integer, intent(in)    :: npx, npy

   real(KIND=r8), dimension(1:nx,1:ny), intent(inout) :: x

   real(KIND=r8), dimension(nx) :: xsend_n, xrecv_n
   real(KIND=r8), dimension(nx) :: xsend_s, xrecv_s
   real(KIND=r8), dimension(ny) :: xsend_w, xrecv_w
   real(KIND=r8), dimension(ny) :: xsend_e, xrecv_e

   integer            :: flag, counts
   integer, dimension(mpi_status_size) :: status
   integer :: i, j, inc

   !--------------------------------
   ! update bounds of 2d-domain
   !--------------------------------

   ! create cartesian grid
   if (first_time ) then
     call cartesian_grid ( npx, npy )
     first_time = .false.
   endif

   !Send my boundary to North
   flag   = 1 
   counts = nx
   inc    = 1
   do i=1,nx
     xsend_n(inc) = x(i,ny-1); inc = inc + 1
   enddo
   call MPI_SEND(xsend_n, counts, MPI_DOUBLE_PRECISION, neighBor(N), flag, comm2d, status, ierror)

   !Recv boundary from North
   call MPI_RECV(xrecv_n, counts, MPI_DOUBLE_PRECISION, neighBor(S), flag, comm2d, status, ierror)
   if ( neighBor(S) /= -1 ) then
     inc = 1
     do i=1,nx
       x(i,1) = xrecv_n(inc); inc = inc + 1
     enddo
   endif

   !Send my boundary to South
   flag   = 2 
   counts = nx
   inc    = 1
   do i=1,nx
     xsend_s(inc) = x(i,2); inc = inc + 1
   enddo
   call MPI_SEND(xsend_s, counts, MPI_DOUBLE_PRECISION, neighBor(S), flag, comm2d, status, ierror)

   !Recv boundary from South
   call MPI_RECV(xrecv_s, counts, MPI_DOUBLE_PRECISION, neighBor(N), flag, comm2d, status, ierror)
   if ( neighBor(N) /= -1 ) then
     inc = 1
     do i=1,nx
       x(i,ny) = xrecv_s(inc); inc = inc + 1
     enddo
   endif

   !Send my boundary to West 
   flag   = 3 
   counts = ny
   inc = 1
   do j=1,ny
     xsend_w(inc) = x(2,j); inc = inc + 1
   enddo
   call MPI_SEND(xsend_w, counts, MPI_DOUBLE_PRECISION, neighBor(W), flag, comm2d, status, ierror)

   !Recv boundary from West 
   call MPI_RECV(xrecv_w, counts, MPI_DOUBLE_PRECISION, neighBor(E), flag, comm2d, status, ierror)
   if ( neighBor(E) /= -1 ) then
     inc = 1
     do j=1,ny
       x(nx,j) = xrecv_w(inc); inc = inc + 1
     enddo
   endif

   !Send my boundary to East 
   flag = 4 
   inc  = 1
   do j=1,ny
     xsend_e(inc) = x(nx-1,j); inc = inc + 1
   enddo
   call MPI_SEND(xsend_e, counts, MPI_DOUBLE_PRECISION, neighBor(E), flag, comm2d, status, ierror)

   !Recv boundary from East 
   call MPI_RECV(xrecv_e, counts, MPI_DOUBLE_PRECISION, neighBor(W), flag, comm2d, status, ierror)
   if ( neighBor(W) /= -1 ) then
     inc = 1
     do j=1,ny
       x(1,j) = xrecv_e(inc); inc = inc + 1
     enddo
   endif

end subroutine update_domain_2d

subroutine cartesian_grid ( npx, npy )

   integer, intent(in) :: npx, npy

   logical, parameter    :: reorder = .false.
   integer, dimension(2) :: dims
   logical, dimension(2) :: periods
   integer :: ndims

! create 2D cartesian grid
   periods(:) = .false.
   ndims      = 2 
   dims(1)    = npy 
   dims(2)    = npx 

   call MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, comm2d, ierror)

   ! Identify neihgbors
   neighBor(:) = MPI_Proc_null

   ! Left/West and right/East neihgbors
   call MPI_Cart_shift(comm2d, 1, 1, neighBor(W), neighBor(E), ierror) 

   ! Bottom/South and Upper/North neihgbors
   call MPI_Cart_shift(comm2d, 0, 1, neighBor(S), neighBor(N), ierror) 

!  print *, rank, neighBor(W), neighBor(E), neighBor(S), neighBor(N)

end subroutine cartesian_grid

end module mpp_update_domain_mod
