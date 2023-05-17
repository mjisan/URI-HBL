module mpp_setup_mod

   implicit none
   private

   include 'mpif.h'

   public :: mpp_init, mpp_end, mpp_barrier
   public :: size, master, rank
   public :: npx, npy

   integer :: size, rank
   integer :: master = 0 ! master MPI-Rank 
   integer :: npx, npy

contains

subroutine mpp_init (npx_in, npy_in)

   integer, intent(IN) :: npx_in, npy_in
   integer :: ierror, size_out, rank_out

   call MPI_INIT(ierror)

   call MPI_COMM_SIZE(MPI_COMM_WORLD, size_out, ierror)
   call MPI_COMM_RANK(MPI_COMM_WORLD, rank_out, ierror)

   size = size_out; rank = rank_out
   npx = npx_in; npy = npy_in
   
end subroutine mpp_init

subroutine mpp_end 

   integer :: ierror

   call MPI_Finalize(ierror)

end subroutine mpp_end 

subroutine mpp_barrier
   integer :: ierror

   call MPI_BARRIER(MPI_COMM_WORLD, ierror)

end subroutine mpp_barrier

end module mpp_setup_mod

