module mpp_io_mod

use mpp_setup_mod, only : size, master, rank
use mpp_setup_mod, only : mpp_barrier

   implicit none
   private

   include 'mpif.h'

   public :: get_unit_num, error_handler, add_rank_to_file_name
   public :: to_upper_case, to_lower_case
   public :: read_diag_table

   integer, public :: stdout=6, stderr=5

   integer :: unit_intern=8
   logical :: first_time=.true.

contains

subroutine get_unit_num ( unit_avail )

! gets an available unit number
   integer, intent(out)  :: unit_avail

   unit_avail  = unit_intern 
   unit_intern = unit_intern + 1

end subroutine get_unit_num

subroutine error_handler ( error_message, site )

! error handler routine

   character(*), intent(in) :: error_message, site

   character(len=128) :: message_upper
   integer :: errcode, ierr

   message_upper = to_upper_case(error_message)

   if ( TRIM(message_upper)  == 'WARNING' ) then
     if ( rank == master ) then
       write (stdout,1000) TRIM(site); call flush(stdout)
     endif
     call mpp_barrier
1000 format (2x, 'WARNING', ' at call site ', a64)
     return
   endif
  
   if ( TRIM(message_upper)  == 'ERROR' ) then
     if ( rank == master ) then
       write (stdout,2000) TRIM(site); call flush(stdout)
     endif
     call mpp_barrier
2000 format (2x, 'ERROR', ' at call site: ', a64)
     call MPI_Abort(MPI_COMM_WORLD, errcode, ierr)
   endif

end subroutine error_handler

function to_upper_case(strIn) result(strOut)

! converts character string to upper case

  character(len=*), intent(in) :: strIn

  character(len=len(strIn)) :: strOut
  integer :: i, j

  do i=1,len(strIn)
    j = iachar(strIn(i:i))
    if (j>= iachar("a") .and. j<=iachar("z") ) then
         strOut(i:i) = achar(iachar(strIn(i:i))-32)
    else
         strOut(i:i) = strIn(i:i)
    end if
  enddo

end function to_upper_case

function to_lower_case(strIn) result(strOut)

! converts character string to lowercase
  character(len=*), intent(in) :: strIn
  character(len=len(strIn)) :: strOut

  character (LEN=26), parameter :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
  character (LEN=26), parameter :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  integer :: i, n

! copy input string
  strOut = strIn

! loop over string elements
  do i=1,len(strOut)
!   find location of letter in upper case constant string
    n = INDEX(UPPER_CASE, strOut(i:i))
!   if current substring is an upper case letter, make it lower case
    if ( n /= 0 ) strOut(i:i) = LOWER_CASE(n:n)
  enddo
end function to_lower_case

subroutine add_rank_to_file_name ( file_name_out )


! adds number of npi-rank to file name

   character(*), intent(inout) :: file_name_out

   character(len=48) :: format_string, file_beg, file_mod

   file_beg = file_name_out

   if ( rank <= 9 ) then
     file_mod = TRIM(file_beg) // '.000'
     format_string = "(a,i1)"
   elseif ( rank <= 99 ) then
     file_mod = TRIM(file_beg) // '.00'
     format_string = "(a,i2)"
   elseif ( rank <= 999 ) then
     file_mod = TRIM(file_beg) // '.0'
     format_string = "(a,i3)"
   elseif ( rank <= 9999 ) then
     file_mod = TRIM(file_beg) // '.'
     format_string = "(a,i4)"
   endif

   write(file_name_out, format_string) TRIM(file_mod), rank

end subroutine add_rank_to_file_name 

subroutine read_diag_table ( diag_table_list, num_diag_table )

   character(LEN=16), allocatable, dimension(:), intent(out) :: diag_table_list
   integer,                                      intent(out) :: num_diag_table
    
   integer :: io_unit, ierr
   integer :: i

   num_diag_table = 0

   open (io_unit, file='diag_table', status='old')
   do i=1,1000
     read (io_unit, 1000, iostat=ierr) 
     if ( ierr .eq. -1 ) exit
     num_diag_table = num_diag_table + 1
   enddo
   close(io_unit)

   allocate (diag_table_list(num_diag_table))

   open (io_unit, file='diag_table', status='old')
   
   if ( first_time ) then
     if ( rank == master ) then
       write(stdout,2000) 
     endif
   endif

   do i=1,num_diag_table
     read (io_unit, 1000, iostat=ierr) diag_table_list(i)
     if ( first_time ) then
       if ( rank == master ) then
         write(stdout,3000) i, TRIM(diag_table_list(i))
       endif
     endif
   enddo
   first_time = .false.
   close(io_unit)

1000 format(a16)
2000 format(2x, 'THE FOLLOWING DIAGNOSTICS ARE LISTED IN DIAG_TABLE:' )
3000 format(2x, i3, 2x, a16)

end subroutine read_diag_table

end module mpp_io_mod
