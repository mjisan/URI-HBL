module tprof_mod

use mpp_io_mod, only : stdout, error_handler

  implicit none
  private

  public :: tstart, tstop, tnull, tprnt
 
  integer, parameter :: max_clocks = 100 
  integer            :: num_clocks = 0

  type clock_type
    integer*8          :: clock_time
    integer*8          :: clock_total
    character (len=32) :: clock_name
  end type clock_type

  type (clock_type), dimension(max_clocks) :: clock_details

  logical   :: tprof_init = .false.
  logical   :: tprof_null = .false.
  integer*8 :: start_clocks, stop_clocks, rate_clock

  integer :: block_cnt
  integer :: block_pos
  logical :: block_set= .false.

!---- version number -----
  character(len=128) :: version = ''
  character(len=128) :: tagname = ''

contains

subroutine tinit
  integer :: n
 
  if ( tprof_init ) return

  do n =1,max_clocks
    clock_details(n)%clock_time  = 0
    clock_details(n)%clock_total = 0
    clock_details(n)%clock_name  = ''
  enddo

  tprof_init = .true.

end subroutine tinit

subroutine tnull 

  tprof_null = .true. 

end subroutine tnull 

subroutine tstart( block_name, block_num )
  character(len=*),  intent(in) :: block_name
  integer, optional, intent(in) :: block_num

  character(len=32) :: newStr
  integer           :: n

  if ( tprof_null ) return
  if ( .NOT. tprof_init ) call tinit( )

  newStr = TRIM( to_upper(block_name) )

! assign first clock
  if ( num_clocks == 0 ) then
    num_clocks = num_clocks + 1

    if ( present (block_num) .and. .NOT. block_set ) then
      block_cnt = block_num
      block_pos = num_clocks
      block_set = .true. 
    endif

    clock_details(num_clocks)%clock_name = newStr
    call system_clock(start_clocks, rate_clock)
    clock_details(num_clocks)%clock_time = start_clocks
    return
  endif

! assign if clock exists 
  do n =1,num_clocks   
    if ( newStr == clock_details(n)%clock_name ) then
      call system_clock(start_clocks, rate_clock)
      clock_details(n)%clock_time = start_clocks
      return
    endif
  enddo

! assign new clock as one does not exist
  num_clocks = num_clocks + 1

  if ( present (block_num) .and. .NOT. block_set ) then
    block_cnt = block_num
    block_pos = num_clocks
    block_set = .true. 
  endif

  clock_details(num_clocks)%clock_name = newStr
  call system_clock(start_clocks, rate_clock)
  clock_details(num_clocks)%clock_time = start_clocks 

  if ( num_clocks+1 == max_clocks ) then
    call error_handler ( 'Increase number of clocks available', 'tprof_mod/tstart' ) 
  endif

end subroutine tstart

subroutine tstop( block_name )
  character(len=*), intent(in) :: block_name
  character(len=32) :: newStr
  integer           :: n

  if ( tprof_null ) return

  newStr = TRIM( to_upper(block_name) )

  do n =1,num_clocks
    if ( newStr == clock_details(n)%clock_name ) then
      call system_clock(stop_clocks, rate_clock)
      clock_details(n)%clock_time  = stop_clocks - clock_details(n)%clock_time
      clock_details(n)%clock_total = clock_details(n)%clock_total + clock_details(n)%clock_time
      return
    endif
  enddo

  call error_handler ( 'Name of region not assigned', 'tprof_mod/tstop' ) 

end subroutine tstop

subroutine tprnt 

  integer :: n

  if( tprof_null ) then
    write(stdout,1000)
1000 format( ' TPROF flag was set to null ')
  else
    do n =1,num_clocks
      write(stdout,2000) TRIM(clock_details(n)%clock_name), clock_details(n)%clock_total/REAL(rate_clock)
2000 format( ' Region Name = ', a32, ' : Time (secs) = ', e13.6) 
    enddo
  endif

  if ( block_set ) then
    n = block_pos 
    write(stdout,3000) (clock_details(n)%clock_total/REAL(rate_clock)/block_cnt)
3000 format( ' Time per iteration (secs) = ', e13.6) 
  endif

end subroutine tprnt

function to_upper(strIn) result(strOut)

     character(len=*), intent(in) :: strIn
     character(len=len(strIn)) :: strOut
     integer :: i, j

     do i = 1,len(strIn)
          j = iachar(strIn(i:i))
          if (j>= iachar("a") .and. j<=iachar("z") ) then
               strOut(i:i) = achar(iachar(strIn(i:i))-32)
          else
               strOut(i:i) = strIn(i:i)
          end if
     end do

end function to_upper

end module tprof_mod
