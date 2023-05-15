module linear_interp_mod

use shr_kind_mod, only : r4 => shr_kind_r4
use shr_kind_mod, only : r8 => shr_kind_r8

use mpp_io_mod,   only : error_handler
implicit none
private

public :: linear_interp

contains 

subroutine linear_interp ( maskin, maskout, datain, dataout, xo, xn, io, in ) 

! performs linear interpolation; masked values are excluded

  integer,                       intent(in)  :: io, in 
  real(KIND=r4) ,                intent(in)  :: maskin, maskout
  real(KIND=r4) , dimension(io), intent(in)  :: xo, datain
  real(KIND=r4) , dimension(in), intent(in)  :: xn
  real(KIND=r4) , dimension(in), intent(out) :: dataout

  integer                      :: i, i1, i2, ik
  real(KIND=r4)                :: x, x1, x2, dx, sumval
  real(KIND=r4) , dimension(2) :: flag
 
  if ( xo(1) > xo(io)) then
    call error_handler ( 'ERROR', 'linear_interp' )
  endif

  do i=1,in
    if ( xn(i) > xo(io) ) then
      dataout(i) = datain(io)
      cycle
    endif
    if ( xn(i) < xo(1) ) then
      dataout(i) = datain(1)
      cycle
    endif

    x=xn(i)
    i1=1
    i2=1 
    x1=xo(i1)
    x2=xo(i2)
    do ik=1,io
      if ( xo(ik) <= x ) then
        i1 = ik
        x1=xo(ik)
      endif
    enddo
    do ik=io,1,-1
      if ( xo(ik) >= x ) then
        i2 = ik
        x2=xo(ik)
      endif
    enddo
       
    if (datain(i1) == maskin) then
      flag(1)=0.0_r4
    else
      flag(1)=1.0_r4
    endif
    if (datain(i2) == maskin) then
      flag(2)=0.0_r4
    else
      flag(2)=1.0_r4
    endif

    if (x2 /= x1) then
      sumval=flag(1) + flag(2)
      if (sumval .eq. 0.) then
        dataout(i)=maskout
      else
        dx = flag(1)*(x2-x) + flag(2)*(x-x1)

!--->> y = [ y1*(x2-x) + y2*(x-x1) ] / [x2-x1]
        dataout(i) = flag(1)*datain(i1)*(x2-x)/dx + &
                     flag(2)*datain(i2)*(x-x1)/dx         
      endif
    else
      sumval=flag(1)
      if (sumval == 0.0_r4) then
        dataout(i)=maskout
      else
        dataout(i) = datain(i1)
      endif
    endif
      
  enddo

end subroutine linear_interp

end module linear_interp_mod
