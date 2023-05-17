module bilinear_interp_mod

use shr_kind_mod, only : r4 => shr_kind_r4

use mpp_io_mod,   only : error_handler

  implicit none
  private

  type pair
     integer   :: lo, up
  end type pair
  
  public :: bilinear_interp, pair

contains

! biliear interpolation using binary search. avoid masked points during
! interpolation

subroutine bilinear_interp (maskin, maskout, datain, dataout, xo, yo, xn, yn, &
                                             io, jo, in, jn, ig1, jg1, flagnum)

   real(kind=r4),                intent(in)    :: maskin
   real(kind=r4),                intent(in)    :: maskout
   real(kind=r4), dimension(io), intent(in)    :: xo
   real(kind=r4), dimension(jo), intent(in)    :: yo
   real(kind=r4), dimension(in), intent(in)    :: xn
   real(kind=r4), dimension(jn), intent(in)    :: yn
   integer,                      intent(in)    :: io
   integer,                      intent(in)    :: jo
   integer,                      intent(in)    :: in
   integer,                      intent(in)    :: jn
   type(pair),                   intent(inout) :: ig1
   type(pair),                   intent(inout) :: jg1
   integer,                      intent(in)    :: flagnum
   
   real(kind=r4), intent(in), dimension(io, jo)  :: datain
   real(kind=r4), intent(out),dimension(in, jn)  :: dataout

   type(pair), dimension(in) :: ig
   type(pair), dimension(jn) :: jg 
     
   integer       :: i, j, i1, i2, j1, j2, ik, jk
   real(kind=r4) :: x, y, x1, x2, y1, y2
   real(kind=r4) :: sum, dxy
   real(kind=r4) :: flag(4)

   if ( in /= 1 .OR. jn /= 1 ) then
      call error_handler ( 'ERROR', 'bilinear_interp: in /= 1 .OR. jn /= 1')
   endif
   if (yo(1) > yo(jo) .OR. xo(1) > xo(io)) then
     call error_handler ( 'ERROR', 'bilinear_interp: input x, y are not in ascending order')
   endif

   if (ig1%lo == flagnum .OR. ig1%up == flagnum) then 
      call search_bin (xo, xn, io, in, ig)
      ig1 = ig(1)
   else
      ig(1) = ig1
   endif 
   if ( jg1%lo == flagnum .OR. jg1%up == flagnum) then
     call search_bin (yo, yn, jo, jn, jg)
     jg1 = jg(1)
   else
     jg(1) = jg1
   endif 
    
   do j=1,jn
     do i=1,in

       i1 = ig(i)%lo
       i2 = ig(i)%up
       j1 = jg(j)%lo
       j2 = jg(j)%up

 !   j=1 is south end point and j=jn is north end point 
       if (i1 == -99 .OR. i2 == -99 .OR. &
           j1 == -99 .OR. j2 == -99 ) then 
           dataout(i,j) = maskout
           cycle
       endif
       x1 = xo(i1)
       x2 = xo(i2)
       y1 = yo(j1)
       y2 = yo(j2)

       x = xn(i)
       y = yn(j)

       if ( x1 == x ) then
         x2 = x1
         i2 = i1
       endif
       if ( x1 /= x .AND. x2 == x ) then
         x1 = x2
         i1 = i2
       endif
       if ( y1 == y ) then
         y2 = y1
         j2 = j1
       endif
       if (y1 /= y .AND. y2 == y) then
         y1 = y2
         j1 = j2
       endif
  
       if (datain(i1,j1) == maskin) then
           flag(1)=0.
       else
           flag(1)=1.
       endif
       if (datain(i2,j1) == maskin) then
           flag(2)=0.
       else
           flag(2)=1.
       endif
       if (datain(i1,j2) == maskin) then
           flag(3)=0.
       else
           flag(3)=1.
       endif
       if (datain(i2,j2) == maskin) then
           flag(4)=0.
       else
           flag(4)=1.
       endif
  
       if ( x2 /= x1 .AND. y2 /= y1) then
          sum=flag(1) + flag(2) + flag(3) + flag(4)
          if (sum == 0.) then
            dataout(i,j)=maskin
          else
            dxy = flag(1)*(x2-x)*(y2-y) + flag(2)*(x-x1)*(y2-y)       & 
                + flag(3)*(x2-x)*(y-y1) + flag(4)*(x-x1)*(y-y1)
  
            dataout(i,j) = flag(1)*datain(i1,j1)*(x2-x)*(y2-y)/dxy    &
                         + flag(2)*datain(i2,j1)*(x-x1)*(y2-y)/dxy    &     
                         + flag(3)*datain(i1,j2)*(x2-x)*(y-y1)/dxy    &
                         + flag(4)*datain(i2,j2)*(x-x1)*(y-y1)/dxy
          endif
       elseif( x2 /= x1 .AND. y2 == y1) then
          sum = flag(1) + flag(2) + flag(3) + flag(4)
          if (sum == 0.) then
             dataout(i,j)=maskin
          else
             dxy = flag(1)*(x2-x) + flag(2)*(x-x1) + flag(3)*(x2-x)  &
                                                   + flag(4)*(x-x1)
             dataout(i,j) = flag(1)*datain(i1,j1)*(x2-x)/dxy         &
                          + flag(2)*datain(i2,j1)*(x-x1)/dxy         &
                          + flag(3)*datain(i1,j2)*(x2-x)/dxy         &
                          + flag(4)*datain(i2,j2)*(x-x1)/dxy
          endif
       elseif( x2 == x1 .AND. y2 /= y1) then
          sum = flag(1) + flag(2) + flag(3) + flag(4)
          if (sum == 0.) then
             dataout(i,j)=maskin
          else
             dxy = flag(1)*(y2-y) + flag(2)*(y2-y) + flag(3)*(y-y1)  &
                                                   + flag(4)*(y-y1)
             dataout(i,j) = flag(1)*datain(i1,j1)*(y2-y)/dxy         &
                          + flag(2)*datain(i2,j1)*(y2-y)/dxy         &
                          + flag(3)*datain(i1,j2)*(y-y1)/dxy         &
                          + flag(4)*datain(i2,j2)*(y-y1)/dxy
          endif
       else
          sum = flag(1)
          if (sum == 0.) then
             dataout(i,j)=maskin
          else
            dataout(i,j) = datain(i1,j1)
          endif
       endif
     
     enddo 
   enddo 
  
end subroutine bilinear_interp

 subroutine search_bin (yo, yn, jo, jn, jg)

! doing binary searching  neighboring points

   real(kind=r4), dimension(jo), intent(in)  :: yo
   real(kind=r4), dimension(jn), intent(in)  :: yn
   integer,                      intent(in)  :: jo
   integer,                      intent(in)  :: jn
   type(pair),    dimension(jn), intent(out) :: jg

   integer :: j, jk

   if (yo(1) > yo(jo) .or. yo(1) > yo(2)) then
     call error_handler ( 'ERROR', 'search_bin: inputs are not ascending' )
   endif
     
   do j=1,jn
      
      if (yn(j) == yo(1)) then
          jg(j)%lo=1
          jg(j)%up=1
      elseif (yn(j) == yo(jo)) then
          jg(j)%lo=jo
          jg(j)%up=jo
      elseif (yn(j) < yo(1)) then
          jg(j)%lo=-99
          jg(j)%up=-99
      elseif (yn(j) > yo(jo)) then
          jg(j)%lo=-99
          jg(j)%up=-99
      else
          jg(j)%lo=1
          jg(j)%up=jo
          do while( jg(j)%up - jg(j)%lo > 1 )
            call binary_search( yo, jo, yn(j), jg(j) )
          enddo
      endif
     
   enddo 
     
end subroutine search_bin

subroutine binary_search(yo, jo, y1, jg)

! doing binary searching  neighboring points

    real(kind=r4), dimension(jo), intent(in)    :: yo
    integer,                      intent(in)    :: jo
    real(kind=r4),                intent(in)    :: y1   
    type(pair),                   intent(inout) :: jg

    integer :: j, jk, jm

    if ( y1 > yo(jo) .OR. y1 < yo(1) ) then
      call error_handler ( 'ERROR', 'binary_search: xinp is out of range' )
    endif

    if ( jg%lo <= jg%up ) then
      jm = INT( 0.5 *( jg%lo + jg%up ) )
      if ( yo(jm) >= y1 ) then
           jg%up = jm
      else
           jg%lo = jm
      endif
    else
      call error_handler ( 'ERROR', 'binary_search: data is not ordered' )
    endif
        
end subroutine binary_search
      
end module bilinear_interp_mod
