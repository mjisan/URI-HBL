module swath_mod

use shr_kind_mod,        only : r4 => shr_kind_r4
use shr_kind_mod,        only : r8 => shr_kind_r8

use bilinear_interp_mod, only : bilinear_interp, pair

implicit none

private

integer, parameter :: flagint  = -9999
real(KIND=r4)      :: flagreal = -9999_r4

type(pair) :: ig1, jg1

public :: create_swath

real(KIND=r8) :: deltx, delty

real(KIND=r4), allocatable, dimension(:,:) :: field_srcf
real(KIND=r4), allocatable, dimension(:)   :: lons_srcf
real(KIND=r4), allocatable, dimension(:)   :: lats_srcf 

contains

subroutine create_swath ( uvbot_src, uvbot_des, rainfall_src, rainfall_des, &
                                    lons_src, lats_src, lons_des, lats_des, &
                                                 nxsrc, nysrc, nxdes, nydes )

   real(KIND=r4), dimension(nxsrc,nysrc), intent(in)  :: uvbot_src
   real(KIND=r4), dimension(nxdes,nydes), intent(out) :: uvbot_des
   real(KIND=r4), dimension(nxsrc,nysrc), intent(in)  :: rainfall_src
   real(KIND=r4), dimension(nxdes,nydes), intent(out) :: rainfall_des
   real(KIND=r4), dimension(nxsrc),       intent(in)  :: lons_src
   real(KIND=r4), dimension(nysrc),       intent(in)  :: lats_src
   real(KIND=r4), dimension(nxdes),       intent(in)  :: lons_des
   real(KIND=r4), dimension(nydes),       intent(in)  :: lats_des
   integer,                               intent(in)  :: nxsrc
   integer,                               intent(in)  :: nysrc
   integer,                               intent(in)  :: nxdes
   integer,                               intent(in)  :: nydes

   integer :: i, j
   integer :: nxsrcf, nysrcf
   integer :: nxwest, nxeast, nysouth, nynorth

   call create_swath_init ( lons_src, lats_src, lons_des, lats_des, &
                         nxsrc, nysrc, nxsrcf, nysrcf, nxdes, nydes )

!  create indices to map from local to global src grid

   nxwest   = NINT(ABS((lons_des(1) - lons_src(1)))/deltx) + 1
   nxeast   = nxwest + nxsrc 
   nysouth  = NINT(ABS((lats_des(1) - lats_src(1)))/delty) + 1
   nynorth  = nysouth + nysrc

!  create full src grid for uvbot

   do j=1,nysrcf
     do i=1,nxsrcf
       field_srcf(i,j) = flagreal
     enddo
   enddo

   do j=nysouth,nynorth
     do i=nxwest,nxeast
       field_srcf(i,j) = uvbot_src(i-nxwest+1,j-nysouth+1)
     enddo
   enddo 

   do j=1,nydes
     do i=1,nxdes
       ig1%lo = flagint
       ig1%up = flagint
       jg1%lo = flagint
       jg1%up = flagint

       call bilinear_interp (flagreal, flagreal, field_srcf, uvbot_des(i,j), lons_srcf, lats_srcf, &
                                  lons_des(i), lats_des(j), nxsrcf, nysrcf, 1, 1, ig1, jg1, flagint)
     enddo
   enddo

   do j=1,nydes
     do i=1,nxdes
       if ( uvbot_des(i,j) == flagreal) then
         uvbot_des(i,j) = 0.0
       endif
     enddo
   enddo

!  create full src grid for rainfall

   do j=1,nysrcf
     do i=1,nxsrcf
       field_srcf(i,j) = flagreal
     enddo
   enddo

   do j=nysouth,nynorth
     do i=nxwest,nxeast
       field_srcf(i,j) = rainfall_src(i-nxwest+1,j-nysouth+1)
     enddo
   enddo 

   do j=1,nydes
     do i=1,nxdes
       ig1%lo = flagint
       ig1%up = flagint
       jg1%lo = flagint
       jg1%up = flagint

       call bilinear_interp (flagreal, flagreal, field_srcf, rainfall_des(i,j), lons_srcf, lats_srcf, &
                                     lons_des(i), lats_des(j), nxsrcf, nysrcf, 1, 1, ig1, jg1, flagint)
     enddo
   enddo

   do j=1,nydes
     do i=1,nxdes
       if ( rainfall_des(i,j) == flagreal) then
         rainfall_des(i,j) = 0.0
       endif
     enddo
   enddo

   deallocate ( field_srcf )
   deallocate ( lons_srcf )
   deallocate ( lats_srcf )

end subroutine create_swath

subroutine create_swath_init ( lons_src, lats_src, lons_des, lats_des, &
                            nxsrc, nysrc, nxsrcf, nysrcf, nxdes, nydes )

   real(KIND=r4), dimension(nxsrc),       intent(in)  :: lons_src
   real(KIND=r4), dimension(nysrc),       intent(in)  :: lats_src
   real(KIND=r4), dimension(nxdes),       intent(in)  :: lons_des
   real(KIND=r4), dimension(nydes),       intent(in)  :: lats_des
   integer,                               intent(in)  :: nxsrc
   integer,                               intent(in)  :: nysrc
   integer,                               intent(out) :: nxsrcf
   integer,                               intent(out) :: nysrcf
   integer,                               intent(in)  :: nxdes
   integer,                               intent(in)  :: nydes

   integer :: i, j

   deltx = ABS((lons_src(2) - lons_src(1)))
   delty = ABS((lats_src(2) - lats_src(1)))

   nxsrcf = ABS((lons_des(nxdes) - lons_des(1)))/deltx + 1
   nysrcf = ABS((lats_des(nydes) - lats_des(1)))/delty + 1

   allocate (  lons_srcf(nxsrcf) )
   allocate (  lats_srcf(nysrcf) )
   allocate ( field_srcf(nxsrcf,nysrcf) )

   do i=1,nxsrcf
     lons_srcf(i) = lons_des(1) + deltx*(i-1)
   enddo

   do j=1,nysrcf
     lats_srcf(j) = lats_des(1) + delty*(j-1)
   enddo

end subroutine create_swath_init

end module swath_mod
