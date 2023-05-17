module roughness_mod

use shr_kind_mod,  only : r8 => shr_kind_r8

implicit none
private

public :: roughness_0, roughness_1, roughness_2

contains

subroutine roughness_2 ( uref, roughnessm )

   real(KIND=r8), intent(in)  :: uref
   real(KIND=r8), intent(out) :: roughnessm

   real(KIND=r8), parameter :: bs0 = -8.3672761723972770e-12
   real(KIND=r8), parameter :: bs1 =  1.7398510865876079e-09
   real(KIND=r8), parameter :: bs2 = -1.3318965783633590e-07
   real(KIND=r8), parameter :: bs3 =  4.5070552944387270e-06
   real(KIND=r8), parameter :: bs4 = -6.5086768819069140e-05
   real(KIND=r8), parameter :: bs5 =  0.00044745137674732834
   real(KIND=r8), parameter :: bs6 = -0.0010745704660847233

   real(KIND=r8), parameter :: cf0 =  2.1151080765239772e-13
   real(KIND=r8), parameter :: cf1 = -3.2260663894433345e-11
   real(KIND=r8), parameter :: cf2 = -3.3297059587519610e-10
   real(KIND=r8), parameter :: cf3 =  1.7648562021709124e-07
   real(KIND=r8), parameter :: cf4 =  7.1076368256941820e-06
   real(KIND=r8), parameter :: cf5 = -0.0013914681964973246
   real(KIND=r8), parameter :: cf6 =  0.0406766967657759

   IF ( uref .LE. 5.0_r8 ) THEN
     roughnessm = (0.0185_r8 / 9.8_r8*(7.59e-4*uref**2+2.46e-2*uref)**2)
   elseif (uref .GT. 5.0_r8 .AND. uref .LT. 10.0_r8) THEN
     roughnessm =.00000235*(uref**2 - 25.0_r8 ) + 3.805129199617346e-05
   elseif ( uref .GE. 10.0_r8  .AND. uref .LT. 60.0_r8) THEN
     roughnessm = bs6 + bs5*uref + bs4*uref**2 + bs3*uref**3 + bs2*uref**4   &
          + bs1*uref**5 + bs0*uref**6
   else
     roughnessm = cf6 + cf5*uref + cf4*uref**2 + cf3*uref**3 + cf2*uref**4   &
          + cf1*uref**5 + cf0*uref**6
   endif 

end subroutine roughness_2

subroutine roughness_1 ( uref, roughnessm )

   real(KIND=r8), intent(in)  :: uref
   real(KIND=r8), intent(out) :: roughnessm

   real(KIND=r8), parameter :: yz =  0.0001344_r8
   real(KIND=r8), parameter :: y1 =  3.015e-05
   real(KIND=r8), parameter :: y2 =  1.517e-06
   real(KIND=r8), parameter :: y3 = -3.567e-08
   real(KIND=r8), parameter :: y4 =  2.046e-10

   if ( uref .LT. 12.5_r8 ) THEN
     roughnessm  = (0.0185_r8 / 9.8_r8*(7.59e-4*uref**2+2.46e-2*uref)**2)
   else if ( uref .GE. 12.5_r8 .AND. uref .LT. 30.0_r8 ) THEN
     roughnessm = (0.0739793_r8 * uref -0.58_r8)/1000.0_r8
   else
     roughnessm = yz + uref*y1 + uref**2*y2 + uref**3*y3 + uref**4*y4
   endif 

end subroutine roughness_1

subroutine roughness_0 (uref, roughnessm )

   real(KIND=r8), intent(in)  :: uref
   real(KIND=r8), intent(out) :: roughnessm

   roughnessm  = (0.0185_r8 / 9.8_r8*(7.59e-4*uref**2+2.46e-2*uref)**2)

end subroutine roughness_0
     
end module roughness_mod
