module shr_kind_mod

implicit none

   integer, parameter :: shr_kind_r4 = selected_real_kind(6)  ! 4 byte real
   integer, parameter :: shr_kind_r8 = selected_real_kind(12) ! 8 byte real

   integer, parameter :: shr_kind_i8 = selected_int_kind(18)  ! 8 byte integer

end module shr_kind_mod 

