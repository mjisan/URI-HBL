module time_mod

use shr_kind_mod,  only : r4 => shr_kind_r4

  integer,       parameter :: nrecs = 60+1  ! number of time steps
  real(KIND=r4), parameter :: delt  = 3600 ! time step in seconds


end module time_mod

