module time_mod

use shr_kind_mod,  only : r8 => shr_kind_r8

  real(KIND=r8) :: delt  = 1.0                         ! model time step (sec)
  
  integer, parameter :: spinup_nsteps     =  0*3600    ! number of spinup steps
  integer, parameter :: forecast_nsteps   =  10*3600    ! number of forecast_step

  integer, parameter :: spinup_nrecs      =  0 + 1    ! number of records in spinup
  integer, parameter :: forecast_nrecs    =  10*3600/3600       ! number of records in forecast

  real(KIND=r8) :: spinup_recs_int        =    3600    ! interval for spinup recs (secs) 
  real(KIND=r8) :: forecast_recs_int      =    3600    ! interval for forecast recs (secs) 

  integer, parameter :: nsteps = spinup_nsteps + forecast_nsteps  ! total number of model steps
  integer, parameter :: nrecs  = spinup_nrecs  + forecast_nrecs   ! total number of records in model

end module time_mod

