program model_main

use shr_kind_mod,   only : r8 => shr_kind_r8
use mpp_setup_mod,  only : mpp_init, mpp_end, mpp_barrier
use mpp_setup_mod,  only : size, master, rank
use mpp_io_mod,     only : stdout
use tprof_mod,      only : tstart, tstop, tnull, tprnt
use model_step_mod, only : model_step
use time_mod,       only : nsteps, delt
use grid_mod,       only : nz, npx, npy

implicit none

  integer :: step
  real(KIND=r8) :: time_sec, time_hrs

  call mpp_init ( npx, npy )

! start time loop

  call tstart('MODEL_MAIN_LOOP') 

! call boundary model

  do step=1,nsteps
    time_sec = (step-1)*delt
    time_hrs  = time_sec/3600.0
 
    if (mod(step,60) == 0) then
      if ( rank == master) then
        write(stdout, 1000) step, time_hrs
1000 format(2x, 'step = ', i8, 5x, 'time in hours = ', f6.3)
      endif
    endif

    call model_step (time_sec) 
!print *, 'model_step called: ', time_sec

  enddo

  call tstop ('MODEL_MAIN_LOOP')

  if ( rank == master ) then
    call tprnt
  endif

  call mpp_end

end program model_main

