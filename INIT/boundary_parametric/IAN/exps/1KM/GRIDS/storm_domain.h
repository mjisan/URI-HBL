integer, parameter :: ifplane = 0

  integer, parameter :: im = 2200
  integer, parameter :: jm = 1800

  real(KIND=r4), parameter :: lonstart = -120.0
  real(KIND=r4), parameter :: lonend   = -28.375 ! lonend = lonstart + (im-1)*resol

  real(KIND=r4), parameter :: latstart =  0.0
  real(KIND=r4), parameter :: latend   =  74.9583 ! latend = latstart + (jm-1)*resol

  real(KIND=r4), parameter :: resol  =  1.0/24.0
