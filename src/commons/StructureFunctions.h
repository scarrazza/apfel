*     -*-fortran-*-
*
*     Structure functions on the grid
*
      double precision F2(3:7,0:ngrid_max,0:nint_max)
      double precision FL(3:7,0:ngrid_max,0:nint_max)
      double precision F3(3:7,0:ngrid_max,0:nint_max)
*
      double precision I2(3:7,0:ngrid_max,0:nint_max)
      double precision I3(3:7,0:ngrid_max,0:nint_max)
*
      common / StructFuncAPFEL / F2,FL,F3
      common / StructFuncTMCAPFEL / I2,I3
