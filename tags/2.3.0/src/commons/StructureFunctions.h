*     -*-fortran-*-
*
*     Structure functions on the grid
*
      double precision Q2DIS
      double precision F2(1:7,0:ngrid_max,0:nint_max)
      double precision FL(1:7,0:ngrid_max,0:nint_max)
      double precision F3(1:7,0:ngrid_max,0:nint_max)
*
      common / ZMStructFuncAPFEL / Q2DIS,F2,FL,F3
