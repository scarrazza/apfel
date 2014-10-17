*     -*-fortran-*-

      double precision SC2zm(0:ngrid_max,3:6,3,0:2,
     1                       0:nint_max,0:nint_max)
      double precision SCLzm(0:ngrid_max,3:6,3,0:2,
     1                       0:nint_max,0:nint_max)
      double precision SC3zm(0:ngrid_max,3:6,3,0:2,
     1                       0:nint_max,0:nint_max)
*
      double precision SC2m(0:ngrid_max,nxi,3,0:2,
     1                      0:nint_max,0:nint_max)
      double precision SCLm(0:ngrid_max,nxi,3,0:2,
     1                      0:nint_max,0:nint_max)
      double precision SC3m(0:ngrid_max,nxi,3,0:2,
     1                      0:nint_max,0:nint_max)
*
      common / IntegralsDISzmAPFEL / SC2zm,SCLzm,SC3zm
      common / IntegralsDISmAPFEL  / SC2m,SCLm,SC3m
