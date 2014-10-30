*     -*-fortran-*-

      double precision SC2zm(0:ngrid_max,3:6,3,0:2,
     1                       0:nint_max,0:nint_max)
      double precision SCLzm(0:ngrid_max,3:6,3,0:2,
     1                       0:nint_max,0:nint_max)
      double precision SC3zm(0:ngrid_max,3:6,3,0:2,
     1                       0:nint_max,0:nint_max)
*
      double precision SC2m(2,0:ngrid_max,nxi,3,0:2,
     1                      0:nint_max,0:nint_max)
      double precision SCLm(2,0:ngrid_max,nxi,3,0:2,
     1                      0:nint_max,0:nint_max)
      double precision SC3m(2,0:ngrid_max,nxi,3,0:2,
     1                      0:nint_max,0:nint_max)
*
      double precision SC2m0(2,0:ngrid_max,nxi,3,0:2,
     1                       0:nint_max,0:nint_max)
      double precision SCLm0(2,0:ngrid_max,nxi,3,0:2,
     1                       0:nint_max,0:nint_max)
      double precision SC3m0(2,0:ngrid_max,nxi,3,0:2,
     1                       0:nint_max,0:nint_max)
*
      common / IntegralsDISzmAPFEL / SC2zm,SCLzm,SC3zm
      common / IntegralsDISmAPFEL  / SC2m,SCLm,SC3m
      common / IntegralsDISm0APFEL / SC2m0,SCLm0,SC3m0
