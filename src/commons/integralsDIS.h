*     -*-fortran-*-

      real SC2zm(0:ngrid_max,3:6,3,0:2,0:nint_max,0:nint_max)
      real SCLzm(0:ngrid_max,3:6,3,0:2,0:nint_max,0:nint_max)
      real SC3zm(0:ngrid_max,3:6,3,0:2,0:nint_max,0:nint_max)
*
      real SC2m(2,0:ngrid_max,nxi,3,0:2,0:nint_max,0:nint_max)
      real SCLm(2,0:ngrid_max,nxi,3,0:2,0:nint_max,0:nint_max)
      real SC3m(2,0:ngrid_max,nxi,3,0:2,0:nint_max,0:nint_max)
*
      real SC2m0(2,0:ngrid_max,nxi,3,0:2,0:nint_max,0:nint_max)
      real SCLm0(2,0:ngrid_max,nxi,3,0:2,0:nint_max,0:nint_max)
      real SC3m0(2,0:ngrid_max,nxi,3,0:2,0:nint_max,0:nint_max)
*
      common / IntegralsDISzmAPFEL / SC2zm,SCLzm,SC3zm
      common / IntegralsDISmAPFEL  / SC2m,SCLm,SC3m
      common / IntegralsDISm0APFEL / SC2m0,SCLm0,SC3m0
