*     -*-fortran-*-

      real SC2zm(0:ngrid_max,3:6,4,0:2,0:nint_max_DIS,0:nint_max_DIS)
      real SCLzm(0:ngrid_max,3:6,4,0:2,0:nint_max_DIS,0:nint_max_DIS)
      real SC3zm(0:ngrid_max,3:6,4,0:2,0:nint_max_DIS,0:nint_max_DIS)
*
      real SC2mNC(0:ngrid_max,nxi,3,0:2,0:nint_max_DIS,0:nint_max_DIS)
      real SCLmNC(0:ngrid_max,nxi,3,0:2,0:nint_max_DIS,0:nint_max_DIS)
      real SC3mNC(0:ngrid_max,nxi,3,0:2,0:nint_max_DIS,0:nint_max_DIS)
*
      real SC2mCC(0:ngrid_max,nxi,3,0:2,0:nint_max_DIS,0:nint_max_DIS)
      real SCLmCC(0:ngrid_max,nxi,3,0:2,0:nint_max_DIS,0:nint_max_DIS)
      real SC3mCC(0:ngrid_max,nxi,3,0:2,0:nint_max_DIS,0:nint_max_DIS)
*
      real SC2m0NC(0:ngrid_max,nxi,3,0:2,0:nint_max_DIS,0:nint_max_DIS)
      real SCLm0NC(0:ngrid_max,nxi,3,0:2,0:nint_max_DIS,0:nint_max_DIS)
      real SC3m0NC(0:ngrid_max,nxi,3,0:2,0:nint_max_DIS,0:nint_max_DIS)
*
      real SC2m0CC(0:ngrid_max,nxi,3,0:2,0:nint_max_DIS,0:nint_max_DIS)
      real SCLm0CC(0:ngrid_max,nxi,3,0:2,0:nint_max_DIS,0:nint_max_DIS)
      real SC3m0CC(0:ngrid_max,nxi,3,0:2,0:nint_max_DIS,0:nint_max_DIS)
*
      common / IntegralsDISzmAPFEL / SC2zm,SCLzm,SC3zm
      common / IntegralsDISmAPFEL  / SC2mNC,SCLmNC,SC3mNC,
     1                               SC2mCC,SCLmCC,SC3mCC
      common / IntegralsDISm0APFEL / SC2m0NC,SCLm0NC,SC3m0NC,
     1                               SC2m0CC,SCLm0CC,SC3m0CC
