*     -*-fortran-*-
*
*     DIS Operators
*
*            SF  grid number | l/c/b/t/tot | PDF index | grid indices
*     real OpF2( 0:jgrid,      3:7,          0:13,       0:nint_max_DIS,0:nint_max_DIS )
*     real OpFL( 0:jgrid,      3:7,          0:13,       0:nint_max_DIS,0:nint_max_DIS )
*     real OpF3( 0:jgrid,      3:7,          0:13,       0:nint_max_DIS,0:nint_max_DIS )
*
      real OpF2(0:ngrid_max,3:7,0:13,0:nint_max_DIS,0:nint_max_DIS)
      real OpFL(0:ngrid_max,3:7,0:13,0:nint_max_DIS,0:nint_max_DIS)
      real OpF3(0:ngrid_max,3:7,0:13,0:nint_max_DIS,0:nint_max_DIS)
*
      common / EvolOpDISAPFEL / OpF2,OpFL,OpF3
