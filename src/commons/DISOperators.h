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
      real OpI2(0:ngrid_max,3:7,0:13,0:nint_max_DIS,0:nint_max_DIS)
      real OpI3(0:ngrid_max,3:7,0:13,0:nint_max_DIS,0:nint_max_DIS)

      real EvOpF2(3:7,0:13,0:nint_max_DIS,0:nint_max_DIS)
      real EvOpFL(3:7,0:13,0:nint_max_DIS,0:nint_max_DIS)
      real EvOpF3(3:7,0:13,0:nint_max_DIS,0:nint_max_DIS)
      real EvOpI2(3:7,0:13,0:nint_max_DIS,0:nint_max_DIS)
      real EvOpI3(3:7,0:13,0:nint_max_DIS,0:nint_max_DIS)
*
      common / OperatorDISAPFEL / OpF2,OpFL,OpF3
      common / OperatorDISTMCAPFEL / OpI2,OpI3
      common / EvOperatorDISAPFEL / EvOpF2,EvOpFL,EvOpF3
      common / EvOperatorDISTMCAPFEL / EvOpI2,EvOpI3
