*     -*-fortran-*-

      real SC2zmRes(ngrid_max,2,0:nint_max,0:nint_max,0:na)
      real SCLzmRes(ngrid_max,2,0:nint_max,0:nint_max,0:na)
      real SC2charmNCRes(ngrid_max,2,0:nint_max,0:nint_max,0:na)
      real SCLcharmNCRes(ngrid_max,2,0:nint_max,0:nint_max,0:na)
      real SC2charm0NCRes(ngrid_max,2,0:nint_max,0:nint_max,0:na)
      real SCLcharm0NCRes(ngrid_max,2,0:nint_max,0:nint_max,0:na)
*
      common / IntegralsRPResAPFEL / SC2zmRes,SCLzmRes,
     1                               SC2charmNCRes,SCLcharmNCRes,
     2                               SC2charm0NCRes,SCLcharm0NCRes
