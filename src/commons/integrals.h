*     -*-fortran-*-

      real SP(0:ngrid_max,3:6,7,0:2,0:nint_max,0:nint_max)
      real SQ(ngrid_max,3:6,2:3,25,0:2,0:nint_max,0:nint_max)
      real SM(ngrid_max,4:6,5,0:2,0:nint_max,0:nint_max)
*
      common / IntegralsRPAPFEL / SP
      common / IntegralsRQAPFEL / SQ
      common / IntegralsRMAPFEL / SM
