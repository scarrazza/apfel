*     -*-fortran-*-

      double precision SP(ngrid_max,3:6,7,0:2,0:nint_max,0:nint_max)
      double precision SQ(ngrid_max,3:6,10,0:nint_max,0:nint_max)
      double precision SM(ngrid_max,5,0:2,0:nint_max,0:nint_max)
*
      common / IntegralsRP / SP
      common / IntegralsRQ / SQ
      common / IntegralsRM / SM
