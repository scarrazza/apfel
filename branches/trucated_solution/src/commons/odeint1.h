*     -*-fortran-*-
*
*     Parameter of the ODE solver
*
      integer maxstp
      double precision tiny
      double precision h1
      double precision eps

      parameter(maxstp=1000)
      parameter(tiny=1d-10)
      parameter(h1=1d-3)
      parameter(eps=1d-3)
