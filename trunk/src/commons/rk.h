*     -*-fortran-*-
*
*     Parameters of the Runge-Kutta method
*
      integer nmax
      integer maxstp
      double precision tiny
      double precision safety,pgrow,pshrnk,errcon
*
      parameter(nmax=10000)
      parameter(maxstp=1000)
      parameter(tiny=1d-30)
      parameter(safety=0.9d0,pgrow=-0.2d0,pshrnk=-0.25d0,errcon=1.89d-4)
