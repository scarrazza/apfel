*     -*-fortran-*-
*
*     Parameters of the Q2-grid neede for the caching
*
      integer nQ2g_max
      parameter(nQ2g_max=200)
      double precision Lam2
      parameter(Lam2 = 0.0625d0)
*
      integer nQ2g
      integer inter_degreeQ
      character*4 InQGrid
*
      integer nfin,nffi
      integer nQ(3:7)
      double precision Q2g(-1:nQ2g_max)
*
      common / CacheQgridAPFEL / nQ2g,inter_degreeQ,InQGrid
      common / CacheParametesAPFEL / Q2g,nfin,nffi,nQ
