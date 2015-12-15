*     -*-fortran-*-
*
*     Parameters of the Caching
*
*     Number of points of the Q2-grid
      integer nQ2g
      parameter(nQ2g=100)

*     Lambda parameter
      double precision Lam2
      parameter(Lam2 = 0.0625d0)

*     Interpolation degree
      integer inter_degreeQ
      parameter(inter_degreeQ = 3)

      integer nfin,nffi
      integer nQ(3:7)
      double precision Q2g(-1:nQ2g)
*
      double precision fphxQ(-6:6,0:nint_max,nQ2g)
      double precision fgammaxQ(0:nint_max,nQ2g)
      double precision fleptonxQ(-3:3,0:nint_max,nQ2g)
*
      character*4 InCachePDFs
*
      common / CacheParametesAPFEL / Q2g,nfin,nffi,nQ
      common / CachedPDFsAPFEL / fphxQ,fgammaxQ,fleptonxQ,InCachePDFs
