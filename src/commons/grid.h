*     -*-fortran-*-

      integer nint_max
      parameter(nint_max=160)
      integer ngrid_max
      parameter(ngrid_max=4)
      double precision xmax
      parameter(xmax=1d0)

      integer ngrid,igrid
      integer nin(0:ngrid_max)
      integer inter_degree(0:ngrid_max)
      integer DensityFactor(ngrid_max),TransitionPoint(ngrid_max+1)
      double precision xmin(ngrid_max)
      character*4 InGrid
      logical IsExt(ngrid_max)
      logical ThereAreExtGrids

      double precision xg(0:ngrid_max,0:nint_max)
      double precision step(ngrid_max)

      double precision xgext(0:ngrid_max,0:nint_max)
*
      common / GridParAPFEL / xmin,ngrid,igrid,nin,inter_degree,
     1                        DensityFactor,TransitionPoint,InGrid,
     2                        IsExt,ThereAreExtGrids
      common / GridAPFEL / xg,step
      common / ExtGridAPFEL / xgext
