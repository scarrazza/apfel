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
      double precision xmin(ngrid_max)
      character*4 InGrid

      double precision xg(0:ngrid_max,0:nint_max)
      double precision step(ngrid_max)
*
      common / GridPar / ngrid,igrid,nin,inter_degree,xmin,InGrid
      common / grid / xg,step
