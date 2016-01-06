*     -*-fortran-*-
*
*     Structure functions on the (x,Q2)-grid
*
      double precision SFxQ(0:2,3,3:7,0:nint_max,0:nQ2g_max)
      double precision ISFxQ(0:2,2,3:7,0:nint_max,0:nQ2g_max)
*
      character*4 InCacheSFs
*
      common / CachedSFsAPFEL / SFxQ,ISFxQ,InCacheSFs
