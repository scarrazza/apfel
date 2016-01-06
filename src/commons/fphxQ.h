*     -*-fortran-*-
*
*     PDFs on the (x,Q2)-grid
*
      double precision fphxQ(-6:6,0:nint_max,0:nQ2g_max)
      double precision fgammaxQ(0:nint_max,0:nQ2g_max)
      double precision fleptonxQ(-3:3,0:nint_max,0:nQ2g_max)
*
      character*4 InCachePDFs
*
      common / CachedPDFsAPFEL / fphxQ,fgammaxQ,fleptonxQ,InCachePDFs
