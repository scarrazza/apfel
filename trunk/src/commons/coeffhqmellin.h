*     -*-fortran-*-
*
*     Parameter of the MINIMAX tabulations
*
      integer nxi,n_pol,ncoef,m_coef
      double precision ximin,ximax,xigrid,adler_coef
      double precision coef_p1,coef_p2,coef
*
      parameter (nxi=300,n_pol=15,ncoef=14)
      common / forminimaxAPFEL / ximin,ximax,xigrid(nxi),
     1                           adler_coef(nxi),coef_p1(ncoef),
     2                           coef_p2(ncoef),coef(ncoef,nxi,n_pol),
     3                           m_coef(ncoef)
