*     -*-fortran-*-
*
*     Parameter of the MINIMAX tabulations
*
      integer nxi,n_pol,ncoef,m_coef
      double precision ximin,ximax,xigrid,adler_coef
      double precision coef_p1,coef_p2,coef
*
      parameter (nxi=300,n_pol=15,ncoef=14)
      common /forminimaxAPFEL/ ximin,ximax,xigrid(nxi),adler_coef(nxi),
     1                         coef_p1(ncoef),coef_p2(ncoef),
     2                         coef(ncoef,nxi,n_pol),m_coef(ncoef)
