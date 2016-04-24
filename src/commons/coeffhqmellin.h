*     -*-fortran-*-
*
*     Parameters of the massive coefficient functions tabulation
*
      integer nxi
      double precision xigrid,adler_coef
*
      parameter (nxi=337)
      common / forminimaxAPFEL / xigrid(nxi),adler_coef(nxi)
