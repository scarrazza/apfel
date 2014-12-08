*     -*-fortran-*-

      integer na
      parameter(na=20)
      integer nfg(0:na)
      double precision ag(0:na)
*
      common / gridAlphaAPFEL / ag,nfg
