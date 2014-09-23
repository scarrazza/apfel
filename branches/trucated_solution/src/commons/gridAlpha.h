*     -*-fortran-*-

      integer na
      parameter(na=20)
      integer nfg(0:na)
      double precision ag(0:na)
*
      common / gridAlpha / ag,nfg
