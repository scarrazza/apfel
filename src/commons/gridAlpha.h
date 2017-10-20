*     -*-fortran-*-

      integer na
      parameter(na=20)
      integer nfg(0:na)
      double precision ag(0:na),qag(0:na)
      double precision asmh(4:6)
*
      common / gridAlphaAPFEL / ag,qag,asmh,nfg
