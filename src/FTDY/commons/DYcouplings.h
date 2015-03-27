*     -*-fortran-*-
*
*     Drell-Yan couplings
*
      double precision VV(-6:6,4),AA(-6:6,4)
      double precision CII(-6:6,-6:6,4)
      double precision CFF(-6:6,-6:6,4)
      double precision CIF(-6:6,-6:6,4),CIF_NLO(6,-6:6,4)
*
      common / DrellYanCouplings / VV,AA,CII,CFF,CIF,CIF_NLO
