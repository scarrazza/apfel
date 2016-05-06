*     -*-fortran-*-
*
*     Variables needed in the computation of the O(as) contribution
*     to the intrinsic charm components
*
      double precision m2strange
      parameter(m2strange = 0.01d0)
      double precision one
      parameter(one=0.999999999999999d0)
*
      double precision Q2IC,m12,m22,m1,m2
      double precision Splus,Sminus,Rplus,Rminus
      double precision Del,Del2,Spp,Spm,Smp,eta
      double precision I1,Cplus,C1m,C1p,CRm
      double precision S1,S2,S3,V1,V2,V3
      double precision fact1,fact2,fact3,factL
*
      common / ICwrapScales / Q2IC,m12,m22,m1,m2
      common / ICwrapCouplings / Splus,Sminus,Rplus,Rminus
      common / ICwrapKinVar / Del,Del2,Spp,Spm,Smp,eta
      common / ICwrapVertex / I1,Cplus,C1m,C1p,CRm
      common / ICwrapRealVirt / S1,S2,S3,V1,V2,V3
      common / ICwrapKinFact / fact1,fact2,fact3,factL
