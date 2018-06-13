*     -*-fortran-*-

      double precision m2ph(4:6)
      character*4 InMasses
*
      double precision Q2th(4:6),m2q(4:6)
      character*4 InMassRef
*
      double precision k2th(4:6)
      character*4 InThrRatios
*
      double precision m2th(4:8)
*
      common / MassPhysicalAPFEL / m2ph,InMasses
      common / MassScaleReferAPFEL / Q2th,m2q,InMassRef
      common / MassThreRatiosAPFEL / k2th,InThrRatios
      common / MassThresholdsAPFEL / m2th
