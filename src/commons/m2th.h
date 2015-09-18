*     -*-fortran-*-

      double precision m2th(4:6)
      character*4 InMasses
*
      double precision Q2th(4:6),m2q(4:6)
      character*4 InMassRef
*
      common / MassThresholdsAPFEL / m2th,InMasses
      common / MassScaleReferAPFEL / Q2th,m2q,InMassRef
