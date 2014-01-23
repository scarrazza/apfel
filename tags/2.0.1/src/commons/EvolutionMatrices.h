*     -*-fortran-*-
*
*     Evolution Matrices to Apply to PDFs
*
      integer nfi,nff,sgn

      double precision MQCDsg(3:6,2,2,0:nint_max,0:nint_max)
      double precision MQCDnsp(3:6,0:nint_max,0:nint_max)
      double precision MQCDnsm(3:6,0:nint_max,0:nint_max)
      double precision MQCDnsv(3:6,0:nint_max,0:nint_max)

      double precision MQEDsg(3:6,3,3,0:nint_max,0:nint_max)
      double precision MQEDnsp(3:6,0:nint_max,0:nint_max)
      double precision MQEDnsm(3:6,0:nint_max,0:nint_max)
*
      common / ActiveFlav / nfi,nff,sgn
      common / EvolMatQCD / MQCDsg,MQCDnsp,MQCDnsm,MQCDnsv
      common / EvolMatQED / MQEDsg,MQEDnsp,MQEDnsm
