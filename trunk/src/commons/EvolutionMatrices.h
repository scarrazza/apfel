*     -*-fortran-*-
*
*     Evolution Matrices to Apply to PDFs
*
      integer nfi,nff,sgn

      real MQCDsg(3:6,2,2,0:nint_max,0:nint_max)
      real MQCDnsp(3:6,0:nint_max,0:nint_max)
      real MQCDnsm(3:6,0:nint_max,0:nint_max)
      real MQCDnsv(3:6,0:nint_max,0:nint_max)

      real dMQCDsg(2,2,0:nint_max,0:nint_max)
      real dMQCDnsp(0:nint_max,0:nint_max)
      real dMQCDnsm(0:nint_max,0:nint_max)
      real dMQCDnsv(0:nint_max,0:nint_max)

      real MQEDsg(3:6,3,3,0:nint_max,0:nint_max)
      real MQEDnsp(3:6,0:nint_max,0:nint_max)
      real MQEDnsm(3:6,0:nint_max,0:nint_max)

      real dMQEDsg(3,3,0:nint_max,0:nint_max)
      real dMQEDnsp(0:nint_max,0:nint_max)
      real dMQEDnsm(0:nint_max,0:nint_max)

      real MUnisg1(3:6,4,4,0:nint_max,0:nint_max)
      real MUnisg2(3:6,2,2,0:nint_max,0:nint_max)
      real MUninspu(3:6,0:nint_max,0:nint_max)
      real MUninspd(3:6,0:nint_max,0:nint_max)
      real MUninsmu(3:6,0:nint_max,0:nint_max)
      real MUninsmd(3:6,0:nint_max,0:nint_max)
*
      common / ActiveFlavAPFEL / nfi,nff,sgn
      common / EvolMatQCDAPFEL / MQCDsg,MQCDnsp,MQCDnsm,MQCDnsv,
     1                           dMQCDsg,dMQCDnsp,dMQCDnsm,dMQCDnsv
      common / EvolMatQEDAPFEL / MQEDsg,MQEDnsp,MQEDnsm,
     1                           dMQEDsg,dMQEDnsp,dMQEDnsm
      common / EvolMatUniAPFEL / MUnisg1,MUnisg2,MUninspu,MUninspd,
     1                           MUninsmu,MUninsmd
