*     -*-fortran-*-
*
*     Evolution Matrices to Apply to PDFs
*
      integer nfi,nff,sgn
      integer nfli(2:3),nflf(2:3),nli,nlf

      double precision MQCDsg(3:6,2,2,0:nint_max,0:nint_max)
      double precision MQCDnsp(3:6,0:nint_max,0:nint_max)
      double precision MQCDnsm(3:6,0:nint_max,0:nint_max)
      double precision MQCDnsv(3:6,0:nint_max,0:nint_max)

      double precision dMQCDsg(2,2,0:nint_max,0:nint_max)
      double precision dMQCDnsp(0:nint_max,0:nint_max)
      double precision dMQCDnsm(0:nint_max,0:nint_max)
      double precision dMQCDnsv(0:nint_max,0:nint_max)

      double precision MUnisg1(3:6,2:3,5,5,0:nint_max,0:nint_max)
      double precision MUnisg2(3:6,2:3,2,2,0:nint_max,0:nint_max)
      double precision MUninspu(3:6,2:3,0:nint_max,0:nint_max)
      double precision MUninspd(3:6,2:3,0:nint_max,0:nint_max)
      double precision MUninsmu(3:6,2:3,0:nint_max,0:nint_max)
      double precision MUninsmd(3:6,2:3,0:nint_max,0:nint_max)
      double precision MUninslep(2:3,0:nint_max,0:nint_max)
*
      common / ActiveFlavAPFEL / nfi,nff,nli,nlf,nfli,nflf,sgn
      common / EvolMatQCDAPFEL / MQCDsg,MQCDnsp,MQCDnsm,MQCDnsv,
     1                           dMQCDsg,dMQCDnsp,dMQCDnsm,dMQCDnsv
      common / EvolMatUniAPFEL / MUnisg1,MUnisg2,MUninspu,MUninspd,
     1                           MUninsmu,MUninsmd,MUninslep
