************************************************************************
*
*     MSbarmass.f:
*
*     Function that provides the running of the heavy quark masses in
*     the MSbar scheme.
*
************************************************************************
      FUNCTION MSBARMASS(IM,Q2)
*
      IMPLICIT NONE
*
      include "../commons/m2th.h"
      include "../commons/kren.h"
      include "../commons/Evs.h"
      include "../commons/MaxFlavourAlpha.h"
      include "../commons/Nf_FF.h"
      include "../commons/ThresholdAlphaQCD.h"
**
*     Input Variables
*
      INTEGER IM
      DOUBLE PRECISION Q2
**
*     Internal Variables
*
      DOUBLE PRECISION A_QCD,ASQ,ASI,ASIM,ASTH(4:6),ASTHM(4:6)
      DOUBLE PRECISION INMASS,EVF
      DOUBLE PRECISION EVMASS,DECOUP
      DOUBLE PRECISION EPS
      PARAMETER(EPS=1D-10)
**
*     Output Variables
*
      DOUBLE PRECISION MSBARMASS
*
      IF(IM.LT.4.OR.IM.GT.6)THEN
         WRITE(6,*) "In src/Evolution/MSbarmass.f:"
         WRITE(6,*) "IM out of range, IM =",IM
         CALL EXIT(-10)
      ENDIF
*
      ASQ      = A_QCD(Q2)
*
      ASTH(4)  = asthUp(4)
      ASTH(5)  = asthUp(5)
      ASTH(6)  = asthUp(6)
*
      ASTHM(4) = asthDown(4)
      ASTHM(5) = asthDown(5)
      ASTHM(6) = asthDown(6)
*
      ASI      = asthUp(IM)
      ASIM     = asthDown(IM)
*
      INMASS = DSQRT(M2PH(IM))
*
      IF(EVS.EQ."FF")THEN
         EVF = EVMASS(NF_FF,ASI,ASQ)
      ELSEIF(EVS.EQ."VF")THEN
*        Charm
         IF(IM.EQ.4)THEN
            IF(Q2.GE.M2TH(6).AND.NFMAXALPHA.GE.6)THEN
               EVF = EVMASS(4,ASI,ASTHM(5))
     1              * DECOUP("UP",5,DLOG(K2TH(5)))
     2              * EVMASS(5,ASTH(5),ASTHM(6))
     3              * DECOUP("UP",6,DLOG(K2TH(6)))
     4              * EVMASS(6,ASTH(6),ASQ)
            ELSEIF(Q2.GE.M2TH(5).AND.NFMAXALPHA.GE.5)THEN
               EVF = EVMASS(4,ASI,ASTHM(5))
     1              * DECOUP("UP",5,DLOG(K2TH(5)))
     2              * EVMASS(5,ASTH(5),ASQ)
            ELSEIF(Q2.GE.M2TH(4).AND.NFMAXALPHA.GE.4)THEN
               EVF = EVMASS(4,ASI,ASQ)
            ELSE
               EVF = DECOUP("DW",4,DLOG(K2TH(4))) / EVMASS(3,ASQ,ASIM)
            ENDIF
*        Bottom
         ELSEIF(IM.EQ.5)THEN
            IF(Q2.GE.M2TH(6).AND.NFMAXALPHA.GE.6)THEN
               EVF = EVMASS(5,ASI,ASTHM(6))
     1              * DECOUP("UP",6,DLOG(K2TH(6)))
     2              * EVMASS(6,ASTH(6),ASQ)
            ELSEIF(Q2.GE.M2TH(5).AND.NFMAXALPHA.GE.5)THEN
               EVF = EVMASS(5,ASI,ASQ)
            ELSEIF(Q2.GE.M2TH(4).AND.NFMAXALPHA.GE.4)THEN
               EVF = DECOUP("DW",5,DLOG(K2TH(5))) / EVMASS(4,ASQ,ASIM)
            ELSE
               EVF = DECOUP("DW",4,DLOG(K2TH(4)))
     1              / EVMASS(3,ASQ,ASTHM(4))
     2              * DECOUP("DW",5,DLOG(K2TH(5)))
     3              / EVMASS(4,ASTH(4),ASIM)
            ENDIF
*        Top
         ELSEIF(IM.EQ.6)THEN
            IF(Q2.GE.M2TH(6).AND.NFMAXALPHA.GE.6)THEN
               EVF = EVMASS(6,ASI,ASQ)
            ELSEIF(Q2.GE.M2TH(5).AND.NFMAXALPHA.GE.5)THEN
               EVF = DECOUP("DW",6,DLOG(K2TH(6))) / EVMASS(5,ASQ,ASIM)
            ELSEIF(Q2.GE.M2TH(4).AND.NFMAXALPHA.GE.4)THEN
               EVF = DECOUP("DW",5,DLOG(K2TH(5)))
     1              / EVMASS(4,ASQ,ASTHM(5))
     2              * DECOUP("DW",6,DLOG(K2TH(6)))
     3              / EVMASS(5,ASTH(5),ASIM)
            ELSE
               EVF = DECOUP("DW",4,DLOG(K2TH(4)))
     1              / EVMASS(3,ASQ,ASTHM(4))
     2              * DECOUP("DW",5,DLOG(K2TH(5)))
     3              / EVMASS(4,ASTH(4),ASTHM(5))
     4              * DECOUP("DW",6,DLOG(K2TH(6)))
     5              / EVMASS(5,ASTH(5),ASIM)
            ENDIF
         ENDIF
      ENDIF
*
      MSBARMASS = INMASS * EVF
*
      RETURN
      END
*
****************************************************************
      FUNCTION EVMASS(NF,AS0,AS)
*
      IMPLICIT NONE
*
      include "../commons/ipt.h"
      include "../commons/AlphaEvolution.h"
**
*     Input Variables
*
      INTEGER NF
      DOUBLE PRECISION AS,AS0
**
*     Internal Variables
*
      DOUBLE PRECISION FACT,BT0
      DOUBLE PRECISION BETA0APF,BETA1APF,BETA2APF,B1,B2
      DOUBLE PRECISION GAMMA0APF,GAMMA1APF,GAMMA2APF,C0,C1,C2
      DOUBLE PRECISION INTEGRANDMASSRUNNING,DGAUSS,EPS
      PARAMETER(EPS=1D-7)
      EXTERNAL INTEGRANDMASSRUNNING

      INTEGER CMNF,CMIPT
      COMMON / CMASSRUNNING / CMNF,CMIPT
**
*     Output Variables
*
      DOUBLE PRECISION EVMASS
*
      IF(AlphaEvol(1:8).EQ."expanded")THEN
         BT0 = BETA0APF(NF)
*
         C0  = GAMMA0APF() / BT0
         FACT = DEXP( C0 * DLOG( AS / AS0 ) )
         IF(IPT.EQ.0)THEN
            EVMASS = FACT
         ELSEIF(IPT.EQ.1)THEN
            B1 = BETA1APF(NF)  / BT0
            C1 = GAMMA1APF(NF) / BT0
            EVMASS = FACT * ( 1D0 + ( C1 - B1 * C0 ) * AS )
     1                    / ( 1D0 + ( C1 - B1 * C0 ) * AS0 )
         ELSEIF(IPT.EQ.2)THEN
            B1 = BETA1APF(NF)  / BT0
            B2 = BETA2APF(NF)  / BT0
            C1 = GAMMA1APF(NF) / BT0
            C2 = GAMMA2APF(NF) / BT0
            EVMASS = FACT * ( 1D0 + ( C1 - B1 * C0 ) * AS
     1                    + ( C2 - C1 * B1 - B2 * C0
     2                    + B1**2 * C0 + ( C1
     3                    - B1 * C0 )**2 ) * AS**2 / 2D0 )
     4                    / ( 1D0 + ( C1 - B1 * C0 ) * AS0
     5                    + ( C2 - C1 * B1 - B2 * C0
     6                    + B1**2 * C0 + ( C1
     7                    - B1 * C0 )**2 ) * AS0**2 / 2D0 )
         ENDIF
      ELSE
         CMNF  = NF
         CMIPT = IPT
         EVMASS = DEXP(DGAUSS(INTEGRANDMASSRUNNING,AS0,AS,EPS))
      ENDIF
*
      RETURN
      END
*
****************************************************************
      FUNCTION INTEGRANDMASSRUNNING(AS)
**
*     Input Variables
*
      DOUBLE PRECISION AS
**
*     Internal Variables
*
      DOUBLE PRECISION FGAMMA,FBETA

      INTEGER CMNF,CMIPT
      COMMON / CMASSRUNNING / CMNF,CMIPT
**
*     Output Variables
*
      DOUBLE PRECISION INTEGRANDMASSRUNNING
*
      INTEGRANDMASSRUNNING = FGAMMA(AS,CMNF,CMIPT)
     1                     / FBETA(AS,CMNF,CMIPT)
*
      RETURN
      END
*
****************************************************************
      FUNCTION DECOUP(DIR,IM,LN)
*
      IMPLICIT NONE
*
      include "../commons/ipt.h"
      include "../commons/m2th.h"
      include "../commons/kren.h"
      include "../commons/ThresholdAlphaQCD.h"
**
*     Input Variables
*
      INTEGER IM
      DOUBLE PRECISION LN
      CHARACTER*2 DIR
**
*     Internal Variables
*
      DOUBLE PRECISION ASTH2
      DOUBLE PRECISION EPS
      PARAMETER(EPS=1D-7)
**
*     Output Variables
*
      DOUBLE PRECISION DECOUP
*
      IF(IPT.LE.1)THEN
         DECOUP = 1D0
         RETURN
      ELSE
         IF(DIR.EQ."DW")THEN
            ASTH2  = asthUp(IM)**2
            DECOUP = 1D0 + ASTH2 * ( 89D0 / 27D0 - 20D0 / 9D0 * LN
     1                             + 4D0 / 3D0 * LN**2 )
         ELSEIF(DIR.EQ."UP")THEN
            ASTH2 = asthDown(IM)**2
            DECOUP = 1D0 - ASTH2 * ( 89D0 / 27D0 - 20D0 / 9D0 * LN
     1                             + 4D0 / 3D0 * LN**2 )
         ELSE
            WRITE(6,*) "In src/Evolution/MSbarmass.f:"
            WRITE(6,*) "Unknown direction, DIR =",DIR
            CALL EXIT(-10)
         ENDIF
      ENDIF
*
      RETURN
      END
*
****************************************************************************
*
*     QCD gamma function.
*
****************************************************************************
      function fgamma(a,nf,ipt)
*
      implicit none
**
*     Input Variables
*
      double precision a
      integer nf,ipt
**
*     Internal Variables
*
      double precision gamma0apf,gamma1apf,gamma2apf,gamma3apf
**
*     Output Variables
*
      double precision fgamma
*
      if(ipt.eq.0)then
         fgamma = - a * gamma0apf()
      elseif(ipt.eq.1)then
         fgamma = - a * ( gamma0apf() + a * gamma1apf(nf) )
      elseif(ipt.eq.2)then
         fgamma = - a * ( gamma0apf()
     1        + a * ( gamma1apf(nf) + a * gamma2apf(nf) ) )
      elseif(ipt.eq.3)then
         fgamma = - a * ( gamma0apf()
     1        + a * ( gamma1apf(nf)
     2        + a * ( gamma2apf(nf) + a * gamma3apf(nf) ) ) )
      endif
*
      return
      end
*
****************************************************************************
      function gamma0apf()
*
      implicit none
**
*     Output Variables
*
      double precision gamma0apf
*
      gamma0apf = 4d0
*
      return
      end
*
****************************************************************************
      function gamma1apf(nf)
*
      implicit none
**
*     Input Variables
*
      integer nf
**
*     Output Variables
*
      double precision gamma1apf
*
      gamma1apf = 202d0 / 3d0 - 20d0 / 9d0 * nf
*
      return
      end
*
****************************************************************************
      function gamma2apf(nf)
*
      implicit none
*
      include "../commons/consts.h"
**
*     Input Variables
*
      integer nf
**
*     Output Variables
*
      double precision gamma2apf
*
      gamma2apf = 1249d0 - ( 2216d0 / 27d0 + 160d0 / 3d0 * zeta3 ) * nf 
     1          - 140d0 / 81d0 * nf**2
*
      return
      end
*
****************************************************************************
      function gamma3apf(nf)
*
      implicit none
*
      include "../commons/consts.h"
**
*     Input Variables
*
      integer nf
**
*     Output Variables
*
      double precision gamma3apf
*
      gamma3apf = 4603055d0 / 162d0 + 135680d0 * zeta3 / 28d0
     1     - 8800d0 * zeta5
     2     + ( - 91723d0 / 27d0 - 34192d0 * zeta3 / 9d0 + 880d0 * zeta4
     3     + 18400d0 * zeta5 / 9d0 ) * nf
     5     + ( 5242d0 / 243d0 + 800d0 * zeta3 / 9d0
     6     - 160d0 * zeta4 / 3d0 ) * nf**2
     7     + ( 332d0 / 243d0 + 64d0 * zeta3 / 27d0 ) * nf**3
*
      return
      end
*
****************************************************************************
*
*     Function whose zero is the so-called RG invariant mass.
*     For the MSbar mass, it finds the scale m such that m(m) = m.
*
****************************************************************************
      function MassQSplit(i,Q)
*
      implicit none
*
      include "../commons/m2th.h"
      include "../commons/kren.h"
      include "../commons/Evs.h"
      include "../commons/MaxFlavourAlpha.h"
      include "../commons/Nf_FF.h"
      include "../commons/ThresholdAlphaQCD.h"
**
*     Input Variables
*
      integer i
      double precision Q
**
*     Internal Variables
*
      DOUBLE PRECISION Q2
      DOUBLE PRECISION A_QCD,ASQ,ASI,ASTH(4:6),ASTHM(4:6)
      DOUBLE PRECISION LN
      DOUBLE PRECISION INMASS,EVF
      DOUBLE PRECISION EVMASS,DECOUP
*
*     Output Variables
*
      double precision MassQSplit
*
      Q2  = Q * Q
*
      ASI = A_QCD(Q2)
      ASQ = A_QCD(Q2TH(I))
*
      ASTH(4)  = asthUp(4)
      ASTH(5)  = asthUp(5)
      ASTH(6)  = asthUp(6)
*
      ASTHM(4) = asthDown(4)
      ASTHM(5) = asthDown(5)
      ASTHM(6) = asthDown(6)
*
      INMASS = DSQRT(M2Q(I))
*
      IF(EVS.EQ."FF")THEN
         EVF = EVMASS(NF_FF,ASI,ASQ)
      ELSEIF(EVS.EQ."VF")THEN
         LN = DLOG(KREN)
*        Charm
         IF(I.EQ.4)THEN
            IF(Q2.GE.M2TH(6).AND.NFMAXALPHA.GE.6)THEN
               EVF = EVMASS(6,ASI,ASQ)
            ELSEIF(Q2.GE.M2TH(5).AND.NFMAXALPHA.GE.5)THEN
               IF(Q2TH(4).GE.M2TH(6).AND.NFMAXALPHA.GE.6)THEN
                  EVF = EVMASS(5,ASI,    ASTHM(6)) * DECOUP("UP",6,LN)
     1                * EVMASS(6,ASTH(6),ASQ)
               ELSE
                  EVF = EVMASS(5,ASI,ASQ)
               ENDIF
            ELSEIF(Q2.GE.M2TH(4).AND.NFMAXALPHA.GE.4)THEN
               IF(Q2TH(4).GE.M2TH(6).AND.NFMAXALPHA.GE.6)THEN
                  EVF = EVMASS(4,ASI,     ASTHM(5)) * DECOUP("UP",5,LN)
     1                * EVMASS(5,ASTH(5), ASTHM(6)) * DECOUP("UP",6,LN)
     2                * EVMASS(6,ASTH(6),ASQ)
               ELSEIF(Q2TH(4).GE.M2TH(5).AND.NFMAXALPHA.GE.5)THEN
                  EVF = EVMASS(4,ASI,    ASTHM(5)) * DECOUP("UP",5,LN)
     2                * EVMASS(5,ASTH(5),ASQ)
               ELSE
                  EVF = EVMASS(4,ASI,ASQ)
               ENDIF
            ENDIF
*        Bottom
         ELSEIF(I.EQ.5)THEN
            IF(Q2.GE.M2TH(6).AND.NFMAXALPHA.GE.6)THEN
               EVF = EVMASS(6,ASI,ASQ)
            ELSEIF(Q2.GE.M2TH(5).AND.NFMAXALPHA.GE.5)THEN
               IF(Q2TH(5).GE.M2TH(6).AND.NFMAXALPHA.GE.6)THEN
                  EVF = EVMASS(5,ASI,    ASTHM(6)) * DECOUP("UP",6,LN)
     1                * EVMASS(6,ASTH(6),ASQ)
               ELSE
                  EVF = EVMASS(5,ASI,ASQ)
               ENDIF
            ENDIF
*        Top
         ELSEIF(I.EQ.6)THEN
            EVF = EVMASS(6,ASI,ASQ) 
         ENDIF
      ENDIF
*
      MassQSplit = INMASS / EVF - Q
*
      return
      end
*
****************************************************************************
*
*     Subroutine that computes the RG invariant masses
*
****************************************************************************
      subroutine ComputeRGInvariantMasses
*
      implicit none
*
      include "../commons/m2th.h"
**
*     Internal Variables
*
      integer i
      double precision MassQSplit,zriddr
      double precision x1,x2
      double precision acc
      parameter(acc=1d-10)
      external MassQSplit
*
      do i=6,4,-1
         call ThresholdAlphaQCD
         if(q2th(i).ne.m2q(i))then
            x1 = dsqrt(m2q(i))
            x2 = dsqrt(q2th(i))
            m2ph(i) = zriddr(MassQSplit,i,x1,x2,acc)**2
         endif
         call ComputeHeavyQuarkThresholds
      enddo
*
      return
      end
