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
**
*     Input Variables
*
      INTEGER IM
      DOUBLE PRECISION Q2
**
*     Internal Variables
*
      DOUBLE PRECISION MU2
      DOUBLE PRECISION A_QCD,ASQ,ASTH(4:6),ASTHM(4:6)
      DOUBLE PRECISION LN
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
      ASTH(4)  = A_QCD(M2TH(4))
      ASTH(5)  = A_QCD(M2TH(5))
      ASTH(6)  = A_QCD(M2TH(6))
*
      ASTHM(4) = A_QCD(M2TH(4)-EPS)
      ASTHM(5) = A_QCD(M2TH(5)-EPS)
      ASTHM(6) = A_QCD(M2TH(6)-EPS)
*
      INMASS = DSQRT(M2TH(IM))
*
      MU2 = KREN * Q2
*
      IF(EVS.EQ."FF")THEN
         EVF = EVMASS(NF_FF,ASTH(IM),ASQ)
      ELSEIF(EVS.EQ."VF")THEN
         LN = DLOG(KREN)
*        Charm
         IF(IM.EQ.4)THEN
            IF(MU2.GE.M2TH(6).AND.NFMAXALPHA.GE.6)THEN
               EVF = EVMASS(4,ASTH(4),ASTHM(5)) * DECOUP("UP",5,LN)
     1             * EVMASS(5,ASTH(5),ASTHM(6)) * DECOUP("UP",6,LN)
     2             * EVMASS(6,ASTH(6),ASQ)
            ELSEIF(MU2.GE.M2TH(5).AND.NFMAXALPHA.GE.5)THEN
               EVF = EVMASS(4,ASTH(4),ASTHM(5)) * DECOUP("UP",5,LN)
     1             * EVMASS(5,ASTH(5),ASQ)
            ELSEIF(MU2.GE.M2TH(4).AND.NFMAXALPHA.GE.4)THEN
               EVF = EVMASS(4,ASTH(4),ASQ)
            ELSE
               EVF = DECOUP("DW",4,LN) / EVMASS(3,ASQ,ASTHM(4)) 
            ENDIF
*        Bottom
         ELSEIF(IM.EQ.5)THEN
            IF(MU2.GE.M2TH(6).AND.NFMAXALPHA.GE.6)THEN
               EVF = EVMASS(5,ASTH(5),ASTHM(6)) * DECOUP("UP",6,LN)
     1             * EVMASS(6,ASTH(6),ASQ)
            ELSEIF(MU2.GE.M2TH(5).AND.NFMAXALPHA.GE.5)THEN
               EVF = EVMASS(5,ASTH(5),ASQ) 
            ELSEIF(MU2.GE.M2TH(4).AND.NFMAXALPHA.GE.4)THEN
               EVF = DECOUP("DW",5,LN) / EVMASS(4,ASQ,ASTHM(5)) 
            ELSE
               EVF = DECOUP("DW",4,LN) / EVMASS(3,ASQ,ASTHM(4)) 
     1             * DECOUP("DW",5,LN) / EVMASS(4,ASTH(4),ASTHM(5)) 
            ENDIF
*        Top
         ELSEIF(IM.EQ.6)THEN
            IF(MU2.GE.M2TH(6).AND.NFMAXALPHA.GE.6)THEN
               EVF = EVMASS(6,ASTH(6),ASQ) 
            ELSEIF(MU2.GE.M2TH(5).AND.NFMAXALPHA.GE.5)THEN
               EVF = DECOUP("DW",6,LN) / EVMASS(5,ASQ,ASTHM(6)) 
            ELSEIF(MU2.GE.M2TH(4).AND.NFMAXALPHA.GE.4)THEN
               EVF = DECOUP("DW",5,LN) / EVMASS(4,ASQ,ASTHM(5)) 
     1             * DECOUP("DW",6,LN) / EVMASS(5,ASTH(5),ASTHM(6)) 
            ELSE
               EVF = DECOUP("DW",4,LN) / EVMASS(3,ASQ,ASTHM(4))
     1             * DECOUP("DW",5,LN) / EVMASS(4,ASTH(4),ASTHM(5)) 
     2             * DECOUP("DW",6,LN) / EVMASS(5,ASTH(5),ASTHM(6)) 
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
**
*     Input Variables
*
      INTEGER NF
      DOUBLE PRECISION BETA0APF,BETA1APF,BETA2APF,B1,B2
      DOUBLE PRECISION GAMMA0APF,GAMMA1APF,GAMMA2APF,C0,C1,C2
      DOUBLE PRECISION AS,AS0
**
*     Internal Variables
*
      DOUBLE PRECISION FACT
**
*     Output Variables
*
      DOUBLE PRECISION EVMASS
*
      B1 = BETA1APF(NF)  / BETA0APF(NF)
      B2 = BETA2APF(NF)  / BETA0APF(NF)
      C0 = GAMMA0APF()   / BETA0APF(NF)
      C1 = GAMMA1APF(NF) / BETA0APF(NF)
      C2 = GAMMA2APF(NF) / BETA0APF(NF)
*
      FACT = DEXP( C0 * DLOG( AS / AS0 ) )
      IF(IPT.EQ.0)THEN
         EVMASS = FACT
      ELSEIF(IPT.EQ.1)THEN
         EVMASS = FACT * ( 1D0 + ( C1 - B1 * C0 ) * AS )
     1                 / ( 1D0 + ( C1 - B1 * C0 ) * AS0 )
      ELSEIF(IPT.EQ.2)THEN
         EVMASS = FACT * ( 1D0 + ( C1 - B1 * C0 ) * AS 
     1                 + ( C2 - C1 * B1 - B2 * C0 
     2                 + B1**2D0 * C0 + ( C1 
     3                 - B1 * C0 )**2D0 ) * AS**2D0 / 2D0 )
     4                 / ( 1D0 + ( C1 - B1 * C0 ) * AS0 
     5                 + ( C2 - C1 * B1 - B2 * C0 
     6                 + B1**2D0 * C0 + ( C1 
     7                 - B1 * C0 )**2D0 ) * AS0**2D0 / 2D0 )
      ENDIF
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
**
*     Input Variables
*
      INTEGER IM
      DOUBLE PRECISION LN
      CHARACTER*2 DIR
**
*     Internal Variables
*
      DOUBLE PRECISION A_QCD,ASTH2(4:6),ASTHM2(4:6)
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
         ASTH2(4)  = A_QCD(M2TH(4))**2D0
         ASTH2(5)  = A_QCD(M2TH(5))**2D0
         ASTH2(6)  = A_QCD(M2TH(6))**2D0
*
         ASTHM2(4) = A_QCD(M2TH(4)-EPS)**2D0
         ASTHM2(5) = A_QCD(M2TH(5)-EPS)**2D0
         ASTHM2(6) = A_QCD(M2TH(6)-EPS)**2D0
*
         IF(DIR.EQ."DW")THEN
            DECOUP = 1D0 + ASTH2(IM) * ( 89D0 / 27D0 - 20D0 / 9D0 * LN 
     1                                 + 4D0 / 3D0 * LN**2D0 )
         ELSEIF(DIR.EQ."UP")THEN
            DECOUP = 1D0 - ASTHM2(IM) * ( 89D0 / 27D0 - 20D0 / 9D0 * LN 
     1                                  + 4D0 / 3D0 * LN**2D0 )
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
      double precision gamma0apf,gamma1apf,gamma2apf
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
     1            + a * ( gamma1apf(nf) + a * gamma2apf(nf) ) )
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
     1          - 140d0 / 81d0 * nf**2d0
*
      return
      end
