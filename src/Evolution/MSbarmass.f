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
      DOUBLE PRECISION A_QCD,ASQ,ASI,ASIM,ASTH(4:6),ASTHM(4:6)
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
      ASI      = A_QCD(M2Q(IM))
      ASIM     = A_QCD(M2Q(IM)-EPS)
*
      INMASS = DSQRT(M2Q(IM))
*
      MU2 = KREN * Q2
*
      IF(EVS.EQ."FF")THEN
         EVF = EVMASS(NF_FF,ASI,ASQ)
      ELSEIF(EVS.EQ."VF")THEN
         LN = DLOG(KREN)
*        Charm
         IF(IM.EQ.4)THEN
            IF(MU2.GE.M2TH(6).AND.NFMAXALPHA.GE.6)THEN
               EVF = EVMASS(4,ASI,ASTHM(5)) * DECOUP("UP",5,LN)
     1             * EVMASS(5,ASTH(5),ASTHM(6)) * DECOUP("UP",6,LN)
     2             * EVMASS(6,ASTH(6),ASQ)
            ELSEIF(MU2.GE.M2TH(5).AND.NFMAXALPHA.GE.5)THEN
               EVF = EVMASS(4,ASI,ASTHM(5)) * DECOUP("UP",5,LN)
     1             * EVMASS(5,ASTH(5),ASQ)
            ELSEIF(MU2.GE.M2TH(4).AND.NFMAXALPHA.GE.4)THEN
               EVF = EVMASS(4,ASI,ASQ)
            ELSE
               EVF = DECOUP("DW",4,LN) / EVMASS(3,ASQ,ASIM) 
            ENDIF
*        Bottom
         ELSEIF(IM.EQ.5)THEN
            IF(MU2.GE.M2TH(6).AND.NFMAXALPHA.GE.6)THEN
               EVF = EVMASS(5,ASI,ASTHM(6)) * DECOUP("UP",6,LN)
     1             * EVMASS(6,ASTH(6),ASQ)
            ELSEIF(MU2.GE.M2TH(5).AND.NFMAXALPHA.GE.5)THEN
               EVF = EVMASS(5,ASI,ASQ) 
            ELSEIF(MU2.GE.M2TH(4).AND.NFMAXALPHA.GE.4)THEN
               EVF = DECOUP("DW",5,LN) / EVMASS(4,ASQ,ASIM) 
            ELSE
               EVF = DECOUP("DW",4,LN) / EVMASS(3,ASQ,ASTHM(4)) 
     1             * DECOUP("DW",5,LN) / EVMASS(4,ASTH(4),ASIM) 
            ENDIF
*        Top
         ELSEIF(IM.EQ.6)THEN
            IF(MU2.GE.M2TH(6).AND.NFMAXALPHA.GE.6)THEN
               EVF = EVMASS(6,ASI,ASQ) 
            ELSEIF(MU2.GE.M2TH(5).AND.NFMAXALPHA.GE.5)THEN
               EVF = DECOUP("DW",6,LN) / EVMASS(5,ASQ,ASIM) 
            ELSEIF(MU2.GE.M2TH(4).AND.NFMAXALPHA.GE.4)THEN
               EVF = DECOUP("DW",5,LN) / EVMASS(4,ASQ,ASIM) 
     1             * DECOUP("DW",6,LN) / EVMASS(5,ASTH(5),ASTHM(6)) 
            ELSE
               EVF = DECOUP("DW",4,LN) / EVMASS(3,ASQ,ASIM)
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
         ASTH2(4)  = A_QCD(M2TH(4))**2
         ASTH2(5)  = A_QCD(M2TH(5))**2
         ASTH2(6)  = A_QCD(M2TH(6))**2
*
         ASTHM2(4) = A_QCD(M2TH(4)-EPS)**2
         ASTHM2(5) = A_QCD(M2TH(5)-EPS)**2
         ASTHM2(6) = A_QCD(M2TH(6)-EPS)**2
*
         IF(DIR.EQ."DW")THEN
            DECOUP = 1D0 + ASTH2(IM) * ( 89D0 / 27D0 - 20D0 / 9D0 * LN 
     1                                 + 4D0 / 3D0 * LN**2 )
         ELSEIF(DIR.EQ."UP")THEN
            DECOUP = 1D0 - ASTHM2(IM) * ( 89D0 / 27D0 - 20D0 / 9D0 * LN 
     1                                  + 4D0 / 3D0 * LN**2 )
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
     1          - 140d0 / 81d0 * nf**2
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
**
*     Input Variables
*
      integer i
      double precision Q
**
*     Internal Variables
*
      double precision evmass
      double precision a_QCD,as0,as
**
*     Output Variables
*
      double precision MassQSplit
*
      as0 = a_QCD(q2th(i))
      as  = a_QCD(Q**2)
      MassQSplit = dsqrt(m2q(i)) * evmass(i,as0,as) - Q
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
      integer i,j
      double precision MassQSplit,zriddr
      double precision x1,x2
      double precision acc,window
      parameter(acc=1d-10)
      parameter(window=2d0)
      external MassQSplit
*
*     Iterate more times to make sure that the running of alphas
*     is computed using the correct thresholds.
*
      do j=1,4
         do i=4,6
            if(q2th(i).ne.m2q(i))then
               x1 = dsqrt(q2th(i)) - window
               x2 = dsqrt(q2th(i)) + window
               m2ph(i) = zriddr(MassQSplit,i,x1,x2,acc)**2
            endif
         enddo
         call ComputeHeavyQuarkThresholds
      enddo
*
      return
      end
