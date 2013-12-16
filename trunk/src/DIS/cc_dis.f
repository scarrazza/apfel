************************************************************************
*
*     cc_dis.f:
*
*     x-space computation of CC observables.
*
************************************************************************
      SUBROUTINE CC_DIS(X,QI,QF,RAP,POL,SCHEME,PTO,PDFSET,IREP,TARGET,
     1                  PROJ,F2,F3,FL,SIGMA)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
      include "../commons/scale.h"
      include "../commons/hmass.h"
      include "../commons/neut.h"
      include "../commons/vfns.h"
      include "../commons/heavythr.h"
      include "../commons/beta.h"
      include "../commons/jpt.h"
      include "../commons/asq.h"
      include "../commons/nf.h"
**
*     Input Varibles
*
      INTEGER PTO,IREP
      DOUBLE PRECISION X,QI,QF,RAP,POL
      CHARACTER*5  SCHEME
      CHARACTER*53 PDFSET
      CHARACTER*9  TARGET
      CHARACTER*12 PROJ
**
*     Internal Variables
*
      INTEGER IPDF,IBENCH,IGM,I
      INTEGER HQMASS
      INTEGER F3SGN
      DOUBLE PRECISION XI,Y,YPL,YMN,POLFACT
      DOUBLE PRECISION Q2,LAMBDA,KA,DAMP(4:6)
      DOUBLE PRECISION YP,YM
      DOUBLE PRECISION ASQ2
      DOUBLE PRECISION YPCC,FACT_CC,CONV
      DOUBLE PRECISION MW,MW2,MN,GF,GF2,PI
      DOUBLE PRECISION ZNUCL,ANUCL,FRAT
      DOUBLE PRECISION V_UD,V_US,V_UD2,V_US2
      DOUBLE PRECISION V_CD,V_CS,V_CD2,V_CS2
      DOUBLE PRECISION V_UB,V_CB,V_UB2,V_CB2
      DOUBLE PRECISION V_TD,V_TS,V_TB,V_TD2,V_TS2,V_TB2
      DOUBLE PRECISION DGAUSS,A,B,C,D,EPS,EPS2
      DOUBLE PRECISION XD,XU,XS,XC,XB,XT,XGL
      DOUBLE PRECISION LO,LO1,LO2,LO3
      DOUBLE PRECISION LOP,LOM,DLO,DY
      DOUBLE PRECISION INTEG1(0:6),INTEG2(0:6),INTEG3(0:6)
      DOUBLE PRECISION CONST1,CONST2,CONST3,H1
      DOUBLE PRECISION NLOQ1,NLOG1,NLO1
      DOUBLE PRECISION NLOQ2,NLOG2,NLO2
      DOUBLE PRECISION NLOQ3,NLOG3,NLO3
      DOUBLE PRECISION NNLOQ1,NNLOG1,NNLO1
      DOUBLE PRECISION NNLOQ2,NNLOG2,NNLO2
      DOUBLE PRECISION NNLOQ3,NNLOG3,NNLO3
      DOUBLE PRECISION PLUS_TERM1,PLUS_TERM2,PLUS_TERM3
      DOUBLE PRECISION C2G2C,C2NN2C,CLNN2C,C3NM2C
      DOUBLE PRECISION F1CG(3),F2CG(3),F3CG(3),FLCG(3),SIGMACG(3)
      DOUBLE PRECISION F1BG(3),F2BG(3),F3BG(3),FLBG(3),SIGMABG(3)
      DOUBLE PRECISION F1TG(3),F2TG(3),F3TG(3),FLTG(3),SIGMATG(3)
      DOUBLE PRECISION F2P,F3P,FLP,SIGMAP
      DOUBLE PRECISION F1L,F2L,F3L,FLL,SIGMAL
      DOUBLE PRECISION F1C,F2C,F3C,FLC,SIGMAC
      DOUBLE PRECISION F1B,F2B,F3B,FLB,SIGMAB
      DOUBLE PRECISION F1T,F2T,F3T,FLT,SIGMAT

      DOUBLE PRECISION XGL_XC1G
      EXTERNAL XGL_XC1G
      DOUBLE PRECISION XPS_XC1Q_UNPLUS
      EXTERNAL XPS_XC1Q_UNPLUS
      DOUBLE PRECISION XD_XC1Q_UNPLUS
      EXTERNAL XD_XC1Q_UNPLUS
      DOUBLE PRECISION XD_XC1Q_PLUS
      EXTERNAL XD_XC1Q_PLUS
      DOUBLE PRECISION XU_XC1Q_UNPLUS
      EXTERNAL XU_XC1Q_UNPLUS
      DOUBLE PRECISION XU_XC1Q_PLUS
      EXTERNAL XU_XC1Q_PLUS
      DOUBLE PRECISION XS_XC1Q_UNPLUS
      EXTERNAL XS_XC1Q_UNPLUS
      DOUBLE PRECISION XS_XC1Q_PLUS
      EXTERNAL XS_XC1Q_PLUS
      DOUBLE PRECISION XC_XC1Q_UNPLUS
      EXTERNAL XC_XC1Q_UNPLUS
      DOUBLE PRECISION XC_XC1Q_PLUS
      EXTERNAL XC_XC1Q_PLUS
      DOUBLE PRECISION XB_XC1Q_UNPLUS
      EXTERNAL XB_XC1Q_UNPLUS
      DOUBLE PRECISION XB_XC1Q_PLUS
      EXTERNAL XB_XC1Q_PLUS
      DOUBLE PRECISION XT_XC1Q_UNPLUS
      EXTERNAL XT_XC1Q_UNPLUS
      DOUBLE PRECISION XT_XC1Q_PLUS
      EXTERNAL XT_XC1Q_PLUS

      DOUBLE PRECISION C1Q_PLUS_WRAP
      EXTERNAL C1Q_PLUS_WRAP

      DOUBLE PRECISION XGL_XC2G
      EXTERNAL XGL_XC2G
      DOUBLE PRECISION XPS_XC2Q_UNPLUS
      EXTERNAL XPS_XC2Q_UNPLUS
      DOUBLE PRECISION XD_XC2Q_UNPLUS
      EXTERNAL XD_XC2Q_UNPLUS
      DOUBLE PRECISION XD_XC2Q_PLUS
      EXTERNAL XD_XC2Q_PLUS
      DOUBLE PRECISION XU_XC2Q_UNPLUS
      EXTERNAL XU_XC2Q_UNPLUS
      DOUBLE PRECISION XU_XC2Q_PLUS
      EXTERNAL XU_XC2Q_PLUS
      DOUBLE PRECISION XS_XC2Q_UNPLUS
      EXTERNAL XS_XC2Q_UNPLUS
      DOUBLE PRECISION XS_XC2Q_PLUS
      EXTERNAL XS_XC2Q_PLUS
      DOUBLE PRECISION XC_XC2Q_UNPLUS
      EXTERNAL XC_XC2Q_UNPLUS
      DOUBLE PRECISION XC_XC2Q_PLUS
      EXTERNAL XC_XC2Q_PLUS
      DOUBLE PRECISION XB_XC2Q_UNPLUS
      EXTERNAL XB_XC2Q_UNPLUS
      DOUBLE PRECISION XB_XC2Q_PLUS
      EXTERNAL XB_XC2Q_PLUS
      DOUBLE PRECISION XT_XC2Q_UNPLUS
      EXTERNAL XT_XC2Q_UNPLUS
      DOUBLE PRECISION XT_XC2Q_PLUS
      EXTERNAL XT_XC2Q_PLUS

      DOUBLE PRECISION C2Q_PLUS_WRAP
      EXTERNAL C2Q_PLUS_WRAP

      DOUBLE PRECISION XGL_XC3G
      EXTERNAL XGL_XC3G
      DOUBLE PRECISION XPS_XC3Q_UNPLUS
      EXTERNAL XPS_XC3Q_UNPLUS
      DOUBLE PRECISION XD_XC3Q_UNPLUS
      EXTERNAL XD_XC3Q_UNPLUS
      DOUBLE PRECISION XD_XC3Q_PLUS
      EXTERNAL XD_XC3Q_PLUS
      DOUBLE PRECISION XU_XC3Q_UNPLUS
      EXTERNAL XU_XC3Q_UNPLUS
      DOUBLE PRECISION XU_XC3Q_PLUS
      EXTERNAL XU_XC3Q_PLUS
      DOUBLE PRECISION XS_XC3Q_UNPLUS
      EXTERNAL XS_XC3Q_UNPLUS
      DOUBLE PRECISION XS_XC3Q_PLUS
      EXTERNAL XS_XC3Q_PLUS
      DOUBLE PRECISION XC_XC3Q_UNPLUS
      EXTERNAL XC_XC3Q_UNPLUS
      DOUBLE PRECISION XC_XC3Q_PLUS
      EXTERNAL XC_XC3Q_PLUS
      DOUBLE PRECISION XB_XC3Q_UNPLUS
      EXTERNAL XB_XC3Q_UNPLUS
      DOUBLE PRECISION XB_XC3Q_PLUS
      EXTERNAL XB_XC3Q_PLUS
      DOUBLE PRECISION XT_XC3Q_UNPLUS
      EXTERNAL XT_XC3Q_UNPLUS
      DOUBLE PRECISION XT_XC3Q_PLUS
      EXTERNAL XT_XC3Q_PLUS

      DOUBLE PRECISION C3Q_PLUS_WRAP
      EXTERNAL C3Q_PLUS_WRAP

      CHARACTER*5 WVFNS,IVFNS(3)

      PARAMETER(EPS  = 1D-5)
      PARAMETER(EPS2 = 1D-10)
      PARAMETER(DY   = 1D-5)
      PARAMETER(PI   = 3.14159265358979D0)
      PARAMETER(MN   = 0.9389D0)
      PARAMETER(MW   = 80.399D0)
      PARAMETER(GF   = 1.16637D-5)
      PARAMETER(CONV = 3.893793D10)
      PARAMETER(V_UD = 0.97428D0,  V_US = 0.22530D0)
      PARAMETER(V_CD = 0.22520D0,  V_CS = 0.97345D0)
      PARAMETER(V_UB = 0.003470D0, V_CB = 0.041000D0)
      PARAMETER(V_TD = 0.00862D0,  V_TS = 0.04030D0, V_TB = 0.999152D0)
**
*     Output Variables
*
      DOUBLE PRECISION F2(3:7),F3(3:7),FL(3:7),SIGMA(3:7)
*
*     Beta function coefficients
*
      DO I=3,6
         BETA0(I) = ( 33D0 - 2D0 * I ) / 3D0
         BETA1(I) = 102D0 - 38D0 / 3D0 * I
         BETA2(I) = 2857D0 / 2D0 - 5033D0 / 18D0 * I 
     1            + 325D0 / 54D0 * I**2D0
         B1(I)    = BETA1(I) / BETA0(I)
         B2(I)    = BETA2(I) / BETA0(I)
      ENDDO
*
*     CKM matrix elements
*
      V_UD2 = V_UD**2D0
      V_US2 = V_US**2D0
      V_CD2 = V_CD**2D0
      V_CS2 = V_CS**2D0
      V_UB2 = V_UB**2D0
      V_CB2 = V_CB**2D0
      V_TD2 = V_TD**2D0
      V_TS2 = V_TS**2D0
      V_TB2 = V_TB**2D0
*
      IVFNS(1) = "FFNS"
      IVFNS(2) = "ZMVN"
      IVFNS(3) = "FFN0"
*
*     Set parameters
*
      Q  = QF
      Q2 = QF * QF
      Q0 = QI - EPS2
*
      YM = 1D0 - ( 1D0 - RAP )**2D0
      YP = 1D0 + ( 1D0 - RAP )**2D0
*
*     Define Observable
*
      IF(TARGET(1:6).EQ."PROTON")THEN
         IBENCH = 2
         ZNUCL  = 1D0
         ANUCL  = 1D0
      ELSEIF(TARGET(1:7).EQ."NEUTRON")THEN
         IBENCH = 2
         ZNUCL  = 0D0
         ANUCL  = 1D0
      ELSEIF(TARGET(1:9).EQ."ISOSCALAR")THEN
         IBENCH = 1
         ZNUCL  = 1D0
         ANUCL  = 2D0
      ELSEIF(TARGET(1:4).EQ."IRON")THEN
         IBENCH = 0
         ZNUCL  = 23.403D0
         ANUCL  = 49.618D0
      ELSE
         WRITE(6,*) "In cc_dis.f:"
         WRITE(6,*) "Unknown value target = ",TARGET
         CALL EXIT(-10)
      ENDIF
      FRAT = ZNUCL / ANUCL
*
      IF(PROJ(1:8).EQ."ELECTRON")THEN
         NEUT  = - 1
         F3SGN = 1
      ELSEIF(PROJ(1:8).EQ."POSITRON")THEN
         NEUT  = 1
         F3SGN = 1
      ELSEIF(PROJ(1:8).EQ."NEUTRINO")THEN
         NEUT  = 1
         F3SGN = - 1
      ELSEIF(PROJ(1:12).EQ."ANTINEUTRINO")THEN
         NEUT  = - 1
         F3SGN = - 1
      ELSE
         WRITE(6,*) "In cc_dis.f:"
         WRITE(6,*) "Unknown value projectile = ",PROJ
         CALL EXIT(-10)
      ENDIF
*
      IF(SCHEME(1:4).EQ."ZMVN") THEN
         VFNS = "ZMVN"
         WVFNS = VFNS
      ELSEIF(SCHEME(1:4).EQ."FFNS") THEN
         VFNS = "FFNS"
         WVFNS = VFNS         
      ELSEIF(SCHEME(1:4).EQ."FFN0") THEN
         VFNS = "FFN0"
         WVFNS = VFNS        
      ELSEIF(SCHEME(1:5).EQ."FONLL") THEN
         VFNS = "FONLL"
         WVFNS = VFNS          
      ELSE
         WRITE(6,*) "In nc_dis.f:"
         WRITE(6,*) "Invalid evolution scheme: VFNS = ",VFNS
         CALL EXIT(-10)
      ENDIF
*
      IPT = PTO
      IF(IPT.NE.0.AND.IPT.NE.1.AND.PTO.NE.2)THEN
         WRITE(6,*) "In cc_dis.f:"
         WRITE(6,*) "Invalid perturbative order: IPT =",IPT
         CALL EXIT(-10)
      ENDIF
*
      IF(IBENCH.EQ.0)THEN
         MW2     = MW**2D0
         FACT_CC = 100D0 / ( 2D0 * ( 1D0 +  Q2 / MW2 )**2D0 )
         YPCC    = YP - 2D0 * ( MN * X * RAP )**2D0 / Q2
      ELSEIF(IBENCH.EQ.1)THEN
         MW2     = MW**2D0
         GF2     = GF**2D0
         FACT_CC = CONV * GF2 * MN 
     1           / ( 2D0 * PI * ( 1D0 + ( Q2 / MW2 ) )**2D0)
         YPCC    = YP - 2D0 * ( MN * X * RAP )**2D0 / Q2
      ELSEIF(IBENCH.EQ.2)THEN
         FACT_CC = 1D0 / 4D0
         YPCC    = YP
      ENDIF
*
      HQMASS = 0 ! Only pole mass for the time being
*
*     Initialize PDFs and parameters
*
      CALL PDFPARAM(PDFSET,IREP)
*
      ASQ2 = ASQ * ASQ
*
      DO I=4,6
         IF(Q2.GT.Q2TH(I))THEN
            DAMP(I) = ( 1D0 - Q2TH(I) / Q2 )**2D0
         ELSE
            DAMP(I) = 0D0
         ENDIF
      ENDDO
*
************************************************************************
*     LIGHT
************************************************************************
*
*     Always ZM-VFNS
*
      VFNS = "ZMVN"
*
      Y      = X
      CONST1 = - 4D0 / 3D0 * ( PI**2D0 / 3D0 + 9D0 / 2D0 )
      CONST2 = CONST1
      CONST3 = CONST1
      XBJ    = X
*
*     LO
*
      IH = 4
      LO = 2D0 * V_UD2 * ( FRAT * XD(Y,Q) + (1D0-FRAT) * XU(Y,Q) )
     1   + 2D0 * V_US2 * XS(Y,Q)
*
      IH = 5
      LO1 = LO + 2D0 * ( V_UD2 + V_US2 ) * ( FRAT * XU(Y,Q) 
     1    + ( 1D0 - FRAT ) * XD(Y,Q) )
      LO2 = LO + 2D0 * ( V_UD2 + V_US2 ) * ( FRAT * XU(Y,Q) 
     1    + ( 1D0 - FRAT ) * XD(Y,Q) )
      LO3 = LO - 2D0 * ( V_UD2 + V_US2 ) * ( FRAT * XU(Y,Q) 
     1    + ( 1D0 - FRAT ) * XD(Y,Q) )
*
*     Bounds of the integrals
*
      A = Y
      B = 1D0
      C = 0D0
      D = Y
*
*     NLO
*
      NLO1 = 0D0
      NLO2 = 0D0
      NLO3 = 0D0
      IF(PTO.GE.1)THEN
*
         IPT = 1
*
*        Gluon
*
         INTEG1(0) = 2D0 * ( V_UD2 + V_US2 ) * DGAUSS(XGL_XC1G,A,B,EPS)
         INTEG2(0) = 2D0 * ( V_UD2 + V_US2 ) * DGAUSS(XGL_XC2G,A,B,EPS)
         INTEG3(0) = 2D0 * ( V_UD2 + V_US2 ) * DGAUSS(XGL_XC3G,A,B,EPS)
*
*        Quark
*
*        Contribution coming from the plus prescription in an integral between x and 1
*
         PLUS_TERM1 = DGAUSS(C1Q_PLUS_WRAP,C,D,EPS)
         PLUS_TERM2 = DGAUSS(C2Q_PLUS_WRAP,C,D,EPS)
         PLUS_TERM3 = DGAUSS(C3Q_PLUS_WRAP,C,D,EPS)
*
*        F1
*
         IH = 4
         INTEG1(1) = 2D0 * V_UD2 * FRAT 
     1             * ( DGAUSS(XD_XC1Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XD_XC1Q_PLUS,A,B,EPS) 
     3             + ( CONST1 - PLUS_TERM1 ) * XD(Y,Q) )
         INTEG1(2) = 2D0 * V_UD2 * ( 1D0 - FRAT )
     1             * ( DGAUSS(XU_XC1Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XU_XC1Q_PLUS,A,B,EPS) 
     3             + ( CONST1 - PLUS_TERM1 ) * XU(Y,Q) )
         INTEG1(3) = 2D0 * V_US2
     1             * ( DGAUSS(XS_XC1Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XS_XC1Q_PLUS,A,B,EPS) 
     3             + ( CONST1 - PLUS_TERM1 ) * XS(Y,Q) )
         IH = 5
         INTEG1(4) = 2D0 * ( V_UD2 + V_US2 ) * FRAT 
     1             * ( DGAUSS(XU_XC1Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XU_XC1Q_PLUS,A,B,EPS) 
     3             + ( CONST1 - PLUS_TERM1 ) * XU(Y,Q) )
         INTEG1(5) = 2D0 * ( V_UD2 + V_US2 ) * ( 1D0 - FRAT )
     1             * ( DGAUSS(XD_XC1Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XD_XC1Q_PLUS,A,B,EPS) 
     3             + ( CONST1 - PLUS_TERM1 ) * XD(Y,Q) )
         INTEG1(6) = 0D0
*
*        F2
*
         IH = 4
         INTEG2(1) = 2D0 * V_UD2 * FRAT 
     1             * ( DGAUSS(XD_XC2Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XD_XC2Q_PLUS,A,B,EPS) 
     3             + ( CONST2 - PLUS_TERM2 ) * XD(Y,Q) )
         INTEG2(2) = 2D0 * V_UD2 * ( 1D0 - FRAT )
     1             * ( DGAUSS(XU_XC2Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XU_XC2Q_PLUS,A,B,EPS) 
     3             + ( CONST2 - PLUS_TERM2 ) * XU(Y,Q) )
         INTEG2(3) = 2D0 * V_US2
     1             * ( DGAUSS(XS_XC2Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XS_XC2Q_PLUS,A,B,EPS) 
     3             + ( CONST2 - PLUS_TERM2 ) * XS(Y,Q) )
         IH = 5
         INTEG2(4) = 2D0 * ( V_UD2 + V_US2 ) * FRAT 
     1             * ( DGAUSS(XU_XC2Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XU_XC2Q_PLUS,A,B,EPS) 
     3             + ( CONST2 - PLUS_TERM2 ) * XU(Y,Q) )
         INTEG2(5) = 2D0 * ( V_UD2 + V_US2 ) * ( 1D0 - FRAT )
     1             * ( DGAUSS(XD_XC2Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XD_XC2Q_PLUS,A,B,EPS) 
     3             + ( CONST2 - PLUS_TERM2 ) * XD(Y,Q) )
         INTEG2(6) = 0D0
*
*        F3
*
         IH = 4
         INTEG3(1) = 2D0 * V_UD2 * FRAT 
     1             * ( DGAUSS(XD_XC3Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XD_XC3Q_PLUS,A,B,EPS) 
     3             + ( CONST3 - PLUS_TERM3 ) * XD(Y,Q) )
         INTEG3(2) = 2D0 * V_UD2 * ( 1D0 - FRAT )
     1             * ( DGAUSS(XU_XC3Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XU_XC3Q_PLUS,A,B,EPS) 
     3             + ( CONST3 - PLUS_TERM3 ) * XU(Y,Q) )
         INTEG3(3) = 2D0 * V_US2
     1             * ( DGAUSS(XS_XC3Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XS_XC3Q_PLUS,A,B,EPS) 
     3             + ( CONST3 - PLUS_TERM3 ) * XS(Y,Q) )
         IH = 5
         INTEG3(4) = - 2D0 * ( V_UD2 + V_US2 ) * FRAT 
     1             * ( DGAUSS(XU_XC3Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XU_XC3Q_PLUS,A,B,EPS) 
     3             + ( CONST3 - PLUS_TERM3 ) * XU(Y,Q) )
         INTEG3(5) = - 2D0 * ( V_UD2 + V_US2 ) * ( 1D0 - FRAT )
     1             * ( DGAUSS(XD_XC3Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XD_XC3Q_PLUS,A,B,EPS) 
     3             + ( CONST3 - PLUS_TERM3 ) * XD(Y,Q) )
         INTEG3(6) = 0D0
*
         NLOG1 = 0D0
         NLOG2 = 0D0
         NLOG3 = 0D0
*
         NLOQ1 = 0D0
         NLOQ2 = 0D0
         NLOQ3 = 0D0
*
         IF(IPT.GE.1)THEN
            NLOG1 = 2D0 * INTEG1(0)
            NLOG2 = 2D0 * INTEG2(0)
            NLOG3 = 2D0 * INTEG3(0)
            DO IPDF=1,6
               NLOQ1 = NLOQ1 + INTEG1(IPDF)
               NLOQ2 = NLOQ2 + INTEG2(IPDF)
               NLOQ3 = NLOQ3 + INTEG3(IPDF)
            ENDDO
         ENDIF
*
         NLO1 = 2D0 * ( NLOG1 + NLOQ1 )
         NLO2 = 2D0 * ( NLOG2 + NLOQ2 )
         NLO3 = 2D0 * ( NLOG3 + NLOQ3 )
      ENDIF
*
*     NNLO
*
*     Initialize
*
      NNLO1 = 0D0
      NNLO2 = 0D0
      NNLO3 = 0D0
      IF(PTO.GE.2)THEN
*
         IPT = 2
*
*        Gluon + Pure Singlet
*
         INTEG1(0) = ( V_UD2 + V_US2 ) 
     1             * ( DGAUSS(XGL_XC1G,A,B,EPS) + C2G2C(X,1) * XGL(X,Q)
     2             +   DGAUSS(XPS_XC1Q_UNPLUS,A,B,EPS) )
         INTEG2(0) = ( V_UD2 + V_US2 ) 
     1             * ( DGAUSS(XGL_XC2G,A,B,EPS) + C2G2C(X,1) * XGL(X,Q)
     2             +   DGAUSS(XPS_XC2Q_UNPLUS,A,B,EPS) )
         INTEG3(0) = 0D0
*
*        Quark
*
*        Contribution coming from the plus prescription in an integral between x and 1
*
         PLUS_TERM1 = C2NN2C(X,NF) - CLNN2C(X,NF)
         PLUS_TERM2 = C2NN2C(X,NF)
         PLUS_TERM3 = C3NM2C(X,NF)
*
*        F1
*
         IH = 4
         INTEG1(1) = V_UD2 * FRAT 
     1             * ( DGAUSS(XD_XC1Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XD_XC1Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM1 * XD(Y,Q) )
         INTEG1(2) = V_UD2 * ( 1D0 - FRAT )
     1             * ( DGAUSS(XU_XC1Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XU_XC1Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM1 * XU(Y,Q) )
         INTEG1(3) = V_US2
     1             * ( DGAUSS(XS_XC1Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XS_XC1Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM1 * XS(Y,Q) )
         IH = 5
         INTEG1(4) = ( V_UD2 + V_US2 ) * FRAT 
     1             * ( DGAUSS(XU_XC1Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XU_XC1Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM1 * XU(Y,Q) )
         INTEG1(5) = 2D0 * ( V_UD2 + V_US2 ) * ( 1D0 - FRAT )
     1             * ( DGAUSS(XD_XC1Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XD_XC1Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM1 * XD(Y,Q) )
         INTEG1(6) = 0D0
*
*        F2
*
         IH = 4
         INTEG2(1) = V_UD2 * FRAT 
     1             * ( DGAUSS(XD_XC2Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XD_XC2Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM2 * XD(Y,Q) )
         INTEG2(2) = V_UD2 * ( 1D0 - FRAT )
     1             * ( DGAUSS(XU_XC2Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XU_XC2Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM2 * XU(Y,Q) )
         INTEG2(3) = V_US2
     1             * ( DGAUSS(XS_XC2Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XS_XC2Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM2 * XS(Y,Q) )
         IH = 5
         INTEG2(4) = ( V_UD2 + V_US2 ) * FRAT 
     1             * ( DGAUSS(XU_XC2Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XU_XC2Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM2 * XU(Y,Q) )
         INTEG2(5) = ( V_UD2 + V_US2 ) * ( 1D0 - FRAT )
     1             * ( DGAUSS(XD_XC2Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XD_XC2Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM2 * XD(Y,Q) )
         INTEG2(6) = 0D0
*
*        F3
*
         IH = 4
         INTEG3(1) = V_UD2 * FRAT 
     1             * ( DGAUSS(XD_XC3Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XD_XC3Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM3 * XD(Y,Q) )
         INTEG3(2) = V_UD2 * ( 1D0 - FRAT )
     1             * ( DGAUSS(XU_XC3Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XU_XC3Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM3 * XU(Y,Q) )
         INTEG3(3) = V_US2
     1             * ( DGAUSS(XS_XC3Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XS_XC3Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM3 * XS(Y,Q) )
         IH = 5
         INTEG3(4) = - ( V_UD2 + V_US2 ) * FRAT 
     1             * ( DGAUSS(XU_XC3Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XU_XC3Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM3 * XU(Y,Q) )
         INTEG3(5) = - ( V_UD2 + V_US2 ) * ( 1D0 - FRAT )
     1             * ( DGAUSS(XD_XC3Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XD_XC3Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM3 * XD(Y,Q) )
         INTEG3(6) = 0D0
*
         NNLOG1 = 0D0
         NNLOG2 = 0D0
         NNLOG3 = 0D0
*
         NNLOQ1 = 0D0
         NNLOQ2 = 0D0
         NNLOQ3 = 0D0
*
         IF(IPT.GE.1)THEN
            NNLOG1 = INTEG1(0)
            NNLOG2 = INTEG2(0)
            NNLOG3 = INTEG3(0)
            DO IPDF=1,6
               NNLOQ1 = NNLOQ1 + INTEG1(IPDF)
               NNLOQ2 = NNLOQ2 + INTEG2(IPDF)
               NNLOQ3 = NNLOQ3 + INTEG3(IPDF)
            ENDDO
         ENDIF
*
         NNLO1 = 2D0 * ( NNLOG1 + NNLOQ1 )
         NNLO2 = 2D0 * ( NNLOG2 + NNLOQ2 )
         NNLO3 = 2D0 * ( NNLOG3 + NNLOQ3 )
      ENDIF
*
      F1L = LO1 + ASQ * NLO1 + ASQ2 * NNLO1
      F2L = LO2 + ASQ * NLO2 + ASQ2 * NNLO2
      F3L = F3SGN * ( LO3 + ASQ * NLO3  + ASQ2 * NNLO3 )
*
*     FL
*
      FLL = F2L - F1L
*
*     Evaluation of the cross section
*
      SIGMAL = FACT_CC * ( YPCC * F2L - RAP**2D0 * FLL - YM * F3L )
*
      VFNS = WVFNS
*
************************************************************************
*     CHARM
************************************************************************
*
      IH  = 4
      Q2H = Q2TH(IH)
*
*     Charm contribution evaluated only for energies over the threshold
*
      IF(Q2.LT.Q2H) GOTO 200
*
*     In case GM scheme is chosen
*
      IGM = 1
 405  IF(WVFNS.EQ."FONLL")THEN
         VFNS = IVFNS(IGM)
      ENDIF
*
      LAMBDA = 1D0 / ( 1D0 + Q2H / Q2 )
      XI     = X / LAMBDA
      IF(XI.GE.1D0) XI = 0.99999999D0
      IF(VFNS.EQ."ZMVN".OR.VFNS.EQ."FFN0")THEN
         Y      = X
         CONST1 = - 4D0 / 3D0 * ( PI**2D0 / 3D0 + 9D0 / 2D0 )
         CONST2 = CONST1
         CONST3 = CONST1
         XBJ    = X
      ELSEIF(VFNS.EQ."FFNS")THEN
         Y      = XI
         KA     = ( 1D0 - LAMBDA ) * DLOG( 1D0 - LAMBDA ) / LAMBDA
         CONST1 = - 4D0 / 3D0 * ( 4D0 + 1D0 / ( 2D0 * LAMBDA )
     1          + PI**2D0 / 3D0 + ( 1D0 + 3D0 * LAMBDA ) * KA 
     2          / (2D0 * LAMBDA ) )
         CONST2 = - 4D0 / 3D0 * ( 4D0 + 1D0 / ( 2D0 * LAMBDA )
     1          + PI**2D0 / 3D0 + ( 1D0 + LAMBDA ) * KA 
     2          / (2D0 * LAMBDA ) )
         CONST3 = CONST1
         XBJ    = XI
      ENDIF
*
*     LO
*
      LO = 2D0 * V_CD2 * ( FRAT * XD(Y,Q) + (1D0-FRAT) * XU(Y,Q) )
     1   + 2D0 * V_CS2 * XS(Y,Q)
*
      LO1 = LO
      LO2 = LO
      LO3 = LO
      IF(VFNS.EQ."ZMVN")THEN
         LO1 = LO1 + 2D0 * ( V_CD2 + V_CS2 ) * XC(Y,Q)
         LO2 = LO2 + 2D0 * ( V_CD2 + V_CS2 ) * XC(Y,Q)
         LO3 = LO3 - 2D0 * ( V_CD2 + V_CS2 ) * XC(Y,Q)
      ENDIF
*
*     Bounds of the integrals
*
      A = Y
      B = 1D0
      C = 0D0
      D = Y
*
*     NLO
*
      NLO1 = 0D0
      NLO2 = 0D0
      NLO3 = 0D0
      IF(PTO.GE.1)THEN
*
         IPT = 1
*
*        Gluon
*
         INTEG1(0) = 2D0 * ( V_CD2 + V_CS2 ) * DGAUSS(XGL_XC1G,A,B,EPS)
         INTEG2(0) = 2D0 * ( V_CD2 + V_CS2 ) * DGAUSS(XGL_XC2G,A,B,EPS)
         INTEG3(0) = 2D0 * ( V_CD2 + V_CS2 ) * DGAUSS(XGL_XC3G,A,B,EPS)
*
*        Quark
*
*        Contribution coming from the plus prescription in an integral between x and 1
*
         PLUS_TERM1 = DGAUSS(C1Q_PLUS_WRAP,C,D,EPS)
         PLUS_TERM2 = DGAUSS(C2Q_PLUS_WRAP,C,D,EPS)
         PLUS_TERM3 = DGAUSS(C3Q_PLUS_WRAP,C,D,EPS)
*
*        F1
*
         INTEG1(1) = 2D0 * V_CD2 * FRAT 
     1             * ( DGAUSS(XD_XC1Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XD_XC1Q_PLUS,A,B,EPS) 
     3             + ( CONST1 - PLUS_TERM1 ) * XD(Y,Q) )
         INTEG1(2) = 2D0 * V_CD2 * ( 1D0 - FRAT )
     1             * ( DGAUSS(XU_XC1Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XU_XC1Q_PLUS,A,B,EPS) 
     3             + ( CONST1 - PLUS_TERM1 ) * XU(Y,Q) )
         INTEG1(3) = 2D0 * V_CS2
     1             * ( DGAUSS(XS_XC1Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XS_XC1Q_PLUS,A,B,EPS) 
     3             + ( CONST1 - PLUS_TERM1 ) * XS(Y,Q) )
         IF(VFNS.EQ."ZMVN")THEN
            INTEG1(4) = 2D0 * ( V_CD2 + V_CS2 )
     1                * ( DGAUSS(XC_XC1Q_UNPLUS,A,B,EPS)
     2                +   DGAUSS(XC_XC1Q_PLUS,A,B,EPS) 
     3                + ( CONST1 - PLUS_TERM1 ) * XC(Y,Q) )
         ELSE
            INTEG1(4) = 0D0
         ENDIF
         INTEG1(5) = 0D0
         INTEG1(6) = 0D0
*
*        F2
*
         INTEG2(1) = 2D0 * V_CD2 * FRAT 
     1             * ( DGAUSS(XD_XC2Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XD_XC2Q_PLUS,A,B,EPS) 
     3             + ( CONST2 - PLUS_TERM2 ) * XD(Y,Q) )
         INTEG2(2) = 2D0 * V_CD2 * ( 1D0 - FRAT )
     1             * ( DGAUSS(XU_XC2Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XU_XC2Q_PLUS,A,B,EPS) 
     3             + ( CONST2 - PLUS_TERM2 ) * XU(Y,Q) )
         INTEG2(3) = 2D0 * V_CS2
     1             * ( DGAUSS(XS_XC2Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XS_XC2Q_PLUS,A,B,EPS) 
     3             + ( CONST2 - PLUS_TERM2 ) * XS(Y,Q) )
         IF(VFNS.EQ."ZMVN")THEN
            INTEG2(4) = 2D0 * ( V_CD2 + V_CS2 )
     1                * ( DGAUSS(XC_XC2Q_UNPLUS,A,B,EPS)
     2                +   DGAUSS(XC_XC2Q_PLUS,A,B,EPS) 
     3                + ( CONST2 - PLUS_TERM2 ) * XC(Y,Q) )
         ELSE
            INTEG2(4) = 0D0
         ENDIF
         INTEG2(5) = 0D0
         INTEG2(6) = 0D0
*
*        F3
*
         INTEG3(1) = 2D0 * V_CD2 * FRAT 
     1             * ( DGAUSS(XD_XC3Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XD_XC3Q_PLUS,A,B,EPS) 
     3             + ( CONST3 - PLUS_TERM3 ) * XD(Y,Q) )
         INTEG3(2) = 2D0 * V_CD2 * ( 1D0 - FRAT )
     1             * ( DGAUSS(XU_XC3Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XU_XC3Q_PLUS,A,B,EPS) 
     3             + ( CONST3 - PLUS_TERM3 ) * XU(Y,Q) )
         INTEG3(3) = 2D0 * V_CS2
     1             * ( DGAUSS(XS_XC3Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XS_XC3Q_PLUS,A,B,EPS) 
     3             + ( CONST3 - PLUS_TERM3 ) * XS(Y,Q) )
         IF(VFNS.EQ."ZMVN")THEN
            INTEG3(4) = - 2D0 * ( V_CD2 + V_CS2 )
     1                * ( DGAUSS(XC_XC3Q_UNPLUS,A,B,EPS)
     2                +   DGAUSS(XC_XC3Q_PLUS,A,B,EPS) 
     3                + ( CONST3 - PLUS_TERM3 ) * XC(Y,Q) )
         ELSE
            INTEG3(4) = 0D0
         ENDIF
         INTEG3(5) = 0D0
         INTEG3(6) = 0D0
*
         NLOG1 = 0D0
         NLOG2 = 0D0
         NLOG3 = 0D0
*
         NLOQ1 = 0D0
         NLOQ2 = 0D0
         NLOQ3 = 0D0
*
         IF(IPT.GE.1)THEN
            NLOG1 = 2D0 * INTEG1(0)
            NLOG2 = 2D0 * INTEG2(0)
            NLOG3 = 2D0 * INTEG3(0)
*
            DO IPDF=1,6
               NLOQ1 = NLOQ1 + INTEG1(IPDF)
               NLOQ2 = NLOQ2 + INTEG2(IPDF)
               NLOQ3 = NLOQ3 + INTEG3(IPDF)
            ENDDO
         ENDIF
*
         NLO1 = 2D0 * ( NLOG1 + NLOQ1 )
         NLO2 = 2D0 * ( NLOG2 + NLOQ2 )
         NLO3 = 2D0 * ( NLOG3 + NLOQ3 )
*
         IF(VFNS.EQ."FFNS")THEN
            IF(HQMASS.EQ.1)THEN
               H1 = 4D0 * ( 4D0 + 3D0 * DLOG(Q2/Q2H) ) / 3D0
*
               YPL = Y * ( 1D0 + DY )
               YMN = Y * ( 1D0 - DY )
*
               LOP = 2D0 * V_CD2 * ( FRAT * XD(YPL,Q)
     1             + ( 1D0 - FRAT ) * XU(YPL,Q) )
     2             + 2D0 * V_CS2 * XS(YPL,Q)
               LOM = 2D0 * V_CD2 * ( FRAT * XD(YMN,Q)
     1             + ( 1D0 - FRAT ) * XU(YMN,Q) )
     2             + 2D0 * V_CS2 * XS(YMN,Q)
*
               DLO = ( LOP - LOM ) / ( YPL - YMN )
*
               NLO1 = NLO1 + 2D0 * H1 * ( 1D0 - LAMBDA ) * Y * DLO
               NLO2 = NLO2 + 2D0 * H1 * ( 1D0 - LAMBDA ) * Y * DLO
               NLO3 = NLO3 + 2D0 * H1 * ( 1D0 - LAMBDA ) 
     1              * ( - LO + Y * DLO )
            ENDIF
         ENDIF
      ENDIF
*
*     NNLO
*
*     Initialize
*
      NNLO1 = 0D0
      NNLO2 = 0D0
      NNLO3 = 0D0
      IF(PTO.GE.2.AND.VFNS.EQ."ZMVN")THEN
*
         IPT = 2
*
*        Gluon + Pure Singlet
*
         INTEG1(0) = ( V_CD2 + V_CS2 ) 
     1             * ( DGAUSS(XGL_XC1G,A,B,EPS) + C2G2C(X,1) * XGL(X,Q)
     2             +   DGAUSS(XPS_XC1Q_UNPLUS,A,B,EPS) )
         INTEG2(0) = ( V_CD2 + V_CS2 ) 
     1             * ( DGAUSS(XGL_XC2G,A,B,EPS) + C2G2C(X,1) * XGL(X,Q)
     2             +   DGAUSS(XPS_XC2Q_UNPLUS,A,B,EPS) )
         INTEG3(0) = 0D0
*
*        Quark
*
*        Contribution coming from the plus prescription in an integral between x and 1
*
         PLUS_TERM1 = C2NN2C(X,NF) - CLNN2C(X,NF)
         PLUS_TERM2 = C2NN2C(X,NF)
         PLUS_TERM3 = C3NM2C(X,NF)
*
*        F1
*
         INTEG1(1) = V_CD2 * FRAT 
     1             * ( DGAUSS(XD_XC1Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XD_XC1Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM1 * XD(Y,Q) )
         INTEG1(2) = V_CD2 * ( 1D0 - FRAT )
     1             * ( DGAUSS(XU_XC1Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XU_XC1Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM1 * XU(Y,Q) )
         INTEG1(3) = V_CS2
     1             * ( DGAUSS(XS_XC1Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XS_XC1Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM1 * XS(Y,Q) )
         INTEG1(4) = ( V_CD2 + V_CS2 ) * FRAT 
     1             * ( DGAUSS(XU_XC1Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XU_XC1Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM1 * XU(Y,Q) )
         INTEG1(5) = 0D0
         INTEG1(6) = 0D0
*
*        F2
*
         INTEG2(1) = V_CD2 * FRAT 
     1             * ( DGAUSS(XD_XC2Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XD_XC2Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM2 * XD(Y,Q) )
         INTEG2(2) = V_CD2 * ( 1D0 - FRAT )
     1             * ( DGAUSS(XU_XC2Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XU_XC2Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM2 * XU(Y,Q) )
         INTEG2(3) = V_CS2
     1             * ( DGAUSS(XS_XC2Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XS_XC2Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM2 * XS(Y,Q) )
         INTEG2(4) = ( V_CD2 + V_CS2 ) * FRAT 
     1             * ( DGAUSS(XU_XC2Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XU_XC2Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM2 * XU(Y,Q) )
         INTEG2(5) = 0D0
         INTEG2(6) = 0D0
*
*        F3
*
         INTEG3(1) = V_CD2 * FRAT 
     1             * ( DGAUSS(XD_XC3Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XD_XC3Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM3 * XD(Y,Q) )
         INTEG3(2) = V_CD2 * ( 1D0 - FRAT )
     1             * ( DGAUSS(XU_XC3Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XU_XC3Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM3 * XU(Y,Q) )
         INTEG3(3) = V_CS2
     1             * ( DGAUSS(XS_XC3Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XS_XC3Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM3 * XS(Y,Q) )
         INTEG3(4) = - ( V_CD2 + V_CS2 ) * FRAT 
     1             * ( DGAUSS(XU_XC3Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XU_XC3Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM3 * XU(Y,Q) )
         INTEG3(5) = 0D0
         INTEG3(6) = 0D0
*
         NNLOG1 = 0D0
         NNLOG2 = 0D0
         NNLOG3 = 0D0
*
         NNLOQ1 = 0D0
         NNLOQ2 = 0D0
         NNLOQ3 = 0D0
*
         IF(IPT.GE.1)THEN
            NNLOG1 = INTEG1(0)
            NNLOG2 = INTEG2(0)
            NNLOG3 = INTEG3(0)
            DO IPDF=1,6
               NNLOQ1 = NNLOQ1 + INTEG1(IPDF)
               NNLOQ2 = NNLOQ2 + INTEG2(IPDF)
               NNLOQ3 = NNLOQ3 + INTEG3(IPDF)
            ENDDO
         ENDIF
*
         NNLO1 = 2D0 * ( NNLOG1 + NNLOQ1 )
         NNLO2 = 2D0 * ( NNLOG2 + NNLOQ2 )
         NNLO3 = 2D0 * ( NNLOG3 + NNLOQ3 )
      ENDIF
*
      F1C = LO1 + ASQ * NLO1 + ASQ2 * NNLO1
      F2C = LO2 + ASQ * NLO2 + ASQ2 * NNLO2
      F3C = F3SGN * ( LO3 + ASQ * NLO3 + ASQ2 * NNLO3 )
*
*     FL
*
      IF(VFNS.EQ."FFNS")THEN
         FLC = F2C - LAMBDA * F1C
      ELSEIF(VFNS.EQ."ZMVN".OR.VFNS.EQ."FFN0")THEN
         FLC = F2C - F1C
      ENDIF
*
*     Evaluation of the reduced cross section
*
      SIGMAC = FACT_CC * ( YPCC * F2C - RAP**2.D0 * FLC - YM * F3C )
*
      IF(WVFNS.EQ."FONLL")THEN
         F1CG(IGM)    = F1C
         F2CG(IGM)    = F2C
         F3CG(IGM)    = F3C
         FLCG(IGM)    = FLC
         SIGMACG(IGM) = SIGMAC
*
         IF(IGM.LE.2)THEN
            IGM = IGM + 1
            GOTO 405
         ENDIF
*
         F1C    = F1CG(1) + DAMP(4) * ( F1CG(2) - F1CG(3) )
         F2C    = F2CG(1) + DAMP(4) * ( F2CG(2) - F2CG(3) )
         F3C    = F3CG(1) + DAMP(4) * ( F3CG(2) - F3CG(3) )
         FLC    = FLCG(1) + DAMP(4) * ( FLCG(2) - FLCG(3) )
         SIGMAC = SIGMACG(1) + DAMP(4) * ( SIGMACG(2) - SIGMACG(3) )
      ENDIF
*
************************************************************************
*     BOTTOM
************************************************************************
*
      IH  = 5
      Q2H = Q2TH(IH)
*
*     Bottom contribution evaluated only for energies over the threshold
*
      IF(Q2.LT.Q2H) GOTO 200
*
*     In case GM scheme is chosen
*
      IGM = 1
 505  IF(WVFNS.EQ."FONLL")THEN
         VFNS = IVFNS(IGM)
      ENDIF
*
      LAMBDA = 1D0 / ( 1D0 + Q2H / Q2 )
      XI     = X / LAMBDA
      IF(XI.GE.1D0) XI = 0.99999999D0
      IF(VFNS.EQ."ZMVN".OR.VFNS.EQ."FFN0")THEN
         Y      = X
         CONST1 = - 4D0 / 3D0 * ( PI**2D0 / 3D0 + 9D0 / 2D0 )
         CONST2 = CONST1
         CONST3 = CONST1
         XBJ    = X
      ELSEIF(VFNS.EQ."FFNS")THEN
         Y      = XI
         KA     = ( 1D0 - LAMBDA ) * DLOG( 1D0 - LAMBDA ) / LAMBDA
         CONST1 = - 4D0 / 3D0 * ( 4D0 + 1D0 / ( 2D0 * LAMBDA )
     1          + PI**2D0 / 3D0 + ( 1D0 + 3D0 * LAMBDA ) * KA 
     2          / (2D0 * LAMBDA ) )
         CONST2 = - 4D0 / 3D0 * ( 4D0 + 1D0 / ( 2D0 * LAMBDA )
     1          + PI**2D0 / 3D0 + ( 1D0 + LAMBDA ) * KA 
     2          / (2D0 * LAMBDA ) )
         CONST3 = CONST1
         XBJ    = XI
      ENDIF
*
*     LO
*
      LO = 2D0 * V_UB2 * ( FRAT * XU(Y,Q) + ( 1D0 - FRAT ) *XD(Y,Q))
*
      LO1 = LO
      LO2 = LO
      LO3 = LO
      IF(VFNS.EQ."ZMVN")THEN
         LO1 = LO1 + 2D0 * ( V_UB2 + V_CB2 ) * XB(Y,Q) 
     1       + 2D0 * V_CB2 * XC(Y,Q)
         LO2 = LO2 + 2D0 * ( V_UB2 + V_CB2 ) * XB(Y,Q) 
     1       + 2D0 * V_CB2 * XC(Y,Q)
         LO3 = LO3 - 2D0 * ( V_UB2 + V_CB2 ) * XB(Y,Q) 
     1       + 2D0 * V_CB2 * XC(Y,Q)
      ENDIF
*
*     Bounds of the integrals
*
      A = Y
      B = 1D0
      C = 0D0
      D = Y
*
*     NLO
*
      NLO1 = 0D0
      NLO2 = 0D0
      NLO3 = 0D0
      IF(PTO.GE.1)THEN
*
         IPT = 1
*
*        Gluon
*
         INTEG1(0) = 2D0 * ( V_UB2 + V_CB2 ) * DGAUSS(XGL_XC1G,A,B,EPS)
         INTEG2(0) = 2D0 * ( V_UB2 + V_CB2 ) * DGAUSS(XGL_XC2G,A,B,EPS)
         INTEG3(0) = 2D0 * ( V_UB2 + V_CB2 ) * DGAUSS(XGL_XC3G,A,B,EPS)
*
*        Quark
*
*        Contribution coming from the plus prescription in an integral between x and 1
*
         PLUS_TERM1 = DGAUSS(C1Q_PLUS_WRAP,C,D,EPS)
         PLUS_TERM2 = DGAUSS(C2Q_PLUS_WRAP,C,D,EPS)
         PLUS_TERM3 = DGAUSS(C3Q_PLUS_WRAP,C,D,EPS)
*
*        F1
*
         INTEG1(1) = 2D0 * V_UB2 * ( 1D0 - FRAT )
     1             * ( DGAUSS(XD_XC1Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XD_XC1Q_PLUS,A,B,EPS) 
     3             + ( CONST1 - PLUS_TERM1 ) * XD(Y,Q) )
         INTEG1(2) = 2D0 * V_UB2 * FRAT
     1             * ( DGAUSS(XU_XC1Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XU_XC1Q_PLUS,A,B,EPS) 
     3             + ( CONST1 - PLUS_TERM1 ) * XU(Y,Q) )
         INTEG1(3) = 0D0
         IF(VFNS.EQ."ZMVN")THEN
            INTEG1(4) = 2D0 * V_CB2 
     1                * ( DGAUSS(XC_XC1Q_UNPLUS,A,B,EPS)
     2                +   DGAUSS(XC_XC1Q_PLUS,A,B,EPS) 
     3                + ( CONST1 - PLUS_TERM1 ) * XC(Y,Q) )
            INTEG1(5) = 2D0 * ( V_UB2 + V_CB2 )
     1                * ( DGAUSS(XB_XC1Q_UNPLUS,A,B,EPS)
     2                +   DGAUSS(XB_XC1Q_PLUS,A,B,EPS)
     3                + ( CONST1 - PLUS_TERM1 ) * XB(Y,Q) )
         ELSE
            INTEG1(4) = 0D0
            INTEG1(5) = 0D0
         ENDIF
         INTEG1(6) = 0D0
*
*        F2
*
         INTEG2(1) = 2D0 * V_UB2 * ( 1D0 - FRAT )
     1             * ( DGAUSS(XD_XC2Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XD_XC2Q_PLUS,A,B,EPS) 
     3             + ( CONST2 - PLUS_TERM2 ) * XD(Y,Q) )
         INTEG2(2) = 2D0 * V_UB2 * FRAT
     1             * ( DGAUSS(XU_XC2Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XU_XC2Q_PLUS,A,B,EPS) 
     3             + ( CONST2 - PLUS_TERM2 ) * XU(Y,Q) )
         INTEG2(3) = 0D0
         IF(VFNS.EQ."ZMVN")THEN
            INTEG2(4) = 2D0 * V_CB2
     1                * ( DGAUSS(XC_XC2Q_UNPLUS,A,B,EPS)
     2                +   DGAUSS(XC_XC2Q_PLUS,A,B,EPS) 
     3                + ( CONST2 - PLUS_TERM2 ) * XC(Y,Q) )
            INTEG2(5) = 2D0 * ( V_UB2 + V_CB2 )
     1                * ( DGAUSS(XB_XC2Q_UNPLUS,A,B,EPS)
     2                +   DGAUSS(XB_XC2Q_PLUS,A,B,EPS) 
     3                + ( CONST2 - PLUS_TERM2 ) * XB(Y,Q) )
         ELSE
            INTEG2(4) = 0D0
            INTEG2(5) = 0D0
         ENDIF
         INTEG2(6) = 0D0
*
*        F3
*
         INTEG3(1) = 2D0 * V_UB2 * ( 1D0 - FRAT )
     1             * ( DGAUSS(XD_XC3Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XD_XC3Q_PLUS,A,B,EPS) 
     3             + ( CONST3 - PLUS_TERM3 ) * XD(Y,Q) )
         INTEG3(2) = 2D0 * V_UB2 * FRAT
     1             * ( DGAUSS(XU_XC3Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XU_XC3Q_PLUS,A,B,EPS) 
     3             + ( CONST3 - PLUS_TERM3 ) * XU(Y,Q) )
         INTEG3(3) = 0D0
         IF(VFNS.EQ."ZMVN")THEN
            INTEG3(4) = 2D0 * V_CB2
     1                * ( DGAUSS(XC_XC3Q_UNPLUS,A,B,EPS)
     2                +   DGAUSS(XC_XC3Q_PLUS,A,B,EPS) 
     3                + ( CONST3 - PLUS_TERM3 ) * XC(Y,Q) )
            INTEG3(5) = - 2D0 * ( V_UB2 + V_CB2 )
     1                * ( DGAUSS(XB_XC3Q_UNPLUS,A,B,EPS)
     2                +   DGAUSS(XB_XC3Q_PLUS,A,B,EPS) 
     3                + ( CONST3 - PLUS_TERM3 ) * XB(Y,Q) )
         ELSE
            INTEG3(4) = 0D0
            INTEG3(5) = 0D0
         ENDIF
         INTEG3(6) = 0D0
*
         NLOG1 = 0D0
         NLOG2 = 0D0
         NLOG3 = 0D0
*
         NLOQ1 = 0D0
         NLOQ2 = 0D0
         NLOQ3 = 0D0
*
         IF(IPT.GE.1)THEN
            NLOG1 = 2D0 * INTEG1(0)
            NLOG2 = 2D0 * INTEG2(0)
            NLOG3 = 2D0 * INTEG3(0)
*
            DO IPDF=1,6
               NLOQ1 = NLOQ1 + INTEG1(IPDF)
               NLOQ2 = NLOQ2 + INTEG2(IPDF)
               NLOQ3 = NLOQ3 + INTEG3(IPDF)
            ENDDO
         ENDIF
*
         NLO1 = 2D0 * ( NLOG1 + NLOQ1 )
         NLO2 = 2D0 * ( NLOG2 + NLOQ2 )
         NLO3 = 2D0 * ( NLOG3 + NLOQ3 )
      ENDIF
*
*     NNLO
*
*     Initialize
*
      NNLO1 = 0D0
      NNLO2 = 0D0
      NNLO3 = 0D0
      IF(PTO.GE.2.AND.VFNS.EQ."ZMVN")THEN
*
         IPT = 2
*
*        Gluon + Pure Singlet
*
         INTEG1(0) = ( V_UB2 + V_CB2 ) 
     1             * ( DGAUSS(XGL_XC1G,A,B,EPS) + C2G2C(X,1) * XGL(X,Q)
     2             +   DGAUSS(XPS_XC1Q_UNPLUS,A,B,EPS) )
         INTEG2(0) = ( V_UB2 + V_CB2 ) 
     1             * ( DGAUSS(XGL_XC2G,A,B,EPS) + C2G2C(X,1) * XGL(X,Q)
     2             +   DGAUSS(XPS_XC2Q_UNPLUS,A,B,EPS) )
         INTEG3(0) = 0D0
*
*        Quark
*
*        Contribution coming from the plus prescription in an integral between x and 1
*
         PLUS_TERM1 = C2NN2C(X,NF) - CLNN2C(X,NF)
         PLUS_TERM2 = C2NN2C(X,NF)
         PLUS_TERM3 = C3NM2C(X,NF)
*
*        F1
*
         INTEG1(1) = V_UB2 * ( 1D0 - FRAT )
     1             * ( DGAUSS(XD_XC1Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XD_XC1Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM1 * XD(Y,Q) )
         INTEG1(2) = V_UB2 * FRAT
     1             * ( DGAUSS(XU_XC1Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XU_XC1Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM1 * XU(Y,Q) )
         INTEG1(3) = 0D0
         INTEG1(4) = V_CB2 
     1             * ( DGAUSS(XC_XC1Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XC_XC1Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM1 * XC(Y,Q) )
         INTEG1(5) = ( V_UB2 + V_CB2 )
     1             * ( DGAUSS(XB_XC1Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XB_XC1Q_PLUS,A,B,EPS)
     3             +   PLUS_TERM1 * XB(Y,Q) )
         INTEG1(6) = 0D0
*
*        F2
*
         INTEG2(1) = V_UB2 * ( 1D0 - FRAT )
     1             * ( DGAUSS(XD_XC2Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XD_XC2Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM2 * XD(Y,Q) )
         INTEG2(2) = V_UB2 * FRAT
     1             * ( DGAUSS(XU_XC2Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XU_XC2Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM2 * XU(Y,Q) )
         INTEG2(3) = 0D0
         INTEG2(4) = V_CB2
     1             * ( DGAUSS(XC_XC2Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XC_XC2Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM2 * XC(Y,Q) )
         INTEG2(5) = ( V_UB2 + V_CB2 )
     1             * ( DGAUSS(XB_XC2Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XB_XC2Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM2 * XB(Y,Q) )
         INTEG2(6) = 0D0
*
*        F3
*
         INTEG3(1) = V_UB2 * ( 1D0 - FRAT )
     1             * ( DGAUSS(XD_XC3Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XD_XC3Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM3 * XD(Y,Q) )
         INTEG3(2) = V_UB2 * FRAT
     1             * ( DGAUSS(XU_XC3Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XU_XC3Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM3 * XU(Y,Q) )
         INTEG3(3) = 0D0
         INTEG3(4) = V_CB2
     1             * ( DGAUSS(XC_XC3Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XC_XC3Q_PLUS,A,B,EPS) 
     3             + ( CONST3 - PLUS_TERM3 ) * XC(Y,Q) )
         INTEG3(5) = - ( V_UB2 + V_CB2 )
     1             * ( DGAUSS(XB_XC3Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XB_XC3Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM3 * XB(Y,Q) )
         INTEG3(6) = 0D0
*
         NNLOG1 = 0D0
         NNLOG2 = 0D0
         NNLOG3 = 0D0
*
         NNLOQ1 = 0D0
         NNLOQ2 = 0D0
         NNLOQ3 = 0D0
*
         IF(IPT.GE.1)THEN
            NNLOG1 = INTEG1(0)
            NNLOG2 = INTEG2(0)
            NNLOG3 = INTEG3(0)
            DO IPDF=1,6
               NNLOQ1 = NNLOQ1 + INTEG1(IPDF)
               NNLOQ2 = NNLOQ2 + INTEG2(IPDF)
               NNLOQ3 = NNLOQ3 + INTEG3(IPDF)
            ENDDO
         ENDIF
*
         NNLO1 = 2D0 * ( NNLOG1 + NNLOQ1 )
         NNLO2 = 2D0 * ( NNLOG2 + NNLOQ2 )
         NNLO3 = 2D0 * ( NNLOG3 + NNLOQ3 )
      ENDIF
*
      F1B = LO1 + ASQ * NLO1 + ASQ2 * NNLO1
      F2B = LO2 + ASQ * NLO2 + ASQ2 * NNLO2
      F3B = - NEUT * F3SGN * ( LO3 + ASQ * NLO3 + ASQ2 * NNLO3 )
*
*     FL
*
      IF(VFNS.EQ."FFNS")THEN
         FLB = F2B - LAMBDA * F1B
      ELSEIF(VFNS.EQ."ZMVN".OR.VFNS.EQ."FFN0")THEN
         FLB = F2B - F1B
      ENDIF
*
*     Evaluation of the reduced cross section
*
      SIGMAB = FACT_CC * ( YPCC * F2B - RAP**2D0 * FLB - YM * F3B )
*
      IF(WVFNS.EQ."FONLL")THEN
         F1BG(IGM)    = F1B
         F2BG(IGM)    = F2B
         F3BG(IGM)    = F3B
         FLBG(IGM)    = FLB
         SIGMABG(IGM) = SIGMAB
*
         IF(IGM.LE.2)THEN
            IGM = IGM + 1
            GOTO 505
         ENDIF
*
         F1B    = F1BG(1) + DAMP(4) * ( F1BG(2) - F1BG(3) )
         F2B    = F2BG(1) + DAMP(4) * ( F2BG(2) - F2BG(3) )
         F3B    = F3BG(1) + DAMP(4) * ( F3BG(2) - F3BG(3) )
         FLB    = FLBG(1) + DAMP(4) * ( FLBG(2) - FLBG(3) )
         SIGMAB = SIGMABG(1) + DAMP(4) * ( SIGMABG(2) - SIGMABG(3) )
      ENDIF
*
************************************************************************
*     TOP
************************************************************************
*
      IH  = 6
      Q2H = Q2TH(IH)
*
*     Top contribution evaluated only for energies over the threshold
*
      IF(Q2.LT.Q2H) GOTO 200
*
*     In case GM scheme is chosen
*
      IGM = 1
 605  IF(WVFNS.EQ."FONLL")THEN
         VFNS = IVFNS(IGM)
      ENDIF
*
      LAMBDA = 1D0 / ( 1D0 + Q2H / Q2 )
      XI     = X / LAMBDA
      IF(XI.GE.1D0) XI = 0.99999999D0
      IF(VFNS.EQ."ZMVN".OR.VFNS.EQ."FFN0")THEN
         Y      = X
         CONST1 = - 4D0 / 3D0 * ( PI**2D0 / 3D0 + 9D0 / 2D0 )
         CONST2 = CONST1
         CONST3 = CONST1
         XBJ    = X
      ELSEIF(VFNS.EQ."FFNS")THEN
         Y      = XI
         KA     = ( 1D0 - LAMBDA ) * DLOG( 1D0 - LAMBDA ) / LAMBDA
         CONST1 = - 4D0 / 3D0 * ( 4D0 + 1D0 / ( 2D0 * LAMBDA )
     1          + PI**2D0 / 3D0 + ( 1D0 + 3D0 * LAMBDA ) * KA 
     2          / (2D0 * LAMBDA ) )
         CONST2 = - 4D0 / 3D0 * ( 4D0 + 1D0 / ( 2D0 * LAMBDA )
     1          + PI**2D0 / 3D0 + ( 1D0 + LAMBDA ) * KA 
     2          / (2D0 * LAMBDA ) )
         CONST3 = CONST1
         XBJ    = XI
      ENDIF
*
*     LO
*
      LO = 2D0 * V_TD2 * ( FRAT * XD(Y,Q) + (1D0-FRAT) * XU(Y,Q) )
     1   + 2D0 * V_TS2 * XS(Y,Q)
*
      LO1 = LO
      LO2 = LO
      LO3 = LO
      IF(VFNS.EQ."ZMVN")THEN
         LO1 = LO1 + 2D0 * ( V_TD2 + V_TS2 + V_TB2 ) * XT(Y,Q)
     1       + 2D0 * V_TB2 * XB(Y,Q)
         LO2 = LO2 + 2D0 * ( V_TD2 + V_TS2 + V_TB2 ) * XT(Y,Q)
     1       + 2D0 * V_TB2 * XB(Y,Q)
         LO3 = LO3 - 2D0 * ( V_TD2 + V_TS2 + V_TB2 ) * XT(Y,Q)
     1       + 2D0 * V_TB2 * XB(Y,Q)
      ENDIF
*
*     Bounds of the integrals
*
      A = Y
      B = 1D0
      C = 0D0
      D = Y
*
*     NLO
*
      NLO1 = 0D0
      NLO2 = 0D0
      NLO3 = 0D0
      IF(PTO.GE.1)THEN
*
         IPT = 1
*
*        Gluon
*
         INTEG1(0) = 2D0 * ( V_TD2 + V_TS2 + V_TB2 ) 
     1                    * DGAUSS(XGL_XC1G,A,B,EPS)
         INTEG2(0) = 2D0 * ( V_TD2 + V_TS2 + V_TB2 ) 
     1                    * DGAUSS(XGL_XC2G,A,B,EPS)
         INTEG3(0) = 2D0 * ( V_TD2 + V_TS2 + V_TB2 ) 
     1                    * DGAUSS(XGL_XC3G,A,B,EPS)
*
*        Quark
*
*        Contribution coming from the plus prescription in an integral between x and 1
*
         PLUS_TERM1 = DGAUSS(C1Q_PLUS_WRAP,C,D,EPS)
         PLUS_TERM2 = DGAUSS(C2Q_PLUS_WRAP,C,D,EPS)
         PLUS_TERM3 = DGAUSS(C3Q_PLUS_WRAP,C,D,EPS)
*
*        F1
*
         INTEG1(1) = 2D0 * V_TD2 * FRAT 
     1             * ( DGAUSS(XD_XC1Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XD_XC1Q_PLUS,A,B,EPS) 
     3             + ( CONST1 - PLUS_TERM1 ) * XD(Y,Q) )
         INTEG1(2) = 2D0 * V_TD2 * ( 1D0 - FRAT )
     1             * ( DGAUSS(XU_XC1Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XU_XC1Q_PLUS,A,B,EPS) 
     3             + ( CONST1 - PLUS_TERM1 ) * XU(Y,Q) )
         INTEG1(3) = 2D0 * V_TS2
     1             * ( DGAUSS(XS_XC1Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XS_XC1Q_PLUS,A,B,EPS) 
     3             + ( CONST1 - PLUS_TERM1 ) * XS(Y,Q) )
         IF(VFNS.EQ."ZMVN")THEN
            INTEG1(4) = 0D0
            INTEG1(5) = 2D0 * V_TB2
     1                * ( DGAUSS(XB_XC1Q_UNPLUS,A,B,EPS)
     2                +   DGAUSS(XB_XC1Q_PLUS,A,B,EPS) 
     3                + ( CONST1 - PLUS_TERM1 ) * XB(Y,Q) )
            INTEG1(6) = 2D0 * ( V_TD2 + V_TS2 + V_TB2 )
     1                * ( DGAUSS(XT_XC1Q_UNPLUS,A,B,EPS)
     2                +   DGAUSS(XT_XC1Q_PLUS,A,B,EPS) 
     3                + ( CONST1 - PLUS_TERM1 ) * XT(Y,Q) )
         ELSE
            INTEG1(4) = 0D0
            INTEG1(5) = 0D0
            INTEG1(6) = 0D0
         ENDIF
*
*        F2
*
         INTEG2(1) = 2D0 * V_TD2 * FRAT 
     1             * ( DGAUSS(XD_XC2Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XD_XC2Q_PLUS,A,B,EPS) 
     3             + ( CONST2 - PLUS_TERM2 ) * XD(Y,Q) )
         INTEG2(2) = 2D0 * V_TD2 * ( 1D0 - FRAT )
     1             * ( DGAUSS(XU_XC2Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XU_XC2Q_PLUS,A,B,EPS) 
     3             + ( CONST2 - PLUS_TERM2 ) * XU(Y,Q) )
         INTEG2(3) = 2D0 * V_TS2
     1             * ( DGAUSS(XS_XC2Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XS_XC2Q_PLUS,A,B,EPS) 
     3             + ( CONST2 - PLUS_TERM2 ) * XS(Y,Q) )
         IF(VFNS.EQ."ZMVN")THEN
            INTEG2(4) = 0D0
            INTEG2(5) = 2D0 * V_TB2
     1                * ( DGAUSS(XB_XC2Q_UNPLUS,A,B,EPS)
     2                +   DGAUSS(XB_XC2Q_PLUS,A,B,EPS) 
     3                + ( CONST2 - PLUS_TERM2 ) * XB(Y,Q) )
            INTEG2(6) = 2D0 * ( V_TD2 + V_TS2 + V_TB2 )
     1                * ( DGAUSS(XT_XC2Q_UNPLUS,A,B,EPS)
     2                +   DGAUSS(XT_XC2Q_PLUS,A,B,EPS) 
     3                + ( CONST2 - PLUS_TERM2 ) * XT(Y,Q) )
         ELSE
            INTEG2(4) = 0D0
            INTEG2(5) = 0D0
            INTEG2(6) = 0D0
         ENDIF
*
*        F3
*
         INTEG3(1) = 2D0 * V_TD2 * FRAT 
     1             * ( DGAUSS(XD_XC3Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XD_XC3Q_PLUS,A,B,EPS) 
     3             + ( CONST3 - PLUS_TERM3 ) * XD(Y,Q) )
         INTEG3(2) = 2D0 * V_TD2 * ( 1D0 - FRAT )
     1             * ( DGAUSS(XU_XC3Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XU_XC3Q_PLUS,A,B,EPS) 
     3             + ( CONST3 - PLUS_TERM3 ) * XU(Y,Q) )
         INTEG3(3) = 2D0 * V_TS2
     1             * ( DGAUSS(XS_XC3Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XS_XC3Q_PLUS,A,B,EPS) 
     3             + ( CONST3 - PLUS_TERM3 ) * XS(Y,Q) )
         IF(VFNS.EQ."ZMVN")THEN
            INTEG3(4) = 0D0
            INTEG3(5) = 2D0 * V_TB2
     1                * ( DGAUSS(XB_XC3Q_UNPLUS,A,B,EPS)
     2                +   DGAUSS(XB_XC3Q_PLUS,A,B,EPS) 
     3                + ( CONST3 - PLUS_TERM3 ) * XB(Y,Q) )
            INTEG3(6) = - 2D0 * ( V_TD2 + V_TS2 + V_TB2 )
     1                * ( DGAUSS(XT_XC3Q_UNPLUS,A,B,EPS)
     2                +   DGAUSS(XT_XC3Q_PLUS,A,B,EPS) 
     3                + ( CONST3 - PLUS_TERM3 ) * XT(Y,Q) )
         ELSE
            INTEG3(4) = 0D0
            INTEG3(5) = 0D0
            INTEG3(6) = 0D0
         ENDIF
*
         NLOG1 = 0D0
         NLOG2 = 0D0
         NLOG3 = 0D0
*
         NLOQ1 = 0D0
         NLOQ2 = 0D0
         NLOQ3 = 0D0
*
         IF(IPT.GE.1)THEN
            NLOG1 = 2D0 * INTEG1(0)
            NLOG2 = 2D0 * INTEG2(0)
            NLOG3 = 2D0 * INTEG3(0)
*
            DO IPDF=1,6
               NLOQ1 = NLOQ1 + INTEG1(IPDF)
               NLOQ2 = NLOQ2 + INTEG2(IPDF)
               NLOQ3 = NLOQ3 + INTEG3(IPDF)
            ENDDO
         ENDIF
*
         NLO1 = 2D0 * ( NLOG1 + NLOQ1 )
         NLO2 = 2D0 * ( NLOG2 + NLOQ2 )
         NLO3 = 2D0 * ( NLOG3 + NLOQ3 )
      ENDIF
*
*     NNLO
*
*     Initialize
*
      NNLO1 = 0D0
      NNLO2 = 0D0
      NNLO3 = 0D0
      IF(PTO.GE.2.AND.VFNS.EQ."ZMVN")THEN
*
         IPT = 2
*
*        Gluon + Pure Singlet
*
         INTEG1(0) = ( V_TD2 + V_TS2 + V_TB2 )
     1             * ( DGAUSS(XGL_XC1G,A,B,EPS) + C2G2C(X,1) * XGL(X,Q)
     2             +   DGAUSS(XPS_XC1Q_UNPLUS,A,B,EPS) )
         INTEG2(0) = ( V_TD2 + V_TS2 + V_TB2 ) 
     1             * ( DGAUSS(XGL_XC2G,A,B,EPS) + C2G2C(X,1) * XGL(X,Q)
     2             +   DGAUSS(XPS_XC2Q_UNPLUS,A,B,EPS) )
         INTEG3(0) = 0D0
*
*        Quark
*
*        Contribution coming from the plus prescription in an integral between x and 1
*
         PLUS_TERM1 = C2NN2C(X,NF) - CLNN2C(X,NF)
         PLUS_TERM2 = C2NN2C(X,NF)
         PLUS_TERM3 = C3NM2C(X,NF)
*
*        F1
*
         INTEG1(1) = V_TD2 * FRAT 
     1             * ( DGAUSS(XD_XC1Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XD_XC1Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM1 * XD(Y,Q) )
         INTEG1(2) = V_TD2 * ( 1D0 - FRAT )
     1             * ( DGAUSS(XU_XC1Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XU_XC1Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM1 * XU(Y,Q) )
         INTEG1(3) = V_TS2
     1             * ( DGAUSS(XS_XC1Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XS_XC1Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM1 * XS(Y,Q) )
         INTEG1(4) = 0D0
         INTEG1(5) = V_TB2
     1             * ( DGAUSS(XB_XC1Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XB_XC1Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM1 * XB(Y,Q) )
         INTEG1(6) = ( V_TD2 + V_TS2 + V_TB2 )
     1             * ( DGAUSS(XT_XC1Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XT_XC1Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM1 * XT(Y,Q) )
*
*        F2
*
         INTEG2(1) = V_TD2 * FRAT 
     1             * ( DGAUSS(XD_XC2Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XD_XC2Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM2 * XD(Y,Q) )
         INTEG2(2) = V_TD2 * ( 1D0 - FRAT )
     1             * ( DGAUSS(XU_XC2Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XU_XC2Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM2 * XU(Y,Q) )
         INTEG2(3) = V_TS2
     1             * ( DGAUSS(XS_XC2Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XS_XC2Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM2 * XS(Y,Q) )
         INTEG2(4) = 0D0
         INTEG2(5) = V_TB2
     1             * ( DGAUSS(XB_XC2Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XB_XC2Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM2 * XB(Y,Q) )
         INTEG2(6) = ( V_TD2 + V_TS2 + V_TB2 )
     1             * ( DGAUSS(XT_XC2Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XT_XC2Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM2 * XT(Y,Q) )
*
*        F3
*
         INTEG3(1) = V_TD2 * FRAT 
     1             * ( DGAUSS(XD_XC3Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XD_XC3Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM3 * XD(Y,Q) )
         INTEG3(2) = V_TD2 * ( 1D0 - FRAT )
     1             * ( DGAUSS(XU_XC3Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XU_XC3Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM3 * XU(Y,Q) )
         INTEG3(3) = V_TS2
     1             * ( DGAUSS(XS_XC3Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XS_XC3Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM3 * XS(Y,Q) )
         INTEG3(4) = 0D0
         INTEG3(5) = V_TB2
     1             * ( DGAUSS(XB_XC3Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XB_XC3Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM3 * XB(Y,Q) )
         INTEG3(6) = - ( V_TD2 + V_TS2 + V_TB2 )
     1             * ( DGAUSS(XT_XC3Q_UNPLUS,A,B,EPS)
     2             +   DGAUSS(XT_XC3Q_PLUS,A,B,EPS) 
     3             +   PLUS_TERM3 * XT(Y,Q) )
*
         NNLOG1 = 0D0
         NNLOG2 = 0D0
         NNLOG3 = 0D0
*
         NNLOQ1 = 0D0
         NNLOQ2 = 0D0
         NNLOQ3 = 0D0
*
         IF(IPT.GE.1)THEN
            NNLOG1 = INTEG1(0)
            NNLOG2 = INTEG2(0)
            NNLOG3 = INTEG3(0)
            DO IPDF=1,6
               NNLOQ1 = NNLOQ1 + INTEG1(IPDF)
               NNLOQ2 = NNLOQ2 + INTEG2(IPDF)
               NNLOQ3 = NNLOQ3 + INTEG3(IPDF)
            ENDDO
         ENDIF
*
         NNLO1 = 2D0 * ( NNLOG1 + NNLOQ1 )
         NNLO2 = 2D0 * ( NNLOG2 + NNLOQ2 )
         NNLO3 = 2D0 * ( NNLOG3 + NNLOQ3 )
      ENDIF
*
      F1T = LO1 + ASQ * NLO1 + ASQ2 * NNLO1
      F2T = LO2 + ASQ * NLO2 + ASQ2 * NNLO2
      F3T = F3SGN * ( LO3 + ASQ * NLO3 + ASQ2 * NNLO3 )
*
*     FL
*
      IF(VFNS.EQ."FFNS")THEN
         FLT = F2T - LAMBDA * F1T
      ELSEIF(VFNS.EQ."ZMVN".OR.VFNS.EQ."FFN0")THEN
         FLT = F2T - F1T
      ENDIF
*
*     Evaluation of the reduced cross section
*
      SIGMAT = FACT_CC * ( YPCC * F2T - RAP**2D0 * FLT - YM * F3T )
*
      IF(WVFNS.EQ."FONLL")THEN
         F1TG(IGM)    = F1T
         F2TG(IGM)    = F2T
         F3TG(IGM)    = F3T
         FLTG(IGM)    = FLT
         SIGMATG(IGM) = SIGMAT
*
         IF(IGM.LE.2)THEN
            IGM = IGM + 1
            GOTO 605
         ENDIF
*
         F1T    = F1TG(1) + DAMP(4) * ( F1TG(2) - F1TG(3) )
         F2T    = F2TG(1) + DAMP(4) * ( F2TG(2) - F2TG(3) )
         F3T    = F3TG(1) + DAMP(4) * ( F3TG(2) - F3TG(3) )
         FLT    = FLTG(1) + DAMP(4) * ( FLTG(2) - FLTG(3) )
         SIGMAT = SIGMATG(1) + DAMP(4) * ( SIGMATG(2) - SIGMATG(3) )
      ENDIF
*
*     Gather pieces
*
 200  IF(Q2.GT.Q2TH(6))THEN
         F2P = F2L + F2C + F2B + F2T
         F3P = F3L + F3C + F3B + F3T
         FLP = FLL + FLC + FLB + FLT
         SIGMAP = SIGMAL + SIGMAC + SIGMAB + SIGMAT
      ELSEIF(Q2.GT.Q2TH(5))THEN
         F2P = F2L + F2C + F2B
         F3P = F3L + F3C + F3B
         FLP = FLL + FLC + FLB
         SIGMAP = SIGMAL + SIGMAC + SIGMAB
      ELSEIF(Q2.GT.Q2TH(4))THEN
         F2P = F2L + F2C
         F3P = F3L + F3C
         FLP = FLL + FLC
         SIGMAP = SIGMAL + SIGMAC
      ELSEIF(Q2.LT.Q2TH(4))THEN
         F2P = F2L
         F3P = F3L
         FLP = FLL
         SIGMAP = SIGMAL
      ENDIF
*
*     Put results into vectors
*
      F2(3) = F2L
      F2(4) = F2C
      F2(5) = F2B
      F2(6) = F2T
      F2(7) = F2P
*
      F3(3) = F3L
      F3(4) = F3C
      F3(5) = F3B
      F3(6) = F3T
      F3(7) = F3P
*
      FL(3) = FLL
      FL(4) = FLC
      FL(5) = FLB
      FL(6) = FLT
      FL(7) = FLP
*
*     ... include polarization
*
      POLFACT = ( 1D0 + NEUT * POL )
      SIGMA(3) = POLFACT * SIGMAL
      SIGMA(4) = POLFACT * SIGMAC
      SIGMA(5) = POLFACT * SIGMAB
      SIGMA(6) = POLFACT * SIGMAT
      SIGMA(7) = POLFACT * SIGMAP
*
      RETURN
      END
