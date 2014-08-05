************************************************************************
*
*     nc_dis.f:
*
*     x-space computation of NC observables.
*
************************************************************************
      SUBROUTINE NC_DIS(X,QI,QF,Y,POL,PROC,SCHEME,PTO,PDFSET,IREP,
     1                  TARGET,PROJ,F2,F3,FL,SIGMA)
*
      IMPLICIT NONE
*
      include "../commons/hmass.h"
      include "../commons/nucleon.h"
      include "../commons/vfns.h"
      include "../commons/heavythr.h"
      include "../commons/asq.h"
      include "../commons/xbj.h"
      include "../commons/scale.h"
      include "../commons/beta.h"
      include "../commons/jpt.h"
      include "../commons/nf.h"
**
*     Input Varibles
*
      INTEGER PTO,IREP
      DOUBLE PRECISION X,QF,QI,Y,POL
      CHARACTER*2  PROC
      CHARACTER*5  SCHEME
      CHARACTER*53 PDFSET
      CHARACTER*9  TARGET
      CHARACTER*12 PROJ
**
*     Internal Variables
*
      INTEGER IPDF,IE
      INTEGER I,IGM,INC
      DOUBLE PRECISION Q2,W2,DAMP(4:6)
      DOUBLE PRECISION YP,YM
      DOUBLE PRECISION EQ(6),EQ2(6),VQ(6),AQ(6),VE,AE
      DOUBLE PRECISION PZ
      DOUBLE PRECISION BQ(6),DQ(6)
      DOUBLE PRECISION XDP,XUP,XSP,XCP,XBP,XTP
      DOUBLE PRECISION XDM,XUM,XSM,XCM,XBM,XTM
      DOUBLE PRECISION XGL
      DOUBLE PRECISION C2G2C,C2NN2C
      DOUBLE PRECISION CLNN2C
      DOUBLE PRECISION C3NM2C
      DOUBLE PRECISION INTEG2(3),INTEG3(3),INTEGL(3)
      DOUBLE PRECISION CONST
      DOUBLE PRECISION DGAUSS,A,B,C,D,EPS,EPS2
      DOUBLE PRECISION PLUS_TERM2,PLUS_TERM3
      DOUBLE PRECISION LOCAL_TERM2,LOCAL_TERM3,LOCAL_TERML
      DOUBLE PRECISION ADSR,ADLERSR
      DOUBLE PRECISION ASQ2
      DOUBLE PRECISION PI,SW,MZ
      DOUBLE PRECISION LO2L,LO3L,LOLL
      DOUBLE PRECISION LO2C,LO3C,LOLC
      DOUBLE PRECISION LO2B,LO3B,LOLB
      DOUBLE PRECISION LO2T,LO3T,LOLT
      DOUBLE PRECISION NLO2L,NLOG2L,NLOQ2L
      DOUBLE PRECISION NLO2C,NLOG2C,NLOQ2C
      DOUBLE PRECISION NLO2B,NLOG2B,NLOQ2B
      DOUBLE PRECISION NLO2T,NLOG2T,NLOQ2T
      DOUBLE PRECISION NLO3L,NLOG3L,NLOQ3L
      DOUBLE PRECISION NLO3C,NLOG3C,NLOQ3C
      DOUBLE PRECISION NLO3B,NLOG3B,NLOQ3B
      DOUBLE PRECISION NLO3T,NLOG3T,NLOQ3T
      DOUBLE PRECISION NLOLL,NLOGLL,NLOQLL
      DOUBLE PRECISION NLOLC,NLOGLC,NLOQLC
      DOUBLE PRECISION NLOLB,NLOGLB,NLOQLB
      DOUBLE PRECISION NLOLT,NLOGLT,NLOQLT
      DOUBLE PRECISION NNLO2L,NNLOG2L,NNLOQ2L
      DOUBLE PRECISION NNLO2C,NNLOG2C,NNLOQ2C
      DOUBLE PRECISION NNLO2B,NNLOG2B,NNLOQ2B
      DOUBLE PRECISION NNLO2T,NNLOG2T,NNLOQ2T
      DOUBLE PRECISION NNLO3L,NNLOG3L,NNLOQ3L
      DOUBLE PRECISION NNLO3C,NNLOG3C,NNLOQ3C
      DOUBLE PRECISION NNLO3B,NNLOG3B,NNLOQ3B
      DOUBLE PRECISION NNLO3T,NNLOG3T,NNLOQ3T
      DOUBLE PRECISION NNLOLL,NNLOGLL,NNLOQLL
      DOUBLE PRECISION NNLOLC,NNLOGLC,NNLOQLC
      DOUBLE PRECISION NNLOLB,NNLOGLB,NNLOQLB
      DOUBLE PRECISION NNLOLT,NNLOGLT,NNLOQLT
      DOUBLE PRECISION F2CG(3),F3CG(3),FLCG(3)
      DOUBLE PRECISION F2BG(3),F3BG(3),FLBG(3)
      DOUBLE PRECISION F2TG(3),F3TG(3),FLTG(3)
      DOUBLE PRECISION F2P,F3P,FLP,SIGMAP
      DOUBLE PRECISION F2L,F3L,FLL,SIGMAL
      DOUBLE PRECISION F2C,F3C,FLC,SIGMAC
      DOUBLE PRECISION F2B,F3B,FLB,SIGMAB
      DOUBLE PRECISION F2T,F3T,FLT,SIGMAT

      DOUBLE PRECISION XGL_XC2G_NC
      EXTERNAL XGL_XC2G_NC
      DOUBLE PRECISION XPS_XC2Q_UNPLUS_NC
      EXTERNAL XPS_XC2Q_UNPLUS_NC
      DOUBLE PRECISION XDP_XC2Q_UNPLUS_NC
      EXTERNAL XDP_XC2Q_UNPLUS_NC
      DOUBLE PRECISION XDP_XC2Q_PLUS_NC
      EXTERNAL XDP_XC2Q_PLUS_NC
      DOUBLE PRECISION XUP_XC2Q_UNPLUS_NC
      EXTERNAL XUP_XC2Q_UNPLUS_NC
      DOUBLE PRECISION XUP_XC2Q_PLUS_NC
      EXTERNAL XUP_XC2Q_PLUS_NC
      DOUBLE PRECISION XSP_XC2Q_UNPLUS_NC
      EXTERNAL XSP_XC2Q_UNPLUS_NC
      DOUBLE PRECISION XSP_XC2Q_PLUS_NC
      EXTERNAL XSP_XC2Q_PLUS_NC
      DOUBLE PRECISION XCP_XC2Q_UNPLUS_NC
      EXTERNAL XCP_XC2Q_UNPLUS_NC
      DOUBLE PRECISION XCP_XC2Q_PLUS_NC
      EXTERNAL XCP_XC2Q_PLUS_NC
      DOUBLE PRECISION XBP_XC2Q_UNPLUS_NC
      EXTERNAL XBP_XC2Q_UNPLUS_NC
      DOUBLE PRECISION XBP_XC2Q_PLUS_NC
      EXTERNAL XBP_XC2Q_PLUS_NC
      DOUBLE PRECISION XTP_XC2Q_UNPLUS_NC
      EXTERNAL XTP_XC2Q_UNPLUS_NC
      DOUBLE PRECISION XTP_XC2Q_PLUS_NC
      EXTERNAL XTP_XC2Q_PLUS_NC

      DOUBLE PRECISION XDM_XC2Q_UNPLUS_NC
      EXTERNAL XDM_XC2Q_UNPLUS_NC
      DOUBLE PRECISION XDM_XC2Q_PLUS_NC
      EXTERNAL XDM_XC2Q_PLUS_NC
      DOUBLE PRECISION XUM_XC2Q_UNPLUS_NC
      EXTERNAL XUM_XC2Q_UNPLUS_NC
      DOUBLE PRECISION XUM_XC2Q_PLUS_NC
      EXTERNAL XUM_XC2Q_PLUS_NC
      DOUBLE PRECISION XSM_XC2Q_UNPLUS_NC
      EXTERNAL XSM_XC2Q_UNPLUS_NC
      DOUBLE PRECISION XSM_XC2Q_PLUS_NC
      EXTERNAL XSM_XC2Q_PLUS_NC
      DOUBLE PRECISION XCM_XC2Q_UNPLUS_NC
      EXTERNAL XCM_XC2Q_UNPLUS_NC
      DOUBLE PRECISION XCM_XC2Q_PLUS_NC
      EXTERNAL XCM_XC2Q_PLUS_NC
      DOUBLE PRECISION XBM_XC2Q_UNPLUS_NC
      EXTERNAL XBM_XC2Q_UNPLUS_NC
      DOUBLE PRECISION XBM_XC2Q_PLUS_NC
      EXTERNAL XBM_XC2Q_PLUS_NC
      DOUBLE PRECISION XTM_XC2Q_UNPLUS_NC
      EXTERNAL XTM_XC2Q_UNPLUS_NC
      DOUBLE PRECISION XTM_XC2Q_PLUS_NC
      EXTERNAL XTM_XC2Q_PLUS_NC

      DOUBLE PRECISION C2Q_PLUS_NC_WRAP
      EXTERNAL C2Q_PLUS_NC_WRAP

      DOUBLE PRECISION XGL_XC3G_NC
      EXTERNAL XGL_XC3G_NC
      DOUBLE PRECISION XDP_XC3Q_UNPLUS_NC
      EXTERNAL XDP_XC3Q_UNPLUS_NC
      DOUBLE PRECISION XDP_XC3Q_PLUS_NC
      EXTERNAL XDP_XC3Q_PLUS_NC
      DOUBLE PRECISION XUP_XC3Q_UNPLUS_NC
      EXTERNAL XUP_XC3Q_UNPLUS_NC
      DOUBLE PRECISION XUP_XC3Q_PLUS_NC
      EXTERNAL XUP_XC3Q_PLUS_NC
      DOUBLE PRECISION XSP_XC3Q_UNPLUS_NC
      EXTERNAL XSP_XC3Q_UNPLUS_NC
      DOUBLE PRECISION XSP_XC3Q_PLUS_NC
      EXTERNAL XSP_XC3Q_PLUS_NC
      DOUBLE PRECISION XCP_XC3Q_UNPLUS_NC
      EXTERNAL XCP_XC3Q_UNPLUS_NC
      DOUBLE PRECISION XCP_XC3Q_PLUS_NC
      EXTERNAL XCP_XC3Q_PLUS_NC
      DOUBLE PRECISION XBP_XC3Q_UNPLUS_NC
      EXTERNAL XBP_XC3Q_UNPLUS_NC
      DOUBLE PRECISION XBP_XC3Q_PLUS_NC
      EXTERNAL XBP_XC3Q_PLUS_NC
      DOUBLE PRECISION XTP_XC3Q_UNPLUS_NC
      EXTERNAL XTP_XC3Q_UNPLUS_NC
      DOUBLE PRECISION XTP_XC3Q_PLUS_NC
      EXTERNAL XTP_XC3Q_PLUS_NC

      DOUBLE PRECISION XDM_XC3Q_UNPLUS_NC
      EXTERNAL XDM_XC3Q_UNPLUS_NC
      DOUBLE PRECISION XDM_XC3Q_PLUS_NC
      EXTERNAL XDM_XC3Q_PLUS_NC
      DOUBLE PRECISION XUM_XC3Q_UNPLUS_NC
      EXTERNAL XUM_XC3Q_UNPLUS_NC
      DOUBLE PRECISION XUM_XC3Q_PLUS_NC
      EXTERNAL XUM_XC3Q_PLUS_NC
      DOUBLE PRECISION XSM_XC3Q_UNPLUS_NC
      EXTERNAL XSM_XC3Q_UNPLUS_NC
      DOUBLE PRECISION XSM_XC3Q_PLUS_NC
      EXTERNAL XSM_XC3Q_PLUS_NC
      DOUBLE PRECISION XCM_XC3Q_UNPLUS_NC
      EXTERNAL XCM_XC3Q_UNPLUS_NC
      DOUBLE PRECISION XCM_XC3Q_PLUS_NC
      EXTERNAL XCM_XC3Q_PLUS_NC
      DOUBLE PRECISION XBM_XC3Q_UNPLUS_NC
      EXTERNAL XBM_XC3Q_UNPLUS_NC
      DOUBLE PRECISION XBM_XC3Q_PLUS_NC
      EXTERNAL XBM_XC3Q_PLUS_NC
      DOUBLE PRECISION XTM_XC3Q_UNPLUS_NC
      EXTERNAL XTM_XC3Q_UNPLUS_NC
      DOUBLE PRECISION XTM_XC3Q_PLUS_NC
      EXTERNAL XTM_XC3Q_PLUS_NC

      DOUBLE PRECISION C3Q_PLUS_NC_WRAP
      EXTERNAL C3Q_PLUS_NC_WRAP

      DOUBLE PRECISION XGL_XCLG_NC
      EXTERNAL XGL_XCLG_NC
      DOUBLE PRECISION XPS_XCLQ_UNPLUS_NC
      EXTERNAL XPS_XCLQ_UNPLUS_NC
      DOUBLE PRECISION XDP_XCLQ_UNPLUS_NC
      EXTERNAL XDP_XCLQ_UNPLUS_NC
      DOUBLE PRECISION XDP_XCLQ_PLUS_NC
      EXTERNAL XDP_XCLQ_PLUS_NC
      DOUBLE PRECISION XUP_XCLQ_UNPLUS_NC
      EXTERNAL XUP_XCLQ_UNPLUS_NC
      DOUBLE PRECISION XUP_XCLQ_PLUS_NC
      EXTERNAL XUP_XCLQ_PLUS_NC
      DOUBLE PRECISION XSP_XCLQ_UNPLUS_NC
      EXTERNAL XSP_XCLQ_UNPLUS_NC
      DOUBLE PRECISION XSP_XCLQ_PLUS_NC
      EXTERNAL XSP_XCLQ_PLUS_NC
      DOUBLE PRECISION XCP_XCLQ_UNPLUS_NC
      EXTERNAL XCP_XCLQ_UNPLUS_NC
      DOUBLE PRECISION XCP_XCLQ_PLUS_NC
      EXTERNAL XCP_XCLQ_PLUS_NC
      DOUBLE PRECISION XBP_XCLQ_UNPLUS_NC
      EXTERNAL XBP_XCLQ_UNPLUS_NC
      DOUBLE PRECISION XBP_XCLQ_PLUS_NC
      EXTERNAL XBP_XCLQ_PLUS_NC
      DOUBLE PRECISION XTP_XCLQ_UNPLUS_NC
      EXTERNAL XTP_XCLQ_UNPLUS_NC
      DOUBLE PRECISION XTP_XCLQ_PLUS_NC
      EXTERNAL XTP_XCLQ_PLUS_NC

      DOUBLE PRECISION XDM_XCLQ_UNPLUS_NC
      EXTERNAL XDM_XCLQ_UNPLUS_NC
      DOUBLE PRECISION XDM_XCLQ_PLUS_NC
      EXTERNAL XDM_XCLQ_PLUS_NC
      DOUBLE PRECISION XUM_XCLQ_UNPLUS_NC
      EXTERNAL XUM_XCLQ_UNPLUS_NC
      DOUBLE PRECISION XUM_XCLQ_PLUS_NC
      EXTERNAL XUM_XCLQ_PLUS_NC
      DOUBLE PRECISION XSM_XCLQ_UNPLUS_NC
      EXTERNAL XSM_XCLQ_UNPLUS_NC
      DOUBLE PRECISION XSM_XCLQ_PLUS_NC
      EXTERNAL XSM_XCLQ_PLUS_NC
      DOUBLE PRECISION XCM_XCLQ_UNPLUS_NC
      EXTERNAL XCM_XCLQ_UNPLUS_NC
      DOUBLE PRECISION XCM_XCLQ_PLUS_NC
      EXTERNAL XCM_XCLQ_PLUS_NC
      DOUBLE PRECISION XBM_XCLQ_UNPLUS_NC
      EXTERNAL XBM_XCLQ_UNPLUS_NC
      DOUBLE PRECISION XBM_XCLQ_PLUS_NC
      EXTERNAL XBM_XCLQ_PLUS_NC
      DOUBLE PRECISION XTM_XCLQ_UNPLUS_NC
      EXTERNAL XTM_XCLQ_UNPLUS_NC
      DOUBLE PRECISION XTM_XCLQ_PLUS_NC
      EXTERNAL XTM_XCLQ_PLUS_NC

      DOUBLE PRECISION CLQ_PLUS_NC_WRAP
      EXTERNAL CLQ_PLUS_NC_WRAP

      CHARACTER*5 IVFNS(3)
      CHARACTER*5 WVFNS,VFNS_TMP

      PARAMETER(EPS = 1D-5)
      PARAMETER(EPS2= 1D-10)
      PARAMETER(PI  = 3.14159265358979D0)
      PARAMETER(SW  = 0.2312D0)
      PARAMETER(MZ  = 91.1876D0)
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
*
         B1(I) = BETA1(I) / BETA0(I)
         B2(I) = BETA2(I) / BETA0(I)
      ENDDO
*
*     Couplings
*
*     Electric Charges
*
      EQ(1) = - 1D0 / 3D0
      EQ(2) = 2D0 / 3D0
      EQ(3) = - 1D0 / 3D0
      EQ(4) = 2D0 / 3D0
      EQ(5) = - 1D0 / 3D0
      EQ(6) = 2D0 / 3D0
*
*     Squared Charges
*
      EQ2(1) = 1D0 / 9D0
      EQ2(2) = 4D0 / 9D0
      EQ2(3) = 1D0 / 9D0
      EQ2(4) = 4D0 / 9D0
      EQ2(5) = 1D0 / 9D0
      EQ2(6) = 4D0 / 9D0
*
*     Vector Couplings
*
      VQ(1) = - 0.5D0 + 2D0 / 3D0 * SW
      VQ(2) = + 0.5D0 - 4D0 / 3D0 * SW
      VQ(3) = - 0.5D0 + 2D0 / 3D0 * SW
      VQ(4) = + 0.5D0 - 4D0 / 3D0 * SW
      VQ(5) = - 0.5D0 + 2D0 / 3D0 * SW
      VQ(6) = + 0.5D0 - 4D0 / 3D0 * SW
*
*     Axial Couplings
*
      AQ(1) = - 0.5D0
      AQ(2) = + 0.5D0
      AQ(3) = - 0.5D0
      AQ(4) = + 0.5D0
      AQ(5) = - 0.5D0
      AQ(6) = + 0.5D0
*
*     Vector and Axial Electron Couplings
*
      VE = - 0.5D0 + 2D0 * SW
      AE = - 0.5D0
*
*     Subschemes of the GM-VFNS
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
      W2 = Q2 * ( 1D0 - X ) / X
      YM = 1D0 - ( 1D0 - Y )**2D0
      YP = 1D0 + ( 1D0 - Y )**2D0
*
*     Define Observable and check the consistency of the
*     input parameters.
*
      IF(TARGET(1:6).EQ."PROTON")THEN
         INUCL = 0
      ELSEIF(TARGET(1:7).EQ."NEUTRON")THEN
         INUCL = 1
      ELSEIF(TARGET(1:9).EQ."ISOSCALAR")THEN
         INUCL = 2
      ELSE
         WRITE(6,*) "In nc_dis.f:"
         WRITE(6,*) "Unknown value target = ",TARGET
         CALL EXIT(-10)
      ENDIF
*
      IF(PROJ(1:8).EQ."ELECTRON")THEN
         IE  = - 1
      ELSEIF(PROJ(1:8).EQ."POSITRON")THEN
         IE  = 1
      ELSE
         WRITE(6,*) "In nc_dis.f:"
         WRITE(6,*) "Unknown value projectile = ",PROJ
         CALL EXIT(-10)
      ENDIF
*
      IF(PROC(1:2).EQ."EM")THEN
         INC = 0
         DO I=1,6
            BQ(I) = EQ2(I) 
            DQ(I) = 0.0D0
         ENDDO
      ELSEIF(PROC(1:2).EQ."NC")THEN
         INC = 1
         PZ = Q2 / ( Q2 + MZ**2D0 ) / ( 4D0 * SW * ( 1D0 - SW ) )
         DO I=1,6
            BQ(I) = EQ2(I) 
     1            - 2D0 * EQ(I) * VQ(I) * ( VE + IE * POL * AE ) * PZ
     2            + ( VE**2D0 + AE**2D0 ) * ( VQ(I)**2D0 + AQ(I)**2D0 
     3            + IE * POL * 2D0 * VE * AE ) * PZ**2D0 
            DQ(I) = - 2D0 * EQ(I) * AQ(I) * ( AE + IE * POL * VE ) * PZ
     1            + 2D0 * VQ(I) * AQ(I) * ( 2D0 * VE * AE 
     2            + IE * POL * ( VE**2D0 + AE**2D0 ) ) * PZ**2D0
         ENDDO
      ELSE
         WRITE(6,*) "In nc_dis.f:"
         WRITE(6,*) "Unknown value process = ",PROC
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
      IF(PTO.NE.0.AND.PTO.NE.1.AND.PTO.NE.2)THEN
         WRITE(6,*) "In nc_dis.f:"
         WRITE(6,*) "Invalid perturbative order: PTO =",PTO
         CALL EXIT(-10)
      ENDIF
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
*     Initialize MINIMAX tables (FFNS coefficient functions)
*
      CALL INITFFNSTABLES
*
************************************************************************
*
*     Light Part
*
************************************************************************
*
      VFNS_TMP = VFNS
      VFNS     = "ZMVN"
      XBJ      = X
*
*     LO
*
      LO2L = BQ(1) * XDP(X,Q) + BQ(2) * XUP(X,Q) + BQ(3) * XSP(X,Q) 
      LO3L = DQ(1) * XDM(X,Q) + DQ(2) * XUM(X,Q) + DQ(3) * XSM(X,Q) 
      LOLL = 0D0
*
*     Integration bounds
*
      A = X
      B = 1D0
      C = 0D0
      D = X
*
*     NLO
*
*     Initialize
*
      NLO2L = 0D0
      NLO3L = 0D0
      NLOLL = 0D0
*
      IF(PTO.GE.1)THEN
*
         IPT = 1
*
*        Gluon
*
         NLOG2L = 2D0 * ( BQ(1) + BQ(2) + BQ(3) ) 
     1          * DGAUSS(XGL_XC2G_NC,A,B,EPS)
         NLOG3L = 2D0 * ( DQ(1) + DQ(2) + DQ(3) ) 
     1          * DGAUSS(XGL_XC3G_NC,A,B,EPS)
         NLOGLL = 2D0 * ( BQ(1) + BQ(2) + BQ(3) ) 
     1          * DGAUSS(XGL_XCLG_NC,A,B,EPS)
*
*        Quark
*
*        Contribution coming from the plus prescription in an integral between x and 1
*
         CONST = - 8D0 / 3D0 * ( PI**2D0 / 3D0 + 9D0 / 2D0 )
*
         PLUS_TERM2 = DGAUSS(C2Q_PLUS_NC_WRAP,C,D,EPS)
         PLUS_TERM3 = DGAUSS(C3Q_PLUS_NC_WRAP,C,D,EPS)
*
*        F2
*
         INTEG2(1) = BQ(1) * ( DGAUSS(XDP_XC2Q_UNPLUS_NC,A,B,EPS)
     1             + DGAUSS(XDP_XC2Q_PLUS_NC,A,B,EPS) 
     2             + ( CONST - PLUS_TERM2 ) * XDP(X,Q) )
         INTEG2(2) = BQ(2) * ( DGAUSS(XUP_XC2Q_UNPLUS_NC,A,B,EPS)
     1             + DGAUSS(XUP_XC2Q_PLUS_NC,A,B,EPS) 
     2             + ( CONST - PLUS_TERM2 ) * XUP(X,Q) )
         INTEG2(3) = BQ(3) * ( DGAUSS(XSP_XC2Q_UNPLUS_NC,A,B,EPS)
     1             + DGAUSS(XSP_XC2Q_PLUS_NC,A,B,EPS) 
     2             + ( CONST - PLUS_TERM2 ) * XSP(X,Q) )
*
*        F3
*
         INTEG3(1) = DQ(1) * ( DGAUSS(XDM_XC3Q_UNPLUS_NC,A,B,EPS)
     1             + DGAUSS(XDM_XC3Q_PLUS_NC,A,B,EPS) 
     2             + ( CONST - PLUS_TERM3 ) * XDM(X,Q) )
         INTEG3(2) = DQ(2) * ( DGAUSS(XUM_XC3Q_UNPLUS_NC,A,B,EPS)
     1             + DGAUSS(XUM_XC3Q_PLUS_NC,A,B,EPS) 
     2             + ( CONST - PLUS_TERM3 ) * XUM(X,Q) )
         INTEG3(3) = DQ(3) * ( DGAUSS(XSM_XC3Q_UNPLUS_NC,A,B,EPS)
     1             + DGAUSS(XSM_XC3Q_PLUS_NC,A,B,EPS) 
     2             + ( CONST - PLUS_TERM3 ) * XSM(X,Q) )
*
*        FL
*
         INTEGL(1) = BQ(1) * DGAUSS(XDP_XCLQ_UNPLUS_NC,A,B,EPS)
         INTEGL(2) = BQ(2) * DGAUSS(XUP_XCLQ_UNPLUS_NC,A,B,EPS)
         INTEGL(3) = BQ(3) * DGAUSS(XSP_XCLQ_UNPLUS_NC,A,B,EPS)
*
         NLOQ2L = 0D0
         NLOQ3L = 0D0
         NLOQLL = 0D0
         DO IPDF=1,3
            NLOQ2L = NLOQ2L + INTEG2(IPDF)
            NLOQ3L = NLOQ3L + INTEG3(IPDF)
            NLOQLL = NLOQLL + INTEGL(IPDF)
         ENDDO
*
         NLO2L = NLOG2L + NLOQ2L
         NLO3L = NLOG3L + NLOQ3L
         NLOLL = NLOGLL + NLOQLL
      ENDIF
*
*     NNLO
*
*     Initialize
*
      NNLO2L = 0D0
      NNLO3L = 0D0
      NNLOLL = 0D0
*
      IF(PTO.GE.2)THEN
*
         IPT = 2
*
*        Gluon + Pure Singlet
*
         NNLOG2L = ( BQ(1) + BQ(2) + BQ(3) ) 
     1           * ( DGAUSS(XGL_XC2G_NC,A,B,EPS) + C2G2C(X,1) * XGL(X,Q)
     2           +   DGAUSS(XPS_XC2Q_UNPLUS_NC,A,B,EPS) )
         NNLOG3L = ( DQ(1) + DQ(2) + DQ(3) ) 
     1           * DGAUSS(XGL_XC3G_NC,A,B,EPS)
         NNLOGLL = ( BQ(1) + BQ(2) + BQ(3) )
     1           * ( DGAUSS(XGL_XCLG_NC,A,B,EPS)
     2           +   DGAUSS(XPS_XCLQ_UNPLUS_NC,A,B,EPS) )
*
*        Quark (Non-singlet)
*
*        F2
*
         LOCAL_TERM2 = C2NN2C(X,NF)
         INTEG2(1)   = BQ(1) * ( DGAUSS(XDP_XC2Q_UNPLUS_NC,A,B,EPS)
     1               + DGAUSS(XDP_XC2Q_PLUS_NC,A,B,EPS) 
     2               + LOCAL_TERM2 * XDP(X,Q) )
         INTEG2(2)   = BQ(2) * ( DGAUSS(XUP_XC2Q_UNPLUS_NC,A,B,EPS)
     1               + DGAUSS(XUP_XC2Q_PLUS_NC,A,B,EPS) 
     2               + LOCAL_TERM2 * XUP(X,Q) )
         INTEG2(3)   = BQ(3) * ( DGAUSS(XSP_XC2Q_UNPLUS_NC,A,B,EPS)
     1               + DGAUSS(XSP_XC2Q_PLUS_NC,A,B,EPS) 
     2               + LOCAL_TERM2 * XSP(X,Q) )
*
*        F3
*
         LOCAL_TERM3 = C3NM2C(X,NF)
         INTEG3(1)   = DQ(1) * ( DGAUSS(XDM_XC3Q_UNPLUS_NC,A,B,EPS)
     1               + DGAUSS(XDM_XC3Q_PLUS_NC,A,B,EPS) 
     2               + LOCAL_TERM3 * XDM(X,Q) )
         INTEG3(2)   = DQ(2) * ( DGAUSS(XUM_XC3Q_UNPLUS_NC,A,B,EPS)
     1               + DGAUSS(XUM_XC3Q_PLUS_NC,A,B,EPS) 
     2               + LOCAL_TERM3 * XUM(X,Q) )
         INTEG3(3)   = DQ(3) * ( DGAUSS(XSM_XC3Q_UNPLUS_NC,A,B,EPS)
     1               + DGAUSS(XSM_XC3Q_PLUS_NC,A,B,EPS) 
     2               + LOCAL_TERM3 * XSM(X,Q) )
*
*        FL
*
         LOCAL_TERML = CLNN2C(X,NF)
         INTEGL(1)   = BQ(1) * ( DGAUSS(XDP_XCLQ_UNPLUS_NC,A,B,EPS)
     1               + LOCAL_TERML * XDP(X,Q) )
         INTEGL(2)   = BQ(2) * ( DGAUSS(XUP_XCLQ_UNPLUS_NC,A,B,EPS)
     1               + LOCAL_TERML * XUP(X,Q) )
         INTEGL(3)   = BQ(3) * ( DGAUSS(XSP_XCLQ_UNPLUS_NC,A,B,EPS)
     1               + LOCAL_TERML * XSP(X,Q) )
*
*     In case one of the FFN schemes or the GM scheme at NNLO are chosen, 
*     correct the light strucure functions with the gluon raditation term.
*
         IF(VFNS_TMP.EQ."FFNS")THEN
            VFNS = "FFNS"
            DO I=4,6
               Q2H = Q2TH(I)
               IF(W2.GT.Q2H)THEN
                  ADSR = ADLERSR(Q2,Q2H)
                  INTEG2(1) = INTEG2(1)
     1                    + BQ(1) * ( DGAUSS(XDP_XC2Q_UNPLUS_NC,A,B,EPS)
     2                    - ADSR * XDP(X,Q) )
                  INTEG2(2) = INTEG2(2)
     1                    + BQ(2) * ( DGAUSS(XUP_XC2Q_UNPLUS_NC,A,B,EPS)
     2                    - ADSR * XUP(X,Q) )
                  INTEG2(3) = INTEG2(3)
     1                    + BQ(3) * ( DGAUSS(XSP_XC2Q_UNPLUS_NC,A,B,EPS)
     2                    - ADSR * XSP(X,Q) )
*
                  INTEGL(1) = INTEGL(1)
     1                      + BQ(1) * DGAUSS(XDP_XCLQ_UNPLUS_NC,A,B,EPS)
                  INTEGL(2) = INTEGL(2)
     1                      + BQ(2) * DGAUSS(XUP_XCLQ_UNPLUS_NC,A,B,EPS)
                  INTEGL(3) = INTEGL(3)
     1                      + BQ(3) * DGAUSS(XSP_XCLQ_UNPLUS_NC,A,B,EPS)
               ENDIF
            ENDDO
         ELSEIF(VFNS_TMP.EQ."FFN0")THEN
            VFNS = "FFN0"
            DO I=4,6
               Q2H = Q2TH(I)
               IF(Q2.GT.Q2H)THEN
                  ADSR = ADLERSR(Q2,Q2H)
                  INTEG2(1) = INTEG2(1)
     1                    + BQ(1) * ( DGAUSS(XDP_XC2Q_UNPLUS_NC,A,B,EPS)
     2                    + DGAUSS(XDP_XC2Q_PLUS_NC,A,B,EPS) 
     3                    + ADSR * XDP(X,Q) )
                  INTEG2(2) = INTEG2(2)
     1                    + BQ(2) * ( DGAUSS(XUP_XC2Q_UNPLUS_NC,A,B,EPS)
     2                    + DGAUSS(XUP_XC2Q_PLUS_NC,A,B,EPS) 
     3                    + ADSR * XUP(X,Q) )
                  INTEG2(3) = INTEG2(3)
     1                    + BQ(3) * ( DGAUSS(XSP_XC2Q_UNPLUS_NC,A,B,EPS)
     2                    + DGAUSS(XSP_XC2Q_PLUS_NC,A,B,EPS) 
     3                    + ADSR * XSP(X,Q) )
*
                  INTEGL(1) = INTEGL(1)
     1                      + BQ(1) * DGAUSS(XDP_XCLQ_UNPLUS_NC,A,B,EPS)
                  INTEGL(2) = INTEGL(2)
     1                      + BQ(2) * DGAUSS(XUP_XCLQ_UNPLUS_NC,A,B,EPS)
                  INTEGL(3) = INTEGL(3)
     1                      + BQ(3) * DGAUSS(XSP_XCLQ_UNPLUS_NC,A,B,EPS)
               ENDIF
            ENDDO
         ELSEIF(VFNS_TMP.EQ."FONLL")THEN
            DO I=4,6
               Q2H  = Q2TH(I)
               VFNS = "FFN0"
               IF(Q2.GT.Q2H)THEN
                  ADSR = ADLERSR(Q2,Q2H)
                  INTEG2(1) = DAMP(I) * ( INTEG2(1)
     1                    - BQ(1) * ( DGAUSS(XDP_XC2Q_UNPLUS_NC,A,B,EPS)
     2                    + DGAUSS(XDP_XC2Q_PLUS_NC,A,B,EPS) 
     3                    + ADSR * XDP(X,Q) ) )
                  INTEG2(2) = DAMP(I) * ( INTEG2(2)
     1                    - BQ(2) * ( DGAUSS(XUP_XC2Q_UNPLUS_NC,A,B,EPS)
     2                    + DGAUSS(XUP_XC2Q_PLUS_NC,A,B,EPS) 
     3                    + ADSR * XUP(X,Q) ) )
                  INTEG2(3) = DAMP(I) * ( INTEG2(3)
     1                    - BQ(3) * ( DGAUSS(XSP_XC2Q_UNPLUS_NC,A,B,EPS)
     2                    + DGAUSS(XSP_XC2Q_PLUS_NC,A,B,EPS) 
     3                    + ADSR * XSP(X,Q) ) )
*
                  INTEGL(1) = DAMP(I) * ( INTEGL(1)
     1                    - BQ(1) * DGAUSS(XDP_XCLQ_UNPLUS_NC,A,B,EPS) )
                  INTEGL(2) = DAMP(I) * ( INTEGL(2)
     1                    - BQ(2) * DGAUSS(XUP_XCLQ_UNPLUS_NC,A,B,EPS) )
                  INTEGL(3) = DAMP(I) * ( INTEGL(3)
     1                    - BQ(3) * DGAUSS(XSP_XCLQ_UNPLUS_NC,A,B,EPS) )
               ENDIF
*
               VFNS = "FFNS"
               IF(W2.GT.Q2H)THEN
                  ADSR = ADLERSR(Q2,Q2H)
                  INTEG2(1) = INTEG2(1)
     1                    + BQ(1) * ( DGAUSS(XDP_XC2Q_UNPLUS_NC,A,B,EPS)
     2                    - ADSR * XDP(X,Q) )
                  INTEG2(2) = INTEG2(2)
     1                    + BQ(2) * ( DGAUSS(XUP_XC2Q_UNPLUS_NC,A,B,EPS)
     2                    - ADSR * XUP(X,Q) )
                  INTEG2(3) = INTEG2(3)
     1                    + BQ(3) * ( DGAUSS(XSP_XC2Q_UNPLUS_NC,A,B,EPS)
     2                    - ADSR * XSP(X,Q) )
*
                  INTEGL(1) = INTEGL(1)
     1                      + BQ(1) * DGAUSS(XDP_XCLQ_UNPLUS_NC,A,B,EPS)
                  INTEGL(2) = INTEGL(2)
     1                      + BQ(2) * DGAUSS(XUP_XCLQ_UNPLUS_NC,A,B,EPS)
                  INTEGL(3) = INTEGL(3)
     1                      + BQ(3) * DGAUSS(XSP_XCLQ_UNPLUS_NC,A,B,EPS)
               ENDIF
            ENDDO
         ENDIF
*
         NNLOQ2L = 0D0
         NNLOQ3L = 0D0
         NNLOQLL = 0D0
         DO IPDF=1,3
            NNLOQ2L = NNLOQ2L + INTEG2(IPDF)
            NNLOQ3L = NNLOQ3L + INTEG3(IPDF)
            NNLOQLL = NNLOQLL + INTEGL(IPDF)
         ENDDO
*
         NNLO2L = NNLOG2L + NNLOQ2L
         NNLO3L = NNLOG3L + NNLOQ3L
         NNLOLL = NNLOGLL + NNLOQLL
      ENDIF
*
*     Gather all the pieces
*
      F2L = LO2L + ASQ * NLO2L + ASQ2 * NNLO2L
      F3L = LO3L + ASQ * NLO3L + ASQ2 * NNLO3L
      FLL = LOLL + ASQ * NLOLL + ASQ2 * NNLOLL
*
      VFNS = VFNS_TMP
*
*     In case GM scheme is chosen
*
      IGM   = 1
      WVFNS = VFNS
 205  IF(WVFNS.EQ."FONLL")THEN
         VFNS = IVFNS(IGM)
      ENDIF
*
************************************************************************
*
*     Heavy Parts
*
************************************************************************
*
*     LO
*
      IF(VFNS.EQ."ZMVN")THEN
         LO2C = BQ(4) * XCP(X,Q)
         LO2B = BQ(5) * XBP(X,Q)
         LO2T = BQ(6) * XTP(X,Q)
      ELSE
         LO2C = 0D0
         LO2C = 0D0
         LO2C = 0D0
      ENDIF
*
      IF(VFNS.EQ."ZMVN")THEN
         LO3C = DQ(4) * XCM(X,Q)
         LO3B = DQ(5) * XBM(X,Q)
         LO3T = DQ(6) * XTM(X,Q)
      ELSE
         LO3C = 0D0
         LO3C = 0D0
         LO3C = 0D0
      ENDIF
*
      LOLC = 0D0
      LOLB = 0D0
      LOLT = 0D0
*
*     NLO
*
*     Initialize
*
      NLO2C = 0D0
      NLO3C = 0D0
      NLOLC = 0D0
*   
      NLO2B = 0D0
      NLO3B = 0D0
      NLOLB = 0D0
*   
      NLO2T = 0D0
      NLO3T = 0D0
      NLOLT = 0D0
*
      IF(PTO.GE.1)THEN
*
         IPT = 1
*
*        Gluon
*
         Q2H = Q2TH(4)
         NLOG2C = 2D0 * BQ(4) * DGAUSS(XGL_XC2G_NC,A,B,EPS)
         NLOG3C = 2D0 * DQ(4) * DGAUSS(XGL_XC3G_NC,A,B,EPS)
         NLOGLC = 2D0 * BQ(4) * DGAUSS(XGL_XCLG_NC,A,B,EPS)
*
         Q2H = Q2TH(5)
         NLOG2B = 2D0 * BQ(5) * DGAUSS(XGL_XC2G_NC,A,B,EPS)
         NLOG3B = 2D0 * DQ(5) * DGAUSS(XGL_XC3G_NC,A,B,EPS)
         NLOGLB = 2D0 * BQ(5) * DGAUSS(XGL_XCLG_NC,A,B,EPS)
*      
         Q2H = Q2TH(6)
         NLOG2T = 2D0 * BQ(6) * DGAUSS(XGL_XC2G_NC,A,B,EPS)
         NLOG3T = 2D0 * DQ(6) * DGAUSS(XGL_XC3G_NC,A,B,EPS)
         NLOGLT = 2D0 * BQ(6) * DGAUSS(XGL_XCLG_NC,A,B,EPS)
*
*        Quark
*
*        F2
*
         IF(VFNS.EQ."ZMVN")THEN
            NLOQ2C = BQ(4) * ( DGAUSS(XCP_XC2Q_UNPLUS_NC,A,B,EPS)
     1             + DGAUSS(XCP_XC2Q_PLUS_NC,A,B,EPS) 
     2             + ( CONST - PLUS_TERM2 ) * XCP(X,Q) )
            NLOQ2B = BQ(5) * ( DGAUSS(XBP_XC2Q_UNPLUS_NC,A,B,EPS)
     1             + DGAUSS(XBP_XC2Q_PLUS_NC,A,B,EPS) 
     2             + ( CONST - PLUS_TERM2 ) * XBP(X,Q) )
            NLOQ2T = BQ(6) * ( DGAUSS(XTP_XC2Q_UNPLUS_NC,A,B,EPS)
     1             + DGAUSS(XTP_XC2Q_PLUS_NC,A,B,EPS) 
     2             + ( CONST - PLUS_TERM2 ) * XTP(X,Q) )
         ELSE
            NLOQ2C = 0D0
            NLOQ2B = 0D0
            NLOQ2T = 0D0
         ENDIF
*
*        F3
*
         IF(VFNS.EQ."ZMVN")THEN
            NLOQ3C = DQ(4) * ( DGAUSS(XCM_XC3Q_UNPLUS_NC,A,B,EPS)
     1             + DGAUSS(XCM_XC3Q_PLUS_NC,A,B,EPS) 
     2             + ( CONST - PLUS_TERM3 ) * XCM(X,Q) )
            NLOQ3B = DQ(5) * ( DGAUSS(XBM_XC3Q_UNPLUS_NC,A,B,EPS)
     1             + DGAUSS(XBM_XC3Q_PLUS_NC,A,B,EPS) 
     2             + ( CONST - PLUS_TERM3 ) * XBM(X,Q) )
            NLOQ3T = DQ(6) * ( DGAUSS(XTM_XC3Q_UNPLUS_NC,A,B,EPS)
     1             + DGAUSS(XTM_XC3Q_PLUS_NC,A,B,EPS) 
     2             + ( CONST - PLUS_TERM3 ) * XTM(X,Q) )
         ELSE
            NLOQ3C = 0D0
            NLOQ3B = 0D0
            NLOQ3T = 0D0
         ENDIF
*
*        FL
*
         IF(VFNS.EQ."ZMVN")THEN
            NLOQLC = BQ(4) * DGAUSS(XCP_XCLQ_UNPLUS_NC,A,B,EPS)
            NLOQLB = BQ(5) * DGAUSS(XBP_XCLQ_UNPLUS_NC,A,B,EPS)
            NLOQLT = BQ(6) * DGAUSS(XTP_XCLQ_UNPLUS_NC,A,B,EPS)
         ELSE
            NLOQLC = 0D0
            NLOQLB = 0D0
            NLOQLT = 0D0
         ENDIF
*
         NLO2C = NLOG2C + NLOQ2C
         NLO3C = NLOG3C + NLOQ3C
         NLOLC = NLOGLC + NLOQLC
*
         NLO2B = NLOG2B + NLOQ2B
         NLO3B = NLOG3B + NLOQ3B
         NLOLB = NLOGLB + NLOQLB
*
         NLO2T = NLOG2T + NLOQ2T
         NLO3T = NLOG3T + NLOQ3T
         NLOLT = NLOGLT + NLOQLT
      ENDIF
*
*     NNLO
*
*     Initialize
*
      NNLO2C = 0D0
      NNLO3C = 0D0
      NNLOLC = 0D0
*   
      NNLO2B = 0D0
      NNLO3B = 0D0
      NNLOLB = 0D0
*   
      NNLO2T = 0D0
      NNLO3T = 0D0
      NNLOLT = 0D0
*
      IF(PTO.GE.2)THEN
*
         IPT = 2
*
*        Gluon + Pure Singlet
*
         Q2H = Q2TH(4)
         NNLOG2C = BQ(4) 
     1           * ( DGAUSS(XGL_XC2G_NC,A,B,EPS)
     2           +   DGAUSS(XPS_XC2Q_UNPLUS_NC,A,B,EPS) )
         NNLOG3C = DQ(4) * DGAUSS(XGL_XC3G_NC,A,B,EPS)
         NNLOGLC = BQ(4)
     1           * ( DGAUSS(XGL_XCLG_NC,A,B,EPS)
     2           +   DGAUSS(XPS_XCLQ_UNPLUS_NC,A,B,EPS) )
*
         Q2H = Q2TH(5)
         NNLOG2B = BQ(5) 
     1           * ( DGAUSS(XGL_XC2G_NC,A,B,EPS)
     2           +   DGAUSS(XPS_XC2Q_UNPLUS_NC,A,B,EPS) )
         NNLOG3B = DQ(5) * DGAUSS(XGL_XC3G_NC,A,B,EPS)
         NNLOGLB = BQ(5)
     1           * ( DGAUSS(XGL_XCLG_NC,A,B,EPS)
     2           +   DGAUSS(XPS_XCLQ_UNPLUS_NC,A,B,EPS) )
*      
         Q2H = Q2TH(6)
         NNLOG2T = BQ(6) 
     1           * ( DGAUSS(XGL_XC2G_NC,A,B,EPS)
     2           +   DGAUSS(XPS_XC2Q_UNPLUS_NC,A,B,EPS) )
         NNLOG3T = DQ(6) * DGAUSS(XGL_XC3G_NC,A,B,EPS)
         NNLOGLT = BQ(6)
     1           * ( DGAUSS(XGL_XCLG_NC,A,B,EPS)
     2           +   DGAUSS(XPS_XCLQ_UNPLUS_NC,A,B,EPS) )
*
*        In case of ZM-VFNS, add to F2 the articifial local term
*        due to the parametrization
*
         IF(VFNS.EQ."ZMVN")THEN
            NNLOG2C = NNLOG2C + BQ(4) * C2G2C(X,1) * XGL(X,Q)
            NNLOG2B = NNLOG2C + BQ(5) * C2G2C(X,1) * XGL(X,Q)
            NNLOG2T = NNLOG2C + BQ(6) * C2G2C(X,1) * XGL(X,Q)
         ENDIF
*
*        Quark (non-singlet)
*
*        F2
*
         IF(VFNS.EQ."ZMVN")THEN
            NNLOQ2C = BQ(4) * ( DGAUSS(XCP_XC2Q_UNPLUS_NC,A,B,EPS)
     1              + DGAUSS(XCP_XC2Q_PLUS_NC,A,B,EPS) 
     2              + LOCAL_TERM2 * XCP(X,Q) )
            NNLOQ2B = BQ(5) * ( DGAUSS(XBP_XC2Q_UNPLUS_NC,A,B,EPS)
     1              + DGAUSS(XBP_XC2Q_PLUS_NC,A,B,EPS) 
     2              + LOCAL_TERM2 * XBP(X,Q) )
            NNLOQ2T = BQ(6) * ( DGAUSS(XTP_XC2Q_UNPLUS_NC,A,B,EPS)
     1              + DGAUSS(XTP_XC2Q_PLUS_NC,A,B,EPS) 
     2              + LOCAL_TERM2 * XTP(X,Q) )
         ELSE
            NNLOQ2C = 0D0
            NNLOQ2B = 0D0
            NNLOQ2T = 0D0
         ENDIF
*
*        F3
*
         IF(VFNS.EQ."ZMVN")THEN
            NNLOQ3C = DQ(4) * ( DGAUSS(XCM_XC3Q_UNPLUS_NC,A,B,EPS)
     1              + DGAUSS(XCM_XC3Q_PLUS_NC,A,B,EPS) 
     2              + LOCAL_TERM3 * XCM(X,Q) )
            NNLOQ3B = DQ(5) * ( DGAUSS(XBM_XC3Q_UNPLUS_NC,A,B,EPS)
     1              + DGAUSS(XBM_XC3Q_PLUS_NC,A,B,EPS) 
     2              + LOCAL_TERM3 * XBM(X,Q) )
            NNLOQ3T = DQ(6) * ( DGAUSS(XTM_XC3Q_UNPLUS_NC,A,B,EPS)
     1              + DGAUSS(XTM_XC3Q_PLUS_NC,A,B,EPS) 
     2              + LOCAL_TERM3 * XTM(X,Q) )
         ELSE
            NNLOQ3C = 0D0
            NNLOQ3B = 0D0
            NNLOQ3T = 0D0
         ENDIF
*
*        FL
*
         IF(VFNS.EQ."ZMVN")THEN
            NNLOQLC = BQ(4) * ( DGAUSS(XCP_XCLQ_UNPLUS_NC,A,B,EPS)
     1              + LOCAL_TERML * XCP(X,Q) )
            NNLOQLB = BQ(5) * ( DGAUSS(XBP_XCLQ_UNPLUS_NC,A,B,EPS)
     1              + LOCAL_TERML * XBP(X,Q) )
            NNLOQLT = BQ(6) * ( DGAUSS(XTP_XCLQ_UNPLUS_NC,A,B,EPS)
     1              + LOCAL_TERML * XTP(X,Q) )
         ELSE
            NNLOQLC = 0D0
            NNLOQLB = 0D0
            NNLOQLT = 0D0
         ENDIF
*
         NNLO2C = NNLOG2C + NNLOQ2C
         NNLO3C = NNLOG3C + NNLOQ3C
         NNLOLC = NNLOGLC + NNLOQLC
*
         NNLO2B = NNLOG2B + NNLOQ2B
         NNLO3B = NNLOG3B + NNLOQ3B
         NNLOLB = NNLOGLB + NNLOQLB
*
         NNLO2T = NNLOG2T + NNLOQ2T
         NNLO3T = NNLOG3T + NNLOQ3T
         NNLOLT = NNLOGLT + NNLOQLT
      ENDIF
*
      F2C = LO2C + ASQ * NLO2C + ASQ2 * NNLO2C
      F3C = LO3C + ASQ * NLO3C + ASQ2 * NNLO3C
      FLC = LOLC + ASQ * NLOLC + ASQ2 * NNLOLC
*
      F2B = LO2B + ASQ * NLO2B + ASQ2 * NNLO2B 
      F3B = LO3B + ASQ * NLO3B + ASQ2 * NNLO3B 
      FLB = LOLB + ASQ * NLOLB + ASQ2 * NNLOLB 
*
      F2T = LO2T + ASQ * NLO2T + ASQ2 * NNLO2T 
      F3T = LO3T + ASQ * NLO3T + ASQ2 * NNLO3T 
      FLT = LOLT + ASQ * NLOLT + ASQ2 * NNLOLT 
*
************************************************************************
*
*     Join the results in case of GM-VFNS
*
      IF(WVFNS.EQ."FONLL")THEN
         F2CG(IGM) = F2C
         F3CG(IGM) = F3C
         FLCG(IGM) = FLC
*
         F2BG(IGM) = F2B
         F3BG(IGM) = F3B
         FLBG(IGM) = FLB
*
         F2TG(IGM) = F2T
         F3TG(IGM) = F3T
         FLTG(IGM) = FLT
*
         IF(IGM.LE.2)THEN
            IGM = IGM + 1
            GOTO 205
         ENDIF
*
         F2C = F2CG(1) + DAMP(4) * ( F2CG(2) - F2CG(3) )
         F2B = F2BG(1) + DAMP(5) * ( F2BG(2) - F2BG(3) )
         F2T = F2TG(1) + DAMP(6) * ( F2TG(2) - F2TG(3) )
*
         FLC = FLCG(1) + DAMP(4) * ( FLCG(2) - FLCG(3) )
         FLB = FLBG(1) + DAMP(5) * ( FLBG(2) - FLBG(3) )
         FLT = FLTG(1) + DAMP(6) * ( FLTG(2) - FLTG(3) )
*
         F3C = F3CG(1) + DAMP(4) * ( F3CG(2) - F3CG(3) )
         F3B = F3BG(1) + DAMP(5) * ( F3BG(2) - F3BG(3) )
         F3T = F3TG(1) + DAMP(6) * ( F3TG(2) - F3TG(3) )
      ENDIF
*
      IF(IE.EQ.1)THEN       !Positron
         SIGMAL = F2L - ( Y**2D0 / YP ) * FLL - ( YM / YP ) * F3L
         SIGMAC = F2C - ( Y**2D0 / YP ) * FLC - ( YM / YP ) * F3C
         SIGMAB = F2B - ( Y**2D0 / YP ) * FLB - ( YM / YP ) * F3B
         SIGMAT = F2T - ( Y**2D0 / YP ) * FLT - ( YM / YP ) * F3T
      ELSEIF(IE.EQ.-1)THEN   !Electron
         SIGMAL = F2L - ( Y**2D0 / YP ) * FLL + ( YM / YP ) * F3L
         SIGMAC = F2C - ( Y**2D0 / YP ) * FLC + ( YM / YP ) * F3C
         SIGMAB = F2B - ( Y**2D0 / YP ) * FLB + ( YM / YP ) * F3B
         SIGMAT = F2T - ( Y**2D0 / YP ) * FLT + ( YM / YP ) * F3T
      ENDIF
*
      IF(Q2.GT.Q2TH(6))THEN
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
      SIGMA(3) = SIGMAL
      SIGMA(4) = SIGMAC
      SIGMA(5) = SIGMAB
      SIGMA(6) = SIGMAT
      SIGMA(7) = SIGMAP
*
      RETURN
      END
