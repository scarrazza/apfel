************************************************************************
*
*     NEUTRAL CURRENT
*
************************************************************************
*
*     MassiveZeroCoefficientFunctions.f:
*
*     This file contains all the exact asymptotic massive coefficient
*     functions as given in the paper below.
*
*     Reference: hep-ph/9601302 (Appendix D)
*
*     The following return the coefficients of the log dependent terms.
*
************************************************************************
*
*     Order alphas coeficient functions (NLO)
*     Expansion parameter alphas/4*pi
*
************************************************************************
*     Equation (4.1): Gluon coefficient function for FL
************************************************************************
      FUNCTION CLG1AM0_A0(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Output Variables
*
      DOUBLE PRECISION CLG1AM0_A0
*
      CLG1AM0_A0 = TR * 16D0 * X * ( 1D0 - X )
*
      RETURN
      END
*
************************************************************************
*     Equation (4.2): Gluon coefficient function for F2
************************************************************************
      FUNCTION C2G1AM0_AQ(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Output Variables
*
      DOUBLE PRECISION C2G1AM0_AQ
*
      C2G1AM0_AQ = TR * ( 4D0 - 8D0 * X + 8D0 * X**2 )
*
      RETURN
      END
*
************************************************************************
      FUNCTION C2G1AM0_A0(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Output Variables
*
      DOUBLE PRECISION C2G1AM0_A0
*
      C2G1AM0_A0 = TR * ( ( 4D0 - 8D0 * X + 8D0 * X**2 )
     1           * DLOG( ( 1D0 - X ) / X ) - 4D0 + 32D0 * X 
     2           - 32D0 * X**2 )
*
      RETURN
      END
*
************************************************************************
*
*     Order alphas^2 coeficient functions (NNLO)
*     Expansion parameter alphas/4*pi 
*
************************************************************************
*     Equation (4.3): Gluon coefficient function for FL
************************************************************************
      FUNCTION CLG2AM0_AF(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION DLX,DLM
      DOUBLE PRECISION C1
**
*     Output Variables
*
      DOUBLE PRECISION CLG2AM0_AF
*
      DLX=DLOG(X)
      DLM=DLOG(1.0D0-X)
*
      C1=64.0D0*X*(1.0D0-X)*DLM-128.0D0*X*DLX-32.0D0-160.0D0*X
     1+544.0D0*X*X/3.0D0+32.0D0/X/3.0D0
*
      CLG2AM0_AF = TR * CA * C1
*
      RETURN
      END
*
************************************************************************
      FUNCTION CLG2AM0_AQ(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION DLX,DLM
      DOUBLE PRECISION A1,B1
**
*     Output Variables
*
      DOUBLE PRECISION CLG2AM0_AQ
*
      DLX=DLOG(X)
      DLM=DLOG(1.0D0-X)
*
      A1=64.0D0*X*(1.0D0-X)*DLM-128.0D0*X*DLX-32.0D0-160.0D0*X
     1+544.0D0*X*X/3.0D0+32.0D0/X/3.0D0
*
      B1=32.0D0*X*DLX+16.0D0*(1.0D0-2.0D0*X*X+X)
*
      CLG2AM0_AQ = TR * ( CA * A1 + CF * B1 )
*
      RETURN
      END
*
************************************************************************
      FUNCTION CLG2AM0_A0(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
      INCLUDE "../commons/consts.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION DLX,DLX2,DLM,DLM2,DLP,S11,S11M
      DOUBLE PRECISION A2,B2
      DOUBLE PRECISION WGPLG
**
*     Output Variables
*
      DOUBLE PRECISION CLG2AM0_A0
*
      DLX=DLOG(X)
      DLX2=DLX*DLX
      DLM=DLOG(1.0D0-X)
      DLM2=DLM*DLM
      DLP=DLOG(1.0D0+X)
      S11=WGPLG(1,1,1.0D0-X)
      S11M=WGPLG(1,1,-X)
*
      A2=96.0D0*X*DLX2+(64.0D0*X*X-192.0D0*X)*DLX*DLM+(32.0D0
     1-416.0D0*X*X+256.0D0*X)*DLX+32.0D0*X*(1.0D0-X)*DLM2+(
     2-32.0D0+928.0D0*X*X/3.0D0-288.0D0*X+32.0D0/X/3.0D0)*DLM
     3+64.0D0*X*X*ZETA2-128.0D0*X*S11+64.0D0*X*(1.0D0+X)*(S11M
     4+DLX*DLP)+32.0D0/3.0D0-1696*X*X/9.0D0+544.0D0*X/3.0D0
     5-32.0D0/X/9.0D0
*
      B2=-(64.0D0*X*X*X/5.0D0+64.0D0*X/3.0D0)*DLX2+32.0D0*X
     1*(S11+DLX*DLM)+(-208.0D0/15.0D0+192.0D0*X*X/5.0D0
     2-416.0D0*X/5.0D0-64.0D0/X/15.0D0)*DLX+(16.0D0-64.0D0*X*X
     3+48.0D0*X)*DLM+(128.0D0*X*X*X/5.0D0-64.0D0*X/3.0D0)*ZETA2
     4+(128.0D0*X*X*X/5.0D0-64.0D0*X/3.0D0+64.0D0/X/X/15.0D0)
     5*(S11M+DLX*DLP)-256.0D0/15.0D0+672.0D0*X*X/5.0D0-608.0D0
     6*X/5.0D0+64.0D0/X/15.0D0
*
      CLG2AM0_A0 = TR * ( CA * A2 + CF * B2 )
*
      RETURN
      END
*
************************************************************************
*     Equation (4.4): Gluon coefficient function for F2
************************************************************************
      FUNCTION C2G2AM0_AQF(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION DLX,DLM
      DOUBLE PRECISION C1
**
*     Output Variables
*
      DOUBLE PRECISION C2G2AM0_AQF
*
      DLX=DLOG(X)
      DLM=DLOG(1.0D0-X)
*
      C1=-248.0D0*X*X/3.0D0+64.0D0*X+8.0D0+32.0D0/X/3.0D0
     1+(32.0D0*X*X-32.0D0*X+16.0D0)*DLM+(64.0D0*X+16.0D0)*DLX
*
      C2G2AM0_AQF = TR * CA * C1
*
      RETURN
      END
*
************************************************************************
      FUNCTION C2G2AM0_AF(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
      INCLUDE "../commons/consts.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION DLX,DLX2,DLM,DLM2
      DOUBLE PRECISION S11
      DOUBLE PRECISION C2
      DOUBLE PRECISION WGPLG
**
*     Output Variables
*
      DOUBLE PRECISION C2G2AM0_AF
*
      DLX=DLOG(X)
      DLX2=DLX*DLX
      DLM=DLOG(1.0D0-X)
      DLM2=DLM*DLM
      S11=WGPLG(1,1,1.0D0-X)
*
      C2=(-16.0D0+32.0D0*X-32.0D0*X*X)*ZETA2+1124.0D0*X*X/3.0D0
     1-968.0D0*X/3.0D0-172.0D0/3.0D0+16.0D0/X/3.0D0+(16.0D0
     2-32.0D0*X+32.0D0*X*X)*DLM2-(8.0D0+32.0D0*X)*DLX2+(248.0D0
     3*X*X/3.0D0-256.0D0*X-8.0D0)*DLX+(32.0D0/X/3.0D0-8.0D0
     4+192.0D0*X-632.0D0*X*X/3.0D0)*DLM+(96.0D0*X-32.0D0*X*X)
     5*DLX*DLM+(64.0D0*X+16.0D0)*S11
*
      C2G2AM0_AF  = TR * CA * C2
*
      RETURN
      END
*
************************************************************************
      FUNCTION C2G2AM0_AQ2(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION DLX,DLM
      DOUBLE PRECISION A1,B1
**
*     Output Variables
*
      DOUBLE PRECISION C2G2AM0_AQ2
*
      DLX=DLOG(X)
      DLM=DLOG(1.0D0-X)
*
      A1=16.0D0/X/3.0D0-124.0D0*X*X/3.0D0+32.0D0*X+4.0D0
     1+(16.0D0*X*X-16.0D0*X+8.0D0)*DLM+(32.0D0*X+8.0D0)*DLX
*
      B1=-2.0D0+8.0D0*X+(16.0D0*X*X-16.0D0*X+8.0D0)*DLM
     1+(-16.0D0*X*X+8.0D0*X-4.0D0)*DLX
*
      C2G2AM0_AQ2 = TR * ( CA * A1 + CF * B1 )
*
      RETURN
      END
*
************************************************************************
      FUNCTION C2G2AM0_AQ(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
      INCLUDE "../commons/consts.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION DLX,DLX2,DLM,DLM2,DLP
      DOUBLE PRECISION S111MX,S11MX
      DOUBLE PRECISION A2,B2
      DOUBLE PRECISION WGPLG
**
*     Output Variables
*
      DOUBLE PRECISION C2G2AM0_AQ
*
      DLX=DLOG(X)
      DLX2=DLX*DLX
      DLM=DLOG(1.0D0-X)
      DLM2=DLM*DLM
      DLP=DLOG(1.0D0+X)
      S111MX=WGPLG(1,1,1.0D0-X)
      S11MX=WGPLG(1,1,-X)
*
      A2=-(16.0D0+32.0D0*X*X)*ZETA2+1628.0D0*X*X/9.0D0
     1-368.0D0*X/3.0D0-220.0D0/3.0D0+208.0D0/X/9.0D0
     2+(16.0D0*X*X-16.0D0*X+8.0D0)*DLM2-(48.0D0*X+16.0D0)
     3*DLX2+(-536.0D0*X*X/3.0D0+160.0D0*X-8.0D0+32.0D0/X/3.0D0)
     4*DLM+(200.0D0*X*X-192.0D0*X)*DLX+(96.0D0*X-32.0D0*X*X)
     5*DLX*DLM+(64.0D0*X+16.0D0)*S111MX-(32.0D0*X*X+32.0D0*X+16.0D0)
     6*(S11MX+DLX*DLP)
*
      B2=(-64.0D0*X*X+64.0D0*X-32.0D0)*ZETA2+16.0D0*X*X-68.0D0*X
     1+36.0D0+(32.0D0*X*X-32.0D0*X+16.0D0)*DLM2+(32.0D0*X*X
     2-16.0D0*X+8.0D0)*DLX2+(-80.0D0*X*X+96.0D0*X-28.0D0)*DLM
     3+(80.0D0*X*X-48.0D0*X+8.0D0)*DLX+(-64.0D0*X*X+48.0D0*X
     4-24.0D0)*DLX*DLM+(8.0D0-16.0D0*X)*S111MX
*
      C2G2AM0_AQ  = TR * ( CA * A2 + CF * B2 )
*
      RETURN
      END
*
************************************************************************
      FUNCTION C2G2AM0_A0(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
      INCLUDE "../commons/consts.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION DLX,DLX2,DLX3,DLM,DLM2,DLM3,DLP,DLP2
      DOUBLE PRECISION S11,S121MX,S12MX,S211MX,S21MX,S111MX,S11MX
      DOUBLE PRECISION Z,S21Z,S21MZ
      DOUBLE PRECISION A31,A32,A33,A34,A3
      DOUBLE PRECISION B31,B32,B33,B34,B3
      DOUBLE PRECISION WGPLG
**
*     Output Variables
*
      DOUBLE PRECISION C2G2AM0_A0
*
      DLX=DLOG(X)
      DLX2=DLX*DLX
      DLX3=DLX2*DLX
      DLM=DLOG(1.0D0-X)
      DLM2=DLM*DLM
      DLM3=DLM2*DLM
      DLP=DLOG(1.0D0+X)
      DLP2=DLP*DLP
      S11=WGPLG(1,1,1.0D0-X)
      S121MX=WGPLG(1,2,1.0D0-X)
      S12MX=WGPLG(1,2,-X)
      S211MX=WGPLG(2,1,1.0D0-X)
      S21MX=WGPLG(2,1,-X)
      S111MX=WGPLG(1,1,1.0D0-X)
      S11MX=WGPLG(1,1,-X)
      Z=(1.0D0-X)/(1.0D0+X)
      S21Z=WGPLG(2,1,Z)
      S21MZ=WGPLG(2,1,-Z)
*
      A31=DLX3*(16.0D0/3.0D0+16.0D0*X)+
     1DLX2*DLM*(-8.0D0+16.0D0*X*X-64.0D0*X)+
     2DLX2*DLP*(12.0D0+32.0D0*X*X+24.0D0*X)+
     3DLX2*(-114.0D0*X*X+184.0D0*X)+
     4DLX*DLM2*(-16.0D0*X*X+48.0D0*X)+
     5DLX*DLM*DLP*(-16.0D0-32.0D0*X*X-32.0D0*X)+
     6DLX*DLM*(16.0D0+292.0D0*X*X-288.0D0*X)+
     7DLX*DLP2*(8.0D0+16.0D0*X)
      A32=DLX*DLP*(-48.0D0+208.0D0*X*X/3.0D0+16.0D0*X-
     132.0D0/X/3.0D0)+
     2DLX*ZETA2*(-16.0D0+32.0D0*X*X-160.0D0*X)+
     3DLX*S111MX*(+32.0D0*X)+
     4DLX*S11MX*(24.0D0+32.0D0*X*X+48.0D0*X)+
     5DLX*(292.0D0/3.0D0-5780.0D0*X*X/9.0D0+332.0D0*X)+
     6DLM2*(-6.0D0-214.0D0*X*X/3.0D0+64.0D0*X+16.0D0/X/3.0D0)+
     7ZETA2*DLM*(-40.0D0-64.0D0*X*X+48.0D0*X)
      A33=DLM*S111MX*(16.0D0+64.0D0*X)+
     1DLM*S11MX*(-16.0D0-32.0D0*X*X-32.0D0*X)+
     2DLM*(-112.0D0/3.0D0+2996.0D0*X*X/9.0D0-860.0D0*X/3.0D0
     3+208.0D0/X/9.0D0)+
     4ZETA2*DLP*(8.0D0+16.0D0*X)+
     5DLP*S11MX*(16.0D0+32.0D0*X)+
     6S21MZ*(-16.0D0-32.0D0*X*X-32.0D0*X)+
     7S21Z*(16.0D0+32.0D0*X*X+32.0D0*X)
      A34=ZETA2*(-4.0D0+796.0D0*X*X/3.0D0-208.0D0*X-32.0D0/X)+
     1ZETA3*(-12.0D0-8.0D0*X*X-56.0D0*X)+
     2S111MX*(20.0D0+80.0D0*X*X/3.0D0-64.0D0*X+64.0D0/X/3.0D0)+
     3S211MX*(-16.0D0-128.0D0*X)+
     4S121MX*(40.0D0+144.0D0*X)+
     5S11MX*(-48.0D0+208.0D0*X*X/3.0D0+16.0D0*X-32.0D0/X/3.0D0)+
     6S21MX*(-24.0D0-48.0D0*X)+
     7S12MX*(16.0D0+32.0D0*X)+80.0D0/X/9.0D0+
     8466.0D0/9.0D0-878.0D0*X*X/9.0D0+260.0D0*X/9.0D0
      A3=A31+A32+A33+A34
*
      B31=DLX3*(-8.0D0/3.0D0-32.0D0*X*X/3.0D0+16.0D0*X/3.0D0)+
     1DLX2*DLM*(16.0D0+48.0D0*X*X-32.0D0*X)+
     2DLX2*DLP*(16.0D0+16.0D0*X*X+32.0D0*X)+
     3DLX2*(-4.0D0-96.0D0*X*X*X/5.0D0-52.0D0*X*X+8.0D0*X/3.0D0)+
     4DLX*DLM2*(-20.0D0-48.0D0*X*X+40.0D0*X)+
     5DLX*DLM*(24.0D0+168.0D0*X*X-160.0D0*X)+
     6DLX*DLP2*(-32.0D0-32.0D0*X*X-64.0D0*X)+
     7DLX*DLP*(96.0D0+192.0D0*X*X*X/5.0D0+128.0D0*X/3.0D0)
      B32=16.0D0*DLX*DLP/X/X/15.0D0-16.0D0*DLX/X/15.0D0+
     1DLX*ZETA2*(32.0D0+64.0D0*X*X-64.0D0*X)+
     2DLX*S111MX*(32.0D0*X*X)+
     3DLX*S11MX*(-32.0D0-32.0D0*X*X+64.0D0*X)+
     4DLX*(-712.0D0/15.0D0-672.0D0*X*X/5.0D0+136.0D0*X/5.0D0)+
     5DLM3*(8.0D0+16.0D0*X*X-16.0D0*X)+
     6DLM2*(-22.0D0-84.0D0*X*X+88.0D0*X)+
     7ZETA2*DLM*(-32.0D0*X*X)
      B33=DLM*S111MX*(8.0D0-16.0D0*X)+
     1DLM*(28.0D0+96.0D0*X*X-132.0D0*X)+
     2ZETA2*DLP*(-32.0D0-32.0D0*X*X-64.0D0*X)+
     3DLP*S11MX*(-64.0D0-64.0D0*X*X-128.0D0*X)+
     4ZETA2*(48.0D0+192.0D0*X*X*X/5.0D0+104.0D0*X*X-
     5208.0D0*X/3.0D0)+
     6ZETA3*(112.0D0+192.0D0*X*X-96.0D0*X)
      B34=S111MX*(-24.0D0+64.0D0*X*X-48.0D0*X)+
     1S211MX*(-24.0D0-32.0D0*X*X+48.0D0*X)+
     2S121MX*(-32.0D0+64.0D0*X)+
     3S11MX*(96.0D0+192.0D0*X*X*X/5.0D0+128.0D0*X/3.0D0)+
     416.0D0*S11MX/X/X/15.0D0+16.0D0/X/15.0D0+
     5S21MX*(96.0D0+96.0D0*X*X-64.0D0*X)+
     6S12MX*(-64.0D0-64.0D0*X*X-128.0D0*X)-
     7904.0D0/15.0D0+328.0D0*X*X/5.0D0+68.0D0*X/5.0D0
      B3=B31+B32+B33+B34
*
      C2G2AM0_A0  = TR * ( CA * A3 + CF * B3 )
*
      RETURN
      END
*
************************************************************************
*     Equation (4.5): Heavy quark coefficient function for FL
************************************************************************
      FUNCTION CLPS2AM0_AF(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION DLX
      DOUBLE PRECISION C1
**
*     Output Variables
*
      DOUBLE PRECISION CLPS2AM0_AF
*
      DLX=DLOG(X)
*
      C1=32.0D0/3.0D0/X+64.0D0*X*X/3.0D0
     1-32.0D0-32.0D0*X*DLX
*
      CLPS2AM0_AF = CF * TR * C1
*
      RETURN
      END
*
************************************************************************
      FUNCTION CLPS2AM0_AQ(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION DLX
      DOUBLE PRECISION A1
**
*     Output Variables
*
      DOUBLE PRECISION CLPS2AM0_AQ
*
      DLX=DLOG(X)
*
      A1=-32.0D0*X*DLX-32.0D0+64.0D0*X*X/3.0D0+32.0D0/X/3.0D0
*
      CLPS2AM0_AQ = CF * TR * A1
*
      RETURN
      END
*
************************************************************************
      FUNCTION CLPS2AM0_A0(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION DLX,DLX2,DLM,SPX
      DOUBLE PRECISION A2
      DOUBLE PRECISION WGPLG
**
*     Output Variables
*
      DOUBLE PRECISION CLPS2AM0_A0
*
      DLX=DLOG(X)
      DLX2=DLX*DLX
      DLM=DLOG(1.0D0-X)
      SPX=WGPLG(1,1,1.0D0-X)
*
      A2=32.0D0*X*(DLX2-DLX*DLM-SPX)+(-32.0D0+64.0D0*X*X/3.0D0
     1+32.0D0/X/3.0D0)*DLM+(32.0D0-64.0D0*X*X-32.0D0*X)*DLX
     2+32.0D0/3.0D0+320.0D0*X*X/9.0D0-128.0D0*X/3.0D0
     3-32.0D0/X/9.0D0
*
      CLPS2AM0_A0 = CF * TR * A2
*
      RETURN
      END
*
************************************************************************
*     Equation (4.6): Heavy quark coefficient function for F2
************************************************************************
      FUNCTION C2PS2AM0_AQF(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION DLX
      DOUBLE PRECISION C1
**
*     Output Variables
*
      DOUBLE PRECISION C2PS2AM0_AQF
*
      DLX=DLOG(X)
*
      C1=-32.0D0*X*X/3.0D0-8.0D0*X+8.0D0+32.0D0/X/3.0D0
     1+16.0D0*(1.0D0+X)*DLX
*
      C2PS2AM0_AQF = CF * TR * C1
*
      RETURN
      END
*
************************************************************************
      FUNCTION C2PS2AM0_AF(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION DLX,DLX2,DLM,SPX
      DOUBLE PRECISION C2
      DOUBLE PRECISION WGPLG
**
*     Output Variables
*
      DOUBLE PRECISION C2PS2AM0_AF
*
      DLX=DLOG(X)
      DLX2=DLX*DLX
      DLM=DLOG(1.0D0-X)
      SPX=WGPLG(1,1,1.0D0-X)
*
      C2=128.0D0*X*X/3.0D0+16.0D0*X/3.0D0-160.0D0/3.0D0
     1+16.0D0/X/3.0D0+8.0D0*(1.0D0+X)*(-DLX2+2.0D0*DLX*DLM
     2+2.0D0*SPX)+(-32.0D0*X*X/3.0D0-8.0D0*X+8.0D0+32.0D0/X/3.0D0)
     3*DLM+(32.0D0*X*X/3.0D0-40.0D0*X-8.0D0)*DLX
*
      C2PS2AM0_AF  = CF * TR * C2
*
      RETURN
      END
*
************************************************************************
      FUNCTION C2PS2AM0_AQ2(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION DLX
      DOUBLE PRECISION A1
**
*     Output Variables
*
      DOUBLE PRECISION C2PS2AM0_AQ2
*
      DLX=DLOG(X)
*
      A1=-16.0D0*X*X/3.0D0-4.0D0*X+4.0D0+16.0D0/X/3.0D0
     1+8.0D0*(1.0D0+X)*DLX
*
      C2PS2AM0_AQ2 = CF * TR * A1
*
      RETURN
      END
*
************************************************************************
      FUNCTION C2PS2AM0_AQ(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION DLX,DLX2,DLM,S11
      DOUBLE PRECISION A2
      DOUBLE PRECISION WGPLG
**
*     Output Variables
*
      DOUBLE PRECISION C2PS2AM0_AQ
*
      DLX=DLOG(X)
      DLX2=DLX*DLX
      DLM=DLOG(1.0D0-X)
      S11=WGPLG(1,1,1.0D0-X)
*
      A2=-64.0D0*X*X/9.0D0+160.0D0*X/3.0D0-208.0D0/3.0D0
     1+208.0D0/X/9.0D0+32.0D0*X*X*DLX+16.0D0*(1.0D0+X)*(-DLX2+DLX*DLM
     2+S11)+(-32.0D0*X*X/3.0D0-8.0D0*X+8.0D0+32.0D0/X/3.0D0)
     3*DLM
*
      C2PS2AM0_AQ  = CF * TR * A2
*
      RETURN
      END
*
************************************************************************
      FUNCTION C2PS2AM0_A0(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
      INCLUDE "../commons/consts.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION DLX,DLX2,DLX3,DLM,DLM2,DLP,S11,S12,S21,S11M
      DOUBLE PRECISION A3
      DOUBLE PRECISION WGPLG
**
*     Output Variables
*
      DOUBLE PRECISION C2PS2AM0_A0
*
      DLX=DLOG(X)
      DLX2=DLX*DLX
      DLX3=DLX*DLX2
      DLM=DLOG(1.0D0-X)
      DLM2=DLM*DLM
      DLP=DLOG(1.0D0+X)
      S11=WGPLG(1,1,1.0D0-X)
      S12=WGPLG(1,2,1.0D0-X)
      S21=WGPLG(2,1,1.0D0-X)
      S11M=WGPLG(1,1,-X)
*
      A3=(1.0D0+X)*(16.0D0*DLX3/3.0D0-16.0D0*DLX2*DLM+8.0D0
     1*DLX*DLM2-32.0D0*ZETA2*DLX+16.0D0*DLM*S11+32.0D0*S12
     2-16.0D0*S21)+(40.0D0*X-16.0D0*X*X)*DLX2+32.0D0*X*X*DLX
     3*DLM+(280.0D0/3.0D0-704.0D0*X*X/9.0D0-88.0D0*X)*DLX
     4+(4.0D0-16.0D0*X*X/3.0D0-4.0D0*X+16.0D0/X/3.0D0)*DLM2
     5+(-208.0D0/3.0D0-64.0D0*X*X/9.0D0+160.0D0*X/3.0D0
     6+208.0D0/X/9.0D0)*DLM+(-16.0D0+64.0D0*X*X/3.0D0-16.0D0*X
     7-32.0D0/X)*ZETA2+(16.0D0-16.0D0*X+64.0D0/X/3.0D0
     8+32.0D0*X*X/3.0D0)*S11+(-32.0D0-32.0D0*X*X/3.0D0-32.0D0
     9*X-32.0D0/X/3.0D0)*(S11M+DLX*DLP)+304.0D0/9.0D0+832.0D0
     1*X*X/9.0D0-1216.0D0*X/9.0D0+80.0D0/X/9.0D0
*
      C2PS2AM0_A0  = CF * TR * A3
*
      RETURN
      END
*
************************************************************************
*     Equation (4.7): Non-singlet quark coefficient function for FL
************************************************************************
      FUNCTION CLNS2AM0_AQ(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Output Variables
*
      DOUBLE PRECISION CLNS2AM0_AQ
*
      CLNS2AM0_AQ = 16D0 * CF * TR * X / 3D0
*
      RETURN
      END
*
************************************************************************
      FUNCTION CLNS2AM0_A0(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Output Variables
*
      DOUBLE PRECISION CLNS2AM0_A0
*
      CLNS2AM0_A0 = 16D0 * CF * TR * ( X * DLOG( 1D0 - X ) 
     1            - 2D0 * X * DLOG(X) - 25D0 * X / 6D0 + 1D0 ) / 3D0
*
      RETURN
      END
*
************************************************************************
*     Equation (4.8): Light quark coefficient function for F2
*     Regular part
************************************************************************
      FUNCTION C2NS2AM0_AQ2(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Input Variables
*
      DOUBLE PRECISION B1
**
*     Output Variables
*
      DOUBLE PRECISION C2NS2AM0_AQ2
*
      B1=(-1D0-X)*2.0D0
*
      C2NS2AM0_AQ2 = 2D0 * CF * TR * B1 / 3D0
*
      RETURN
      END
*
************************************************************************
      FUNCTION C2NS2AM0_AQ(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Input Variables
*
      DOUBLE PRECISION DLX,DLM
      DOUBLE PRECISION A1
**
*     Output Variables
*
      DOUBLE PRECISION C2NS2AM0_AQ
*
      DLX=DLOG(X)
      DLM=DLOG(1.0D0-X)
*
      A1=((-8.0D0*(1.0D0+X*X)*DLX/(1.0D0-X)+13.0D0*X+1.0D0)
     1  + (-1D0-X)*(4.0D0*DLM-29.0D0/3.0D0))
*
      C2NS2AM0_AQ = 2D0 * CF * TR * A1 / 3D0
*
      RETURN
      END
*
************************************************************************
      FUNCTION C2NS2AM0_A0(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
      INCLUDE "../commons/consts.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Input Variables
*
      DOUBLE PRECISION DLX,DLX2,DLM,DLM2,SPX
      DOUBLE PRECISION A2,B2
      DOUBLE PRECISION WGPLG
**
*     Output Variables
*
      DOUBLE PRECISION C2NS2AM0_A0
*
      DLX=DLOG(X)
      DLX2=DLX*DLX
      DLM=DLOG(1.0D0-X)
      DLM2=DLM*DLM
      SPX=WGPLG(1,1,1.0D0-X)
*
      A2=(1.0D0+13.0D0*X)*DLM-(3.0D0+23.0D0*X)*DLX+29.0D0/6.0D0
     1-295.0D0*X/6.0D0+(1.0D0+X*X)*(-4.0D0*SPX-8.0D0*DLX
     1*DLM+6.0D0*DLX2+67.0D0*DLX/3.0D0)/(1.0D0-X)
*
      B2=(-1D0-X)*(-4.0D0*ZETA2+2.0D0*DLM2
     1-29.0D0*DLM/3.0D0+359D0/18.0D0)
*
      C2NS2AM0_A0 = 2D0 * CF * TR * ( A2 + B2 ) / 3D0
*
      RETURN
      END
*
************************************************************************
*     Equation (4.8): Light quark coefficient function for F2
*     Singular part
************************************************************************
      FUNCTION C2NS2BM0_AQ2(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Input Variables
*
      DOUBLE PRECISION Z
**
*     Output Variables
*
      DOUBLE PRECISION C2NS2BM0_AQ2
*
      Z = 2D0 / ( 1D0 - X )
*
      C2NS2BM0_AQ2 = Z * 4D0 * CF * TR / 3D0
*
      RETURN
      END
*
************************************************************************
      FUNCTION C2NS2BM0_AQ(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Input Variables
*
      DOUBLE PRECISION DLM,Z
**
*     Output Variables
*
      DOUBLE PRECISION C2NS2BM0_AQ
*
      DLM = DLOG( 1D0 - X )
      Z   = 2D0 / ( 1D0 - X )
*
      C2NS2BM0_AQ = Z * 2D0 * CF * TR * ( 4D0 * DLM - 29D0 / 3D0 ) / 3D0
*
      RETURN
      END
*
************************************************************************
      FUNCTION C2NS2BM0_A0(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
      INCLUDE "../commons/consts.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Input Variables
*
      DOUBLE PRECISION DLM,DLM2,Z
**
*     Output Variables
*
      DOUBLE PRECISION C2NS2BM0_A0
*
      DLM  = DLOG( 1D0 - X )
      DLM2 = DLM * DLM
      Z    = 2D0 / ( 1D0 - X )
*
      C2NS2BM0_A0 = Z * 2D0 * CF * TR * ( - 4D0 * ZETA2 + 2D0 * DLM2
     1            - 29D0 * DLM / 3D0 + 359D0 / 18D0 ) / 3D0
*
      RETURN
      END
*
************************************************************************
*     Equation (4.10): Light quark coefficient function for F2
*     Local part
************************************************************************
      FUNCTION C2NS2CM0_AQ2(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Input Variables
*
      DOUBLE PRECISION DLM
**
*     Output Variables
*
      DOUBLE PRECISION C2NS2CM0_AQ2
*
      DLM = DLOG( 1D0 - X )
*
      C2NS2CM0_AQ2 = 8D0 * CF * TR * DLM / 3D0 + 2D0 * CF * TR
*
      RETURN
      END
*
************************************************************************
      FUNCTION C2NS2CM0_AQ(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
      INCLUDE "../commons/consts.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Input Variables
*
      DOUBLE PRECISION DLM,DLM2
**
*     Output Variables
*
      DOUBLE PRECISION C2NS2CM0_AQ
*
      DLM  = DLOG( 1D0 - X )
      DLM2 = DLM * DLM
*
      C2NS2CM0_AQ = 4D0 * CF * TR * ( 2D0 * DLM2 - 29D0 * DLM / 3D0 )
     1            / 3D0 - CF * TR * ( 32D0 * ZETA2 / 3D0 + 38D0 / 3D0 )
*
      RETURN
      END
*
************************************************************************
      FUNCTION C2NS2CM0_A0(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
      INCLUDE "../commons/consts.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Input Variables
*
      DOUBLE PRECISION DLM,DLM2,DLM3
**
*     Output Variables
*
      DOUBLE PRECISION C2NS2CM0_A0
*
      DLM  = DLOG( 1D0 - X )
      DLM2 = DLM * DLM
      DLM3 = DLM2 * DLM
*
      C2NS2CM0_A0 = 4D0 * CF * TR * ( - 4D0 * ZETA2 * DLM
     1            + 2D0 * DLM3 / 3D0 - 29D0 * DLM2 / 6D0 
     2            + 359D0 * DLM / 18D0 ) / 3D0
     3            + CF * TR * ( 268D0 * ZETA2 / 9D0 + 265D0 / 9D0 )
*
      RETURN
      END
*
************************************************************************
*
*     CHARGED CURRENT
*
************************************************************************
      FUNCTION CG1ACCM0_AL(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Output Variables
*
      DOUBLE PRECISION CG1ACCM0_AL
*
      CG1ACCM0_AL = 2D0 * TR * ( X**2 + ( 1D0 - X )**2 )
*
      RETURN
      END
*
************************************************************************
*     Logarithmically divergents term to be added to the MSbar ZM
*     coefficient functions to obtain the massless limit of the IC
*     contributions.
************************************************************************
      function DICa(xi,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
      include "../commons/kfacQ.h"
**
*     Input Variables
*
      double precision xi,z
**
*     Internal Variables
*
      double precision lnk
**
*     Output Variables
*
      double precision DICa
*
      lnk = dlog(xi/kfacQ)                ! ln(muF2/mh2)
*
      DICa = 2d0 * CF * ( 1d0 + z ) * ( - lnk + 1d0
     1     + 2d0 * dlog( 1d0 - z ) )
*
      return
      end
*
************************************************************************
      function DICb(xi,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
      include "../commons/kfacQ.h"
**
*     Input Variables
*
      double precision xi,z
**
*     Internal Variables
*
      double precision lnk
**
*     Output Variables
*
      double precision DICb
*
      lnk = dlog(xi/kfacQ)                ! ln(muF2/mh2)
*
      DICb = 4d0 * CF * ( lnk - 1d0 - 2d0 * dlog( 1d0 - z ) )
     1     / ( 1d0 - z )
*
      return
      end
*
************************************************************************
      function DICc(xi,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
      include "../commons/kfacQ.h"
**
*     Input Variables
*
      double precision xi,z
**
*     Internal Variables
*
      double precision lnk
**
*     Output Variables
*
      double precision DICc
*
      lnk = dlog(xi/kfacQ)                ! ln(muF2/mh2)
*
      DICc = 4d0 * CF * ( 3d0 * lnk / 4d0 + 1d0 
     1     + ( lnk - 1d0 ) * dlog( 1d0 - z ) - dlog( 1d0 - z )**2 )
*
      return
      end
