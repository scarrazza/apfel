************************************************************************
*
*     xcoeffunc_nc.f:
*
*     x-space NC coefficient functions times x.
*
************************************************************************
      FUNCTION XC2G_NC(X)
*
      IMPLICIT NONE
*
      include "../commons/vfns.h"
      include "../commons/scale.h"
      include "../commons/hmass.h"
      include "../commons/ipt.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION EPS,V,A,W2,Q2
      DOUBLE PRECISION C2G2A
      DOUBLE PRECISION CH_FFNS_NC_AB,C2GH,C2GHBAR
      DOUBLE PRECISION C2GM,C2GCA,C2GCF
**
*     Output Variables
*
      DOUBLE PRECISION XC2G_NC
*
      IF(VFNS.EQ."ZMVN")THEN
         IF(IPT.EQ.1)THEN
            XC2G_NC = ( ( 1D0 - X )**2D0 + X**2D0 ) 
     1              * DLOG( ( 1D0 - X ) / X ) 
     2              - 8D0 * X**2D0 + 8D0 * X - 1D0
         ELSEIF(IPT.EQ.2)THEN
            XC2G_NC = C2G2A(X,1)
         ENDIF
      ELSEIF(VFNS.EQ."FFNS")THEN
         Q2 = Q * Q
         EPS = Q2H / Q2
         W2 = Q2 * ( 1D0 - X ) / X
         V = DSQRT( 1D0 - 4D0 * Q2H / W2 )
         A = 1D0 / ( 1D0 + 4D0 * EPS )
*
*     Massive scheme coeffcient function
*
         IF(X.GT.A)THEN
            XC2G_NC = 0D0
            RETURN
         ENDIF
         IF(IPT.EQ.1)THEN
            XC2G_NC = ( X**2D0 + ( 1D0 - X )**2D0 
     1              + 4D0 * EPS * X * ( 1D0 - 3D0 * X ) 
     2              - 8D0 * EPS**2D0 * X**2D0 ) 
     3              * DLOG( ( 1D0 + V ) / ( 1D0 -V ) )
     4              + ( 8D0 * X * (1D0 - X ) - 1D0 
     5              - 4D0 * EPS * X * ( 1D0 - X ) ) * V
         ELSEIF(IPT.EQ.2)THEN
            C2GH    = CH_FFNS_NC_AB(3,X,Q2,Q2H)
            C2GHBAR = CH_FFNS_NC_AB(9,X,Q2,Q2H)
*
            XC2G_NC = C2GH + C2GHBAR * DLOG(1D0/EPS)
         ENDIF
      ELSEIF(VFNS.EQ."FFN0")THEN
         Q2 = Q * Q
         IF(IPT.EQ.1)THEN
            XC2G_NC = ( ( 1D0 - X )**2D0 + X**2D0 ) 
     1              * ( DLOG( ( 1D0 - X ) / X ) + DLOG( Q2 / Q2H ) ) 
     2              - 8D0 * X**2D0 + 8D0 * X - 1D0
         ELSEIF(IPT.EQ.2)THEN
            XC2G_NC = C2GM(X,Q2,Q2,Q2H) + C2GCA(X,Q2,Q2H) 
     1              + C2GCF(X,Q2,Q2H)
         ENDIF
      ENDIF
      XC2G_NC = X * XC2G_NC
*
      RETURN
*
      END
*
C======================================================================
      FUNCTION XC2Q_UNPLUS_NC(X)
*
      IMPLICIT NONE
*
      include "../commons/vfns.h"
      include "../commons/scale.h"
      include "../commons/hmass.h"
      include "../commons/ipt.h"
      include "../commons/nf.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION C2NN2A!,C2NN2A_DNF
      DOUBLE PRECISION Q2
      DOUBLE PRECISION CH_FFNS_NC_AB,C2QH
      DOUBLE PRECISION C2Q2R
**
*     Output Variables
*
      DOUBLE PRECISION XC2Q_UNPLUS_NC
*
      IF(VFNS.EQ."ZMVN")THEN
         IF(IPT.EQ.1)THEN
            XC2Q_UNPLUS_NC = 8D0 * ( - ( 1D0 + X ) * DLOG( 1D0 - X )
     1                       - ( 1D0 + X**2D0 ) * DLOG(X) / ( 1D0 -  X ) 
     2                       + 3D0 + 2D0 * X ) / 3D0
         ELSEIF(IPT.EQ.2)THEN
            XC2Q_UNPLUS_NC = C2NN2A(X,NF)
         ENDIF
      ELSEIF(VFNS.EQ."FFNS")THEN
         IF(IPT.EQ.2)THEN
            Q2   = Q * Q
            C2QH = CH_FFNS_NC_AB(7,X,Q2,Q2H)
*
            XC2Q_UNPLUS_NC = C2QH !- C2NN2A_DNF(X)
         ELSE
            XC2Q_UNPLUS_NC = 0D0
         ENDIF
      ELSEIF(VFNS.EQ."FFN0")THEN
         IF(IPT.EQ.2)THEN
            Q2   = Q * Q
            XC2Q_UNPLUS_NC = C2Q2R(X,Q2,Q2H)
         ELSE
            XC2Q_UNPLUS_NC = 0D0
         ENDIF
      ENDIF
      XC2Q_UNPLUS_NC = X * XC2Q_UNPLUS_NC
*
      RETURN
*
      END
*
C======================================================================
      FUNCTION XC2Q_UNPLUS_PS_NC(X)
*
      IMPLICIT NONE
*
      include "../commons/vfns.h"
      include "../commons/scale.h"
      include "../commons/hmass.h"
      include "../commons/ipt.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION C2S2A
      DOUBLE PRECISION Q2,EPS
      DOUBLE PRECISION CH_FFNS_NC_AB,C2QH,C2QHBAR
      DOUBLE PRECISION C2Q1M,C2Q1
**
*     Output Variables
*
      DOUBLE PRECISION XC2Q_UNPLUS_PS_NC
*
      IF(VFNS.EQ."ZMVN")THEN
         IF(IPT.EQ.2)THEN
            XC2Q_UNPLUS_PS_NC = C2S2A(X,1)
         ELSE
            XC2Q_UNPLUS_PS_NC = 0D0
         ENDIF
      ELSEIF(VFNS.EQ."FFNS")THEN
         IF(IPT.EQ.2)THEN
            Q2 = Q * Q
            EPS = Q2H / Q2
*
            C2QH    = CH_FFNS_NC_AB(4,X,Q2,Q2H)
            C2QHBAR = CH_FFNS_NC_AB(10,X,Q2,Q2H)
*
            XC2Q_UNPLUS_PS_NC = C2QH + C2QHBAR * DLOG(1D0/EPS)
         ELSE
            XC2Q_UNPLUS_PS_NC = 0D0
         ENDIF
      ELSEIF(VFNS.EQ."FFN0")THEN
         IF(IPT.EQ.2)THEN
            Q2 = Q * Q
            XC2Q_UNPLUS_PS_NC = C2Q1M(X,Q2,Q2,Q2H) + C2Q1(X,Q2,Q2H)
         ELSE
            XC2Q_UNPLUS_PS_NC = 0D0
         ENDIF
      ENDIF
      XC2Q_UNPLUS_PS_NC = X * XC2Q_UNPLUS_PS_NC
*
      RETURN
*
      END
*
C======================================================================
      FUNCTION XC2Q_PLUS_NC(X)
*
      IMPLICIT NONE
*
      include "../commons/vfns.h"
      include "../commons/scale.h"
      include "../commons/hmass.h"
      include "../commons/ipt.h"
      include "../commons/nf.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION C2NS2B
      DOUBLE PRECISION Q2,C2Q2S
**
*     Output Variables
*
      DOUBLE PRECISION XC2Q_PLUS_NC
*
      IF(IPT.EQ.1)THEN
         XC2Q_PLUS_NC = 8D0 * ( 2D0 * DLOG( 1D0 - X ) / ( 1D0 - X )
     1                - 3D0 / 2D0 / ( 1D0 - X ) ) / 3D0
      ELSEIF(IPT.EQ.2)THEN
         IF(VFNS.EQ."FFN0")THEN
            Q2 = Q * Q
            XC2Q_PLUS_NC = C2Q2S(X,Q2,Q2H)
         ELSE
            XC2Q_PLUS_NC = C2NS2B(X,NF)
         ENDIF
      ENDIF
      XC2Q_PLUS_NC = X * XC2Q_PLUS_NC
*
      RETURN
*
      END
*
C======================================================================
      FUNCTION XC3G_NC(X)
*
      IMPLICIT NONE
*
      include "../commons/ipt.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Output Variables
*
      DOUBLE PRECISION XC3G_NC
*
      XC3G_NC = 0D0
*
      RETURN
*
      END
*
C======================================================================
      FUNCTION XC3Q_UNPLUS_NC(X)
*
      IMPLICIT NONE
*
      include "../commons/ipt.h"
      include "../commons/nf.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION C3NM2A
**
*     Output Variables
*
      DOUBLE PRECISION XC3Q_UNPLUS_NC
*
      IF(IPT.EQ.1)THEN
         XC3Q_UNPLUS_NC = 8D0 * ( - ( 1D0 + X ) * DLOG( 1D0 - X )
     1                  - ( 1D0 + X**2D0 ) * DLOG(X) / ( 1D0 -  X ) 
     2                  + 2D0 + X ) / 3D0
      ELSEIF(IPT.EQ.2)THEN
         XC3Q_UNPLUS_NC = C3NM2A(X,NF)
      ENDIF
      XC3Q_UNPLUS_NC = X * XC3Q_UNPLUS_NC
*
      RETURN
*
      END
*
C======================================================================
      FUNCTION XC3Q_PLUS_NC(X)
*
      IMPLICIT NONE
*
      include "../commons/ipt.h"
      include "../commons/nf.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION C3NS2B
**
*     Output Variables
*
      DOUBLE PRECISION XC3Q_PLUS_NC
*
      IF(IPT.EQ.1)THEN
         XC3Q_PLUS_NC = 8D0 + ( 2D0 * DLOG( 1D0 - X ) / ( 1D0 - X )
     1                - 3D0 / 2D0 / ( 1D0 - X ) ) / 3D0
      ELSEIF(IPT.EQ.2)THEN
         XC3Q_PLUS_NC = C3NS2B(X,NF)
      ENDIF
      XC3Q_PLUS_NC = X * XC3Q_PLUS_NC
*
      RETURN
*
      END
*
C======================================================================
      FUNCTION XCLG_NC(X)
*
      IMPLICIT NONE
*
      include "../commons/vfns.h"
      include "../commons/scale.h"
      include "../commons/hmass.h"
      include "../commons/ipt.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION EPS,V,A,W2,Q2
      DOUBLE PRECISION CLG2A
      DOUBLE PRECISION CH_FFNS_NC_AB,CLGH,CLGHBAR
      DOUBLE PRECISION CLGM,CLGCA,CLGCF
**
*     Output Variables
*
      DOUBLE PRECISION XCLG_NC
*
      IF(VFNS.EQ."ZMVN")THEN
         IF(IPT.EQ.1)THEN
            XCLG_NC = 4D0 * X * ( 1D0 - X )
         ELSEIF(IPT.EQ.2)THEN
            XCLG_NC = CLG2A(X,1)
         ENDIF
      ELSEIF(VFNS.EQ."FFNS")THEN
         Q2 = Q * Q
         EPS = Q2H / Q2
         W2 = Q2 * ( 1D0 - X ) / X
         V = DSQRT( 1D0 - 4D0 * Q2H / W2 )
         A = 1D0 / ( 1D0 + 4D0 * EPS)
*
*     Massive scheme coefficient function
*
         IF(X.GT.A)THEN
            XCLG_NC = 0D0
            RETURN
         ENDIF
         IF(IPT.EQ.1)THEN
            XCLG_NC = - 8D0 * EPS * X**2D0 
     1              * DLOG( ( 1D0 + V ) / ( 1D0 - V ) )
     2              + 4D0 * V * X * ( 1D0 - X )
         ELSEIF(IPT.EQ.2)THEN
            CLGH    = CH_FFNS_NC_AB(5,X,Q2,Q2H)
            CLGHBAR = CH_FFNS_NC_AB(11,X,Q2,Q2H)
*
            XCLG_NC = CLGH + CLGHBAR * DLOG(1D0/EPS)
         ENDIF
      ELSEIF(VFNS.EQ."FFN0")THEN
         IF(IPT.EQ.1)THEN
            XCLG_NC = 4D0 * X * ( 1D0 - X )
         ELSEIF(IPT.EQ.2)THEN
            Q2 = Q * Q
            XCLG_NC = CLGM(X,Q2,Q2,Q2H) + CLGCA(X,Q2,Q2H)
     1              + CLGCF(X,Q2,Q2H)
         ENDIF
      ENDIF
      XCLG_NC = X * XCLG_NC
*
      RETURN
*
      END
*
C======================================================================
      FUNCTION XCLQ_UNPLUS_NC(X)
*
      IMPLICIT NONE
*
      include "../commons/vfns.h"
      include "../commons/scale.h"
      include "../commons/hmass.h"
      include "../commons/ipt.h"
      include "../commons/nf.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION CLNN2A
      DOUBLE PRECISION Q2
      DOUBLE PRECISION CH_FFNS_NC_AB,CLQH
      DOUBLE PRECISION CLQ2
**
*     Output Variables
*
      DOUBLE PRECISION XCLQ_UNPLUS_NC
*
      IF(VFNS.EQ."ZMVN")THEN
         IF(IPT.EQ.1)THEN
            XCLQ_UNPLUS_NC = 16D0 * X / 3D0
         ELSEIF(IPT.EQ.2)THEN
            XCLQ_UNPLUS_NC = CLNN2A(X,NF)
         ENDIF
      ELSEIF(VFNS.EQ."FFNS")THEN
         IF(IPT.EQ.2)THEN
            Q2   = Q * Q
            CLQH = CH_FFNS_NC_AB(8,X,Q2,Q2H)
*
            XCLQ_UNPLUS_NC = CLQH
         ELSE
            XCLQ_UNPLUS_NC = 0D0
         ENDIF
      ELSEIF(VFNS.EQ."FFN0")THEN
         IF(IPT.EQ.2)THEN
            Q2   = Q * Q
            XCLQ_UNPLUS_NC = CLQ2(X,Q2,Q2H)
         ELSE
            XCLQ_UNPLUS_NC = 0D0
         ENDIF
      ENDIF
      XCLQ_UNPLUS_NC = X * XCLQ_UNPLUS_NC
*
      RETURN
*
      END
*
C======================================================================
      FUNCTION XCLQ_UNPLUS_PS_NC(X)
*
      IMPLICIT NONE
*
      include "../commons/vfns.h"
      include "../commons/scale.h"
      include "../commons/hmass.h"
      include "../commons/ipt.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION CLS2A
      DOUBLE PRECISION Q2,EPS
      DOUBLE PRECISION CH_FFNS_NC_AB,CLQH,CLQHBAR
      DOUBLE PRECISION CLQ1M,CLQ1
**
*     Output Variables
*
      DOUBLE PRECISION XCLQ_UNPLUS_PS_NC
*
*
      IF(VFNS.EQ."ZMVN")THEN
         IF(IPT.EQ.2)THEN
            XCLQ_UNPLUS_PS_NC = CLS2A(X,1)
         ELSE
            XCLQ_UNPLUS_PS_NC = 0D0
         ENDIF
      ELSEIF(VFNS.EQ."FFNS")THEN
         IF(IPT.EQ.2)THEN
            Q2 = Q * Q
            EPS = Q2H / Q2
*
            CLQH    = CH_FFNS_NC_AB(6,X,Q2,Q2H)
            CLQHBAR = CH_FFNS_NC_AB(12,X,Q2,Q2H)
*
            XCLQ_UNPLUS_PS_NC = CLQH + CLQHBAR * DLOG(1D0/EPS)
         ELSE
            XCLQ_UNPLUS_PS_NC = 0D0
         ENDIF
      ELSEIF(VFNS.EQ."FFN0")THEN
         IF(IPT.EQ.2)THEN
            Q2 = Q * Q
            XCLQ_UNPLUS_PS_NC = CLQ1M(X,Q2,Q2,Q2H) + CLQ1(X,Q2,Q2H)
         ELSE
            XCLQ_UNPLUS_PS_NC = 0D0
         ENDIF

      ENDIF
      XCLQ_UNPLUS_PS_NC = X * XCLQ_UNPLUS_PS_NC
*
      RETURN
*
      END
*
C======================================================================
      FUNCTION XCLQ_PLUS_NC(X)
*
      IMPLICIT NONE
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Output Variables
*
      DOUBLE PRECISION XCLQ_PLUS_NC
*
      XCLQ_PLUS_NC = 0D0
*
      RETURN
*
      END
C=======================================================================
C
C     Wrap function of one argument (y) needed to be integrated with
C     DGAUSS
C
C=======================================================================
      FUNCTION C2Q_PLUS_NC_WRAP(Y)
*
      IMPLICIT NONE
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION XC2Q_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION C2Q_PLUS_NC_WRAP
*
      C2Q_PLUS_NC_WRAP = XC2Q_PLUS_NC(Y) / Y
*
      RETURN
      END
*
C=======================================================================
      FUNCTION C3Q_PLUS_NC_WRAP(Y)
*
      IMPLICIT NONE
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION XC3Q_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION C3Q_PLUS_NC_WRAP
*
      C3Q_PLUS_NC_WRAP = XC3Q_PLUS_NC(Y) / Y
*
      RETURN
      END
*
C=======================================================================
      FUNCTION CLQ_PLUS_NC_WRAP(Y)
*
      IMPLICIT NONE
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION XCLQ_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION CLQ_PLUS_NC_WRAP
*
      CLQ_PLUS_NC_WRAP = XCLQ_PLUS_NC(Y) / Y
*
      RETURN
      END
*
C=======================================================================
