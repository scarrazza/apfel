************************************************************************
*
*     xcoeffunc_cc.f:
*
*     x-space CC coefficient functions times x.
*
************************************************************************
      FUNCTION XC1G(X)
*
      IMPLICIT NONE
*
      include "../commons/vfns.h"
      include "../commons/hmass.h"
      include "../commons/scale.h"
      include "../commons/jpt.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION LAMBDA,Q2,PQG
      DOUBLE PRECISION C2G2A,CLG2A
**
*     Output Variables
*
      DOUBLE PRECISION XC1G
*
      IF(VFNS.EQ."ZMVN")THEN
         IF(IPT.EQ.1)THEN
            XC1G = X * (1D0 / 2D0 * ( ( ( 1D0 - X )**2D0 + X**2D0 ) 
     1           * DLOG((1D0-X)/X) - 4D0 * X**2D0 + 4D0 * X - 1D0 ) )
         ELSEIF(IPT.EQ.2)THEN
            XC1G = X * ( C2G2A(X,1) - CLG2A(X,1) )
         ENDIF
      ELSEIF(VFNS.EQ."FFNS")THEN
         IF(IPT.EQ.1)THEN
            Q2 = Q * Q
            LAMBDA = 1D0 / ( 1D0 + Q2H / Q2 )
            PQG = 1D0 / 2D0 * ( X**2D0 + ( 1D0 - X )**2D0 )
*
*     Massive scheme coeffcient function
*
            XC1G = PQG * ( - DLOG( LAMBDA * ( 1D0 - LAMBDA ) )
     1           + 2D0 * DLOG( ( 1D0 - X ) / X ) )
     2           + ( 4D0 - 4D0 * ( 1D0 - LAMBDA ) ) * X * ( 1D0 - X )
     3           + ( 1D0 - LAMBDA ) * X / ( 1D0 - LAMBDA * X) - 1D0
     4           + 2D0 * ( 1D0 - LAMBDA ) * X 
     5           * ( 1D0 - 2D0 * LAMBDA * X ) 
     6           * DLOG( ( 1D0 - LAMBDA * X ) / ( ( 1D0 - LAMBDA) * X ))
*     
*     Uncomment below to subtract the colliner divergence
*
C            XC1G =  XC1G - PQG * DLOG(Q2/Q2H)  
*
            XC1G = X * XC1G / 2D0
         ELSE
            XC1G = 0D0
         ENDIF
      ELSEIF(VFNS.EQ."FFN0")THEN
         IF(IPT.EQ.1)THEN
            Q2 = Q * Q
            XC1G = ( ( 1D0 - X )**2D0 + X**2D0 ) 
     1           * ( DLOG( ( 1D0 - X ) / X ) + DLOG( Q2 / Q2H ) / 2D0 ) 
     2           - 4D0 * X**2D0 + 4D0 * X - 1D0
            XC1G = X * XC1G / 2D0
         ELSE
            XC1G = 0D0
         ENDIF
      ENDIF

      RETURN
*
      END
*
************************************************************************
      FUNCTION XC1Q_UNPLUS(X)
*
      IMPLICIT NONE
*
      include "../commons/vfns.h"
      include "../commons/hmass.h"
      include "../commons/scale.h"
      include "../commons/jpt.h"
      include "../commons/nf.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION LAMBDA,Q2
      DOUBLE PRECISION C2NN2A,CLNN2A
**
*     Output Variables
*
      DOUBLE PRECISION XC1Q_UNPLUS
*
      IF(VFNS.EQ."ZMVN")THEN
         IF(IPT.EQ.1)THEN
            XC1Q_UNPLUS = X * ( 4D0/3D0 * ( - (1D0 + X) * DLOG(1D0 - X )
     1                  - (1D0 + X**2D0) * DLOG(X) / (1D0-X) + 3D0 ))
         ELSEIF(IPT.EQ.2)THEN
            XC1Q_UNPLUS = X * ( C2NN2A(X,NF) - CLNN2A(X,NF) )
         ENDIF
      ELSEIF(VFNS.EQ."FFNS")THEN
         IF(IPT.EQ.1)THEN
            Q2 = Q * Q
            LAMBDA = 1D0 / ( 1D0 + Q2H / Q2 )
*
*     Massive scheme coefficient function
*
            XC1Q_UNPLUS = 4D0 / 3D0 * ( - (1D0 + X**2D0) * DLOG(X) 
     1                  / ( 1D0 -  X ) - 2D0 * (1D0 + X) * DLOG(1D0-X)
     2                  + (1D0 + X) * DLOG(1D0-LAMBDA*X) + ( 3D0 - X )
     3                  + X * ( 1D0 - X ) / ( 1D0 - LAMBDA * X ) )
            XC1Q_UNPLUS = X * XC1Q_UNPLUS
         ELSE
            XC1Q_UNPLUS = 0D0
         ENDIF
      ELSEIF(VFNS.EQ."FFN0")THEN
         IF(IPT.EQ.1)THEN
            XC1Q_UNPLUS = X * ( 4D0/3D0 * ( - (1D0 + X) * DLOG(1D0 - X)
     1                  - (1D0 + X**2D0) * DLOG(X) / (1D0-X) + 3D0 ))
         ELSE
            XC1Q_UNPLUS = 0D0
         ENDIF
      ENDIF
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XC1Q_UNPLUS_PS(X)
*
      IMPLICIT NONE
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION C2S2A,CLS2A
**
*     Output Variables
*
      DOUBLE PRECISION XC1Q_UNPLUS_PS
*
      XC1Q_UNPLUS_PS = X * ( C2S2A(X,1) - CLS2A(X,1) )
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XC1Q_PLUS(X)
*
      IMPLICIT NONE
*
      include "../commons/vfns.h"
      include "../commons/hmass.h"
      include "../commons/scale.h"
      include "../commons/jpt.h"
      include "../commons/nf.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION LAMBDA,Q2,PQQ
      DOUBLE PRECISION C2NS2B
**
*     Output Variables
*
      DOUBLE PRECISION XC1Q_PLUS
*
      IF(VFNS.EQ."ZMVN")THEN
         IF(IPT.EQ.1)THEN
            XC1Q_PLUS = X * ( 4D0 / 3D0 * ( 2D0 * DLOG(1D0-X)/(1D0-X)
     1                - 3D0 / 2D0 / ( 1D0 - X ) ) )
         ELSEIF(IPT.EQ.2)THEN
            XC1Q_PLUS = X * C2NS2B(X,NF)
         ENDIF
      ELSEIF(VFNS.EQ."FFNS")THEN
         IF(IPT.EQ.1)THEN
            Q2 = Q * Q
            LAMBDA = 1D0 / ( 1D0 + Q2H / Q2 )
            PQQ = 4D0 / 3D0 * ( 1D0 + X**2D0) / ( 1D0 - X )
*
*     Massive scheme coefficient function
*
            XC1Q_PLUS = - PQQ * DLOG(LAMBDA) + 4D0 / 3D0 
     1                * ( 4D0 * DLOG( 1D0 - X ) / ( 1D0 - X ) 
     2                - 2D0 * DLOG( 1D0 - LAMBDA * X ) / ( 1D0 - X )
     3                - 2D0 / ( 1D0 - X ) + 1D0 / 2D0 * ( ( 1D0 - X ) 
     4                / ( 1D0 - LAMBDA * X )**2D0 ) )
            XC1Q_PLUS = X * XC1Q_PLUS
         ELSE
            XC1Q_PLUS = 0D0
         ENDIF
      ELSEIF(VFNS.EQ."FFN0")THEN
         IF(IPT.EQ.1)THEN
            XC1Q_PLUS = X * ( 4D0 / 3D0 * ( 2D0 * DLOG(1D0-X)/(1D0-X)
     1                - 3D0 / 2D0 / ( 1D0 - X ) ) )
         ELSE
            XC1Q_PLUS = 0D0
         ENDIF
      ENDIF
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XC2G(X)
*
      IMPLICIT NONE
*
      include "../commons/vfns.h"
      include "../commons/hmass.h"
      include "../commons/scale.h"
      include "../commons/jpt.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION LAMBDA,Q2,PQG
      DOUBLE PRECISION C2G2A
**
*     Output Variables
*
      DOUBLE PRECISION XC2G
*
      IF(VFNS.EQ."ZMVN")THEN
         IF(IPT.EQ.1)THEN
            XC2G = X * (1D0 / 2D0 * ( ( ( 1D0 - X )**2D0 + X**2D0 ) 
     1           * DLOG((1D0-X)/X) - 8D0 * X**2D0 + 8D0 * X - 1D0 ))
         ELSEIF(IPT.EQ.2)THEN
            XC2G = X * C2G2A(X,1)
         ENDIF
      ELSEIF(VFNS.EQ."FFNS")THEN
         IF(IPT.EQ.1)THEN
            Q2 = Q * Q
            LAMBDA = 1D0 / ( 1D0 + Q2H / Q2 )
            PQG = 1D0 / 2D0 * ( X**2D0 + ( 1D0 - X )**2D0 )
*
*     Massive scheme coeffcient function
*
            XC2G = PQG * ( - DLOG( LAMBDA * ( 1D0 - LAMBDA ) )
     1           + 2D0 * DLOG( ( 1D0 - X ) / X ) )
     2           + ( 8D0 - 18D0 * ( 1D0 - LAMBDA ) 
     3           + 12D0 * ( 1D0 - LAMBDA )**2D0 ) * X * ( 1D0 - X )
     4           + ( 1D0 - LAMBDA ) / ( 1D0 - LAMBDA * X) - 1D0
     5           + 6D0 * LAMBDA * ( 1D0 - LAMBDA ) * X 
     6           * ( 1D0 - 2D0 * LAMBDA * X ) 
     7           * DLOG( ( 1D0 - LAMBDA * X ) / ( ( 1D0 - LAMBDA) * X ))
*
*     Uncomment below to subtract the colliner divergence
*
C            XC2G =  XC2G - PQG * DLOG( Q2 / Q2H )  
*
            XC2G = X * XC2G / 2D0
         ELSE
            XC2G = 0D0
         ENDIF
      ELSEIF(VFNS.EQ."FFN0")THEN
         IF(IPT.EQ.1)THEN
            Q2 = Q * Q
            XC2G = ( ( 1D0 - X )**2D0 + X**2D0 ) 
     1           * ( DLOG( ( 1D0 - X ) / X ) + DLOG( Q2 / Q2H ) / 2D0 ) 
     2           - 8D0 * X**2D0 + 8D0 * X - 1D0
            XC2G = X * XC2G / 2D0
         ELSE
            XC2G = 0D0
         ENDIF
      ENDIF
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XC2Q_UNPLUS(X)
*
      IMPLICIT NONE
*
      include "../commons/vfns.h"
      include "../commons/hmass.h"
      include "../commons/scale.h"
      include "../commons/jpt.h"
      include "../commons/nf.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION LAMBDA,Q2
      DOUBLE PRECISION C2NN2A
**
*     Output Variables
*
      DOUBLE PRECISION XC2Q_UNPLUS
*
      IF(VFNS.EQ."ZMVN")THEN
         IF(IPT.EQ.1)THEN
            XC2Q_UNPLUS = X * (4D0/3D0 * ( - (1D0 + X) * DLOG(1D0 - X)
     1                  - (1D0 + X**2D0) * DLOG(X) / ( 1D0 -  X ) + 3D0
     2                  + 2D0 * X ) )
         ELSEIF(IPT.EQ.2)THEN
            XC2Q_UNPLUS = X * C2NN2A(X,NF)
         ENDIF
      ELSEIF(VFNS.EQ."FFNS")THEN
         IF(IPT.EQ.1)THEN
            Q2 = Q * Q
            LAMBDA = 1D0 / ( 1D0 + Q2H / Q2 )
*
*     Massive scheme coefficient function
*
            XC2Q_UNPLUS = 4D0 / 3D0 * ( - (1D0 + X**2D0) * DLOG(X) 
     1                  / ( 1D0 -  X ) - 2D0 * (1D0 + X) * DLOG(1D0-X)
     2                  + (1D0 + X) * DLOG(1D0-LAMBDA*X) + ( 2D0 * X 
     3                  + 2D0 - 2D0 / X ) + ( 2D0 / X - 1D0 - X )
     4                  / ( 1D0 - LAMBDA * X ))
            XC2Q_UNPLUS = X * XC2Q_UNPLUS
         ELSE
            XC2Q_UNPLUS = 0D0
         ENDIF
      ELSEIF(VFNS.EQ."FFN0")THEN
         IF(IPT.EQ.1)THEN
            XC2Q_UNPLUS = X * (4D0/3D0 * ( - (1D0 + X) * DLOG(1D0 - X)
     1                  - (1D0 + X**2D0) * DLOG(X) / ( 1D0 -  X ) + 3D0
     2                  + 2D0 * X ) )
         ELSE
            XC2Q_UNPLUS = 0D0
         ENDIF
      ENDIF
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XC2Q_UNPLUS_PS(X)
*
      IMPLICIT NONE
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION C2S2A
**
*     Output Variables
*
      DOUBLE PRECISION XC2Q_UNPLUS_PS
*
      XC2Q_UNPLUS_PS = X * C2S2A(X,1)
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XC2Q_PLUS(X)
*
      IMPLICIT NONE
*
      include "../commons/vfns.h"
      include "../commons/hmass.h"
      include "../commons/scale.h"
      include "../commons/jpt.h"
      include "../commons/nf.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION LAMBDA,Q2,PQQ
      DOUBLE PRECISION C2NS2B
**
*     Output Variables
*
      DOUBLE PRECISION XC2Q_PLUS
*
      IF(VFNS.EQ."ZMVN")THEN
         IF(IPT.EQ.1)THEN
            XC2Q_PLUS = X * ( 4D0 / 3D0 * ( 2D0 * DLOG(1D0-X)/(1D0-X)
     1                - 3D0 / 2D0 / ( 1D0 - X ) ) )
         ELSEIF(IPT.EQ.2)THEN
            XC2Q_PLUS = X * C2NS2B(X,NF)
         ENDIF
      ELSEIF(VFNS.EQ."FFNS")THEN
         IF(IPT.EQ.1)THEN
            Q2 = Q * Q
            LAMBDA = 1D0 / ( 1D0 + Q2H / Q2 )
            PQQ = 4D0 / 3D0 * ( 1D0 + X**2D0) / ( 1D0 - X )
*
*     Massive scheme coefficient function
*
            XC2Q_PLUS = - PQQ * DLOG(LAMBDA) + 4D0 / 3D0 
     1                * ( 4D0 * DLOG( 1D0 - X ) / ( 1D0 - X ) 
     2                - 2D0 * DLOG( 1D0 - LAMBDA * X ) / ( 1D0 - X )
     3                - 2D0 / ( 1D0 - X ) + 1D0 / 2D0 * ( ( 1D0 - X ) 
     4                / ( 1D0 - LAMBDA * X )**2D0 ) )
            XC2Q_PLUS = X * XC2Q_PLUS
         ELSE
            XC2Q_PLUS = 0D0
         ENDIF
      ELSEIF(VFNS.EQ."FFN0")THEN
         IF(IPT.EQ.1)THEN
            XC2Q_PLUS = X * ( 4D0 / 3D0 * ( 2D0 * DLOG(1D0-X)/(1D0-X)
     1                - 3D0 / 2D0 / ( 1D0 - X ) ) )
         ELSE
            XC2Q_PLUS = 0D0
         ENDIF
      ENDIF
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XC3G(X)
*
      IMPLICIT NONE
*
      include "../commons/vfns.h"
      include "../commons/hmass.h"
      include "../commons/scale.h"
      include "../commons/jpt.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION LAMBDA,Q2,PQG
**
*     Output Variables
*
      DOUBLE PRECISION XC3G
*
      IF(VFNS.EQ."ZMVN")THEN
         XC3G = 0D0
      ELSEIF(VFNS.EQ."FFNS")THEN
         IF(IPT.EQ.1)THEN
            Q2 = Q * Q
            LAMBDA = 1D0 / ( 1D0 + Q2H / Q2 )
            PQG = 1D0 / 2D0 * ( X**2D0 + ( 1D0 - X )**2D0 )
*
*     Massive scheme coeffcient function
*
            XC3G = PQG * ( DLOG( ( 1D0 - LAMBDA ) / LAMBDA )
     1           + 2D0 * DLOG( ( 1D0 - X ) / ( 1D0 - LAMBDA * X ) ) )
     2           + 2D0 * ( 1D0 - LAMBDA ) * X * ( 1D0 - X )
     3           + 2D0 * ( 1D0 - LAMBDA ) * X 
     4           * ( ( 1D0 + LAMBDA ) * X - 1D0 ) 
     5           * DLOG( ( 1D0 - LAMBDA * X ) / ( (1D0 - LAMBDA) * X ) )
*
*     Uncomment below to subtract the colliner divergence
*
C            XC3G =  XC3G + PQG * DLOG( Q2 / Q2H )  
*
            XC3G = X * XC3G / 2D0
         ELSE
            XC3G = 0D0
         ENDIF
      ELSEIF(VFNS.EQ."FFN0")THEN
         IF(IPT.EQ.1)THEN
            Q2 = Q * Q
            XC3G = - 1D0 / 2D0 * ( X**2D0 + ( 1D0 - X )**2D0 )
     1           * DLOG( Q2 / Q2H )
            XC3G = X * XC3G / 2D0
         ELSE
            XC3G = 0D0
         ENDIF
      ENDIF

      RETURN
*
      END
*
************************************************************************
      FUNCTION XC3Q_UNPLUS(X)
*
      IMPLICIT NONE
*
      include "../commons/vfns.h"
      include "../commons/hmass.h"
      include "../commons/scale.h"
      include "../commons/jpt.h"
      include "../commons/nf.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION LAMBDA,Q2
      DOUBLE PRECISION C3NM2A
**
*     Output Variables
*
      DOUBLE PRECISION XC3Q_UNPLUS
*
      IF(VFNS.EQ."ZMVN")THEN
         IF(IPT.EQ.1)THEN
            XC3Q_UNPLUS = X * ( 4D0/3D0 * ( - (1D0 + X) * DLOG(1D0 - X)
     1                  - ( 1D0 + X**2D0 ) * DLOG(X) / (1D0 - X ) 
     2                  + ( X + 2D0 ) ) )
         ELSEIF(IPT.EQ.2)THEN
            XC3Q_UNPLUS = X * C3NM2A(X,NF)
         ENDIF
      ELSEIF(VFNS.EQ."FFNS")THEN
         IF(IPT.EQ.1)THEN
            Q2 = Q * Q
            LAMBDA = 1D0 / ( 1D0 + Q2H / Q2 )
*
*     Massive scheme coefficient function
*
            XC3Q_UNPLUS = 4D0 / 3D0 * ( - (1D0 + X**2D0) * DLOG(X) 
     1                  / ( 1D0 -  X ) - 2D0 * (1D0 + X) * DLOG(1D0-X)
     2                  + (1D0 + X) * DLOG(1D0-LAMBDA*X) + ( 1D0 + X )
     3                  + ( 1D0 - X ) / ( 1D0 - LAMBDA * X ) )
            XC3Q_UNPLUS = X * XC3Q_UNPLUS
         ELSE
            XC3Q_UNPLUS = 0D0
         ENDIF
      ELSEIF(VFNS.EQ."FFN0")THEN
         IF(IPT.EQ.1)THEN
            XC3Q_UNPLUS = X * ( 4D0/3D0 * ( - (1D0 + X) * DLOG(1D0 - X)
     1                  - ( 1D0 + X**2D0 ) * DLOG(X) / (1D0 - X ) 
     2                  + ( X + 2D0 ) ) )
         ELSE
            XC3Q_UNPLUS = 0D0
         ENDIF
      ENDIF
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XC3Q_PLUS(X)
*
      IMPLICIT NONE
*
      include "../commons/vfns.h"
      include "../commons/hmass.h"
      include "../commons/scale.h"
      include "../commons/jpt.h"
      include "../commons/nf.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION LAMBDA,Q2,PQQ
      DOUBLE PRECISION C3NS2B
**
*     Output Variables
*
      DOUBLE PRECISION XC3Q_PLUS
*
      IF(VFNS.EQ."ZMVN")THEN
         IF(IPT.EQ.1)THEN
            XC3Q_PLUS = X * ( 4D0 / 3D0 * ( 2D0 * DLOG(1D0-X)/(1D0-X)
     1                - 3D0 / 2D0 / ( 1D0 - X ) ) )
         ELSEIF(IPT.EQ.2)THEN
            XC3Q_PLUS = X * C3NS2B(X,NF)
         ENDIF
      ELSEIF(VFNS.EQ."FFNS")THEN
         IF(IPT.EQ.1)THEN
            Q2 = Q * Q
            LAMBDA = 1D0 / ( 1D0 + Q2H / Q2 )
            PQQ = 4D0 / 3D0 * ( 1D0 + X**2D0) / ( 1D0 - X )
*
*     Massive scheme coefficient function
*
            XC3Q_PLUS = - PQQ * DLOG(LAMBDA) + 4D0 / 3D0 
     1                * ( 4D0 * DLOG( 1D0 - X ) / ( 1D0 - X ) 
     2                - 2D0 * DLOG( 1D0 - LAMBDA * X ) / ( 1D0 - X )
     3                - 2D0 / ( 1D0 - X ) + 1D0 / 2D0 * ( ( 1D0 - X ) 
     4                / ( 1D0 - LAMBDA * X )**2D0 ) )
            XC3Q_PLUS = X * XC3Q_PLUS
         ELSE
            XC3Q_PLUS = 0D0
         ENDIF
      ELSEIF(VFNS.EQ."FFN0")THEN
         IF(IPT.EQ.1)THEN
            XC3Q_PLUS = X * ( 4D0 / 3D0 * ( 2D0 * DLOG(1D0-X)/(1D0-X)
     1                - 3D0 / 2D0 / ( 1D0 - X ) ) )
         ELSE
            XC3Q_PLUS = 0D0
         ENDIF
      ENDIF
*
      RETURN
*
      END
*
************************************************************************
*
*     Wrap function of one argument (y) needed to be integrated with
*     DGAUSS.
*
************************************************************************
      FUNCTION C1Q_PLUS_WRAP(Y)
*
      IMPLICIT NONE
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION XC1Q_PLUS
**
*     Output Variables
*
      DOUBLE PRECISION C1Q_PLUS_WRAP
*
      C1Q_PLUS_WRAP = XC1Q_PLUS(Y) / Y
*
      RETURN
      END
*
************************************************************************
      FUNCTION C2Q_PLUS_WRAP(Y)
*
      IMPLICIT NONE
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION XC2Q_PLUS
**
*     Output Variables
*
      DOUBLE PRECISION C2Q_PLUS_WRAP
*
      C2Q_PLUS_WRAP = XC2Q_PLUS(Y) / Y
*
      RETURN
      END
*
************************************************************************
      FUNCTION C3Q_PLUS_WRAP(Y)
*
      IMPLICIT NONE
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION XC3Q_PLUS
**
*     Output Variables
*
      DOUBLE PRECISION C3Q_PLUS_WRAP
*
      C3Q_PLUS_WRAP = XC3Q_PLUS(Y) / Y
*
      RETURN
      END
