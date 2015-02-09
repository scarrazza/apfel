************************************************************************
*
*     integrands_cc.f:
*
*     Functions to be integrated for the CC observables.
*
************************************************************************
      FUNCTION XGL_XC1G(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XGL,XC1G
**
*     Output Variables
*
      DOUBLE PRECISION XGL_XC1G
*
      X = XBJ
      Z = X / Y
      XGL_XC1G = XGL(Y) * XC1G(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XPS_XC1Q_UNPLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XSG
      DOUBLE PRECISION XC1Q_UNPLUS_PS
**
*     Output Variables
*
      DOUBLE PRECISION XPS_XC1Q_UNPLUS
*
      X = XBJ
      Z = X / Y
      XPS_XC1Q_UNPLUS = XSG(Y) * XC1Q_UNPLUS_PS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XD_XC1Q_UNPLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XD,XC1Q_UNPLUS
**
*     Output Variables
*
      DOUBLE PRECISION XD_XC1Q_UNPLUS
*
      X = XBJ
      Z = X / Y
      XD_XC1Q_UNPLUS = XD(Y) * XC1Q_UNPLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XD_XC1Q_PLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XD,XC1Q_PLUS
**
*     Output Variables
*
      DOUBLE PRECISION XD_XC1Q_PLUS
*
      X = XBJ
      Z = X / Y
      XD_XC1Q_PLUS = ( XD(Y) - XD(X) ) * XC1Q_PLUS(Z) / Y
      RETURN
*
      END
*
************************************************************************
      FUNCTION XU_XC1Q_UNPLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XU,XC1Q_UNPLUS
**
*     Output Variables
*
      DOUBLE PRECISION XU_XC1Q_UNPLUS
*
      X = XBJ
      Z = X / Y
      XU_XC1Q_UNPLUS = XU(Y) * XC1Q_UNPLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XU_XC1Q_PLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XU,XC1Q_PLUS
**
*     Output Variables
*
      DOUBLE PRECISION XU_XC1Q_PLUS
*
      X = XBJ
      Z = X / Y
      XU_XC1Q_PLUS = ( XU(Y) - XU(X) ) * XC1Q_PLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XS_XC1Q_UNPLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XS,XC1Q_UNPLUS
**
*     Output Variables
*
      DOUBLE PRECISION XS_XC1Q_UNPLUS
*
      X = XBJ
      Z = X / Y
      XS_XC1Q_UNPLUS = XS(Y) * XC1Q_UNPLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XS_XC1Q_PLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XS,XC1Q_PLUS
**
*     Output Variables
*
      DOUBLE PRECISION XS_XC1Q_PLUS
*
      X = XBJ
      Z = X / Y
      XS_XC1Q_PLUS = ( XS(Y) - XS(X) ) * XC1Q_PLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XC_XC1Q_UNPLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XC,XC1Q_UNPLUS
**
*     Output Variables
*
      DOUBLE PRECISION XC_XC1Q_UNPLUS
*
      X = XBJ
      Z = X / Y
      XC_XC1Q_UNPLUS = XC(Y) * XC1Q_UNPLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XC_XC1Q_PLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XC,XC1Q_PLUS
**
*     Output Variables
*
      DOUBLE PRECISION XC_XC1Q_PLUS
*
      X = XBJ
      Z = X / Y
      XC_XC1Q_PLUS = ( XC(Y) - XC(X) ) * XC1Q_PLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XB_XC1Q_UNPLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XB,XC1Q_UNPLUS
**
*     Output Variables
*
      DOUBLE PRECISION XB_XC1Q_UNPLUS
*
      X = XBJ
      Z = X / Y
      XB_XC1Q_UNPLUS = XB(Y) * XC1Q_UNPLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XB_XC1Q_PLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XB,XC1Q_PLUS
**
*     Output Variables
*
      DOUBLE PRECISION XB_XC1Q_PLUS
*
      X = XBJ
      Z = X / Y
      XB_XC1Q_PLUS = ( XB(Y) - XB(X) ) * XC1Q_PLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XT_XC1Q_UNPLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XT,XC1Q_UNPLUS
**
*     Output Variables
*
      DOUBLE PRECISION XT_XC1Q_UNPLUS
*
      X = XBJ
      Z = X / Y
      XT_XC1Q_UNPLUS = XT(Y) * XC1Q_UNPLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XT_XC1Q_PLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XT,XC1Q_PLUS
**
*     Output Variables
*
      DOUBLE PRECISION XT_XC1Q_PLUS
*
      X = XBJ
      Z = X / Y
      XT_XC1Q_PLUS = ( XT(Y) - XT(X) ) * XC1Q_PLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XGL_XC2G(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XGL,XC2G
**
*     Output Variables
*
      DOUBLE PRECISION XGL_XC2G
*
      X = XBJ
      Z = X / Y
      XGL_XC2G = XGL(Y) * XC2G(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XPS_XC2Q_UNPLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XSG
      DOUBLE PRECISION XC2Q_UNPLUS_PS
**
*     Output Variables
*
      DOUBLE PRECISION XPS_XC2Q_UNPLUS
*
      X = XBJ
      Z = X / Y
      XPS_XC2Q_UNPLUS = XSG(Y) * XC2Q_UNPLUS_PS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XD_XC2Q_UNPLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XD,XC2Q_UNPLUS
**
*     Output Variables
*
      DOUBLE PRECISION XD_XC2Q_UNPLUS
*
      X = XBJ
      Z = X / Y
      XD_XC2Q_UNPLUS = XD(Y) * XC2Q_UNPLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XD_XC2Q_PLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XD,XC2Q_PLUS
**
*     Output Variables
*
      DOUBLE PRECISION XD_XC2Q_PLUS
*
      X = XBJ
      Z = X / Y
      XD_XC2Q_PLUS = ( XD(Y) - XD(X) ) * XC2Q_PLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XU_XC2Q_UNPLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XU,XC2Q_UNPLUS
**
*     Output Variables
*
      DOUBLE PRECISION XU_XC2Q_UNPLUS
*
      X = XBJ
      Z = X / Y
      XU_XC2Q_UNPLUS = XU(Y) * XC2Q_UNPLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XU_XC2Q_PLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XU,XC2Q_PLUS
**
*     Output Variables
*
      DOUBLE PRECISION XU_XC2Q_PLUS
*
      X = XBJ
      Z = X / Y
      XU_XC2Q_PLUS = ( XU(Y) - XU(X) ) * XC2Q_PLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XS_XC2Q_UNPLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XS,XC2Q_UNPLUS
**
*     Output Variables
*
      DOUBLE PRECISION XS_XC2Q_UNPLUS
*
      X = XBJ
      Z = X / Y
      XS_XC2Q_UNPLUS = XS(Y) * XC2Q_UNPLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XS_XC2Q_PLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XS,XC2Q_PLUS
**
*     Output Variables
*
      DOUBLE PRECISION XS_XC2Q_PLUS
*
      X = XBJ
      Z = X / Y
      XS_XC2Q_PLUS = ( XS(Y) - XS(X) ) * XC2Q_PLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XC_XC2Q_UNPLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XC,XC2Q_UNPLUS
**
*     Output Variables
*
      DOUBLE PRECISION XC_XC2Q_UNPLUS
*
      X = XBJ
      Z = X / Y
      XC_XC2Q_UNPLUS = XC(Y) * XC2Q_UNPLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XC_XC2Q_PLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XC,XC2Q_PLUS
**
*     Output Variables
*
      DOUBLE PRECISION XC_XC2Q_PLUS
*
      X = XBJ
      Z = X / Y
      XC_XC2Q_PLUS = ( XC(Y) - XC(X) ) * XC2Q_PLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XB_XC2Q_UNPLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XB,XC2Q_UNPLUS
**
*     Output Variables
*
      DOUBLE PRECISION XB_XC2Q_UNPLUS
*
      X = XBJ
      Z = X / Y
      XB_XC2Q_UNPLUS = XB(Y) * XC2Q_UNPLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XB_XC2Q_PLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XB,XC2Q_PLUS
**
*     Output Variables
*
      DOUBLE PRECISION XB_XC2Q_PLUS
*
      X = XBJ
      Z = X / Y
      XB_XC2Q_PLUS = ( XB(Y) - XB(X) ) * XC2Q_PLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XT_XC2Q_UNPLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XT,XC2Q_UNPLUS
**
*     Output Variables
*
      DOUBLE PRECISION XT_XC2Q_UNPLUS
*
      X = XBJ
      Z = X / Y
      XT_XC2Q_UNPLUS = XT(Y) * XC2Q_UNPLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XT_XC2Q_PLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XT,XC2Q_PLUS
**
*     Output Variables
*
      DOUBLE PRECISION XT_XC2Q_PLUS
*
      X = XBJ
      Z = X / Y
      XT_XC2Q_PLUS = ( XT(Y) - XT(X) ) * XC2Q_PLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XGL_XC3G(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XGL,XC3G
**
*     Output Variables
*
      DOUBLE PRECISION XGL_XC3G
*
      X = XBJ
      Z = X / Y
      XGL_XC3G = XGL(Y) * XC3G(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XD_XC3Q_UNPLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XD,XC3Q_UNPLUS
**
*     Output Variables
*
      DOUBLE PRECISION XD_XC3Q_UNPLUS
*
      X = XBJ
      Z = X / Y
      XD_XC3Q_UNPLUS = XD(Y) * XC3Q_UNPLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XD_XC3Q_PLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XD,XC3Q_PLUS
**
*     Output Variables
*
      DOUBLE PRECISION XD_XC3Q_PLUS
*
      X = XBJ
      Z = X / Y
      XD_XC3Q_PLUS = ( XD(Y) - XD(X) ) * XC3Q_PLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XU_XC3Q_UNPLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XU,XC3Q_UNPLUS
**
*     Output Variables
*
      DOUBLE PRECISION XU_XC3Q_UNPLUS
*
      X = XBJ
      Z = X / Y
      XU_XC3Q_UNPLUS = XU(Y) * XC3Q_UNPLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XU_XC3Q_PLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XU,XC3Q_PLUS
**
*     Output Variables
*
      DOUBLE PRECISION XU_XC3Q_PLUS
*
      X = XBJ
      Z = X / Y
      XU_XC3Q_PLUS = ( XU(Y) - XU(X) ) * XC3Q_PLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XS_XC3Q_UNPLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XS,XC3Q_UNPLUS
**
*     Output Variables
*
      DOUBLE PRECISION XS_XC3Q_UNPLUS
*
      X = XBJ
      Z = X / Y
      XS_XC3Q_UNPLUS = XS(Y) * XC3Q_UNPLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XS_XC3Q_PLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XS,XC3Q_PLUS
**
*     Output Variables
*
      DOUBLE PRECISION XS_XC3Q_PLUS
*
      X = XBJ
      Z = X / Y
      XS_XC3Q_PLUS = ( XS(Y) - XS(X) ) * XC3Q_PLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XC_XC3Q_UNPLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XC,XC3Q_UNPLUS
**
*     Output Variables
*
      DOUBLE PRECISION XC_XC3Q_UNPLUS
*
      X = XBJ
      Z = X / Y
      XC_XC3Q_UNPLUS = XC(Y) * XC3Q_UNPLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XC_XC3Q_PLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XC,XC3Q_PLUS
**
*     Output Variables
*
      DOUBLE PRECISION XC_XC3Q_PLUS
*
      X = XBJ
      Z = X / Y
      XC_XC3Q_PLUS = ( XC(Y) - XC(X) ) * XC3Q_PLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XB_XC3Q_UNPLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XB,XC3Q_UNPLUS
**
*     Output Variables
*
      DOUBLE PRECISION XB_XC3Q_UNPLUS
*
      X = XBJ
      Z = X / Y
      XB_XC3Q_UNPLUS = XB(Y) * XC3Q_UNPLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XB_XC3Q_PLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XB,XC3Q_PLUS
**
*     Output Variables
*
      DOUBLE PRECISION XB_XC3Q_PLUS
*
      X = XBJ
      Z = X / Y
      XB_XC3Q_PLUS = ( XB(Y) - XB(X) ) * XC3Q_PLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XT_XC3Q_UNPLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XT,XC3Q_UNPLUS
**
*     Output Variables
*
      DOUBLE PRECISION XT_XC3Q_UNPLUS
*
      X = XBJ
      Z = X / Y
      XT_XC3Q_UNPLUS = XT(Y) * XC3Q_UNPLUS(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XT_XC3Q_PLUS(Y)
*
      IMPLICIT NONE
*
      include "commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XT,XC3Q_PLUS
**
*     Output Variables
*
      DOUBLE PRECISION XT_XC3Q_PLUS
*
      X = XBJ
      Z = X / Y
      XT_XC3Q_PLUS = ( XT(Y) - XT(X) ) * XC3Q_PLUS(Z) / Y
*
      RETURN
*
      END
