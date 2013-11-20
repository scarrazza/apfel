************************************************************************
*
*     integrands_nc.f:
*
*     Functions to be integrated for the NC observables.
*
************************************************************************
      FUNCTION XGL_XC2G_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XGL,XC2G_NC
**
*     Output Variables
*
      DOUBLE PRECISION XGL_XC2G_NC
*
      X = XBJ
      Z = X / Y
      XGL_XC2G_NC = XGL(Y) * XC2G_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XPS_XC2Q_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
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
      DOUBLE PRECISION XC2Q_UNPLUS_PS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XPS_XC2Q_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XPS_XC2Q_UNPLUS_NC = XSG(Y) * XC2Q_UNPLUS_PS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XDP_XC2Q_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XDP,XC2Q_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XDP_XC2Q_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XDP_XC2Q_UNPLUS_NC = XDP(Y) * XC2Q_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XDP_XC2Q_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XDP,XC2Q_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XDP_XC2Q_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XDP_XC2Q_PLUS_NC = ( XDP(Y) - XDP(X) ) * XC2Q_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XUP_XC2Q_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XUP,XC2Q_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XUP_XC2Q_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XUP_XC2Q_UNPLUS_NC = XUP(Y) * XC2Q_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XUP_XC2Q_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XUP,XC2Q_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XUP_XC2Q_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XUP_XC2Q_PLUS_NC = ( XUP(Y) - XUP(X) ) * XC2Q_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XSP_XC2Q_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XSP,XC2Q_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XSP_XC2Q_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XSP_XC2Q_UNPLUS_NC = XSP(Y) * XC2Q_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XSP_XC2Q_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XSP,XC2Q_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XSP_XC2Q_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XSP_XC2Q_PLUS_NC = ( XSP(Y) - XSP(X) ) * XC2Q_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XCP_XC2Q_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XCP,XC2Q_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XCP_XC2Q_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XCP_XC2Q_UNPLUS_NC = XCP(Y) * XC2Q_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XCP_XC2Q_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XCP,XC2Q_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XCP_XC2Q_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XCP_XC2Q_PLUS_NC = ( XCP(Y) - XCP(X) ) * XC2Q_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XBP_XC2Q_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XBP,XC2Q_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XBP_XC2Q_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XBP_XC2Q_UNPLUS_NC = XBP(Y) * XC2Q_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XBP_XC2Q_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XBP,XC2Q_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XBP_XC2Q_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XBP_XC2Q_PLUS_NC = ( XBP(Y) - XBP(X) ) * XC2Q_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XTP_XC2Q_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XTP,XC2Q_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XTP_XC2Q_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XTP_XC2Q_UNPLUS_NC = XTP(Y) * XC2Q_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XTP_XC2Q_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XTP,XC2Q_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XTP_XC2Q_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XTP_XC2Q_PLUS_NC = ( XTP(Y) - XTP(X) ) * XC2Q_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XDM_XC2Q_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XDM,XC2Q_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XDM_XC2Q_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XDM_XC2Q_UNPLUS_NC = XDM(Y) * XC2Q_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XDM_XC2Q_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XDM,XC2Q_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XDM_XC2Q_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XDM_XC2Q_PLUS_NC = ( XDM(Y) - XDM(X) ) * XC2Q_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XUM_XC2Q_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XUM,XC2Q_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XUM_XC2Q_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XUM_XC2Q_UNPLUS_NC = XUM(Y) * XC2Q_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XUM_XC2Q_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XUM,XC2Q_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XUM_XC2Q_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XUM_XC2Q_PLUS_NC = ( XUM(Y) - XUM(X) ) * XC2Q_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XSM_XC2Q_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XSM,XC2Q_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XSM_XC2Q_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XSM_XC2Q_UNPLUS_NC = XSM(Y) * XC2Q_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XSM_XC2Q_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XSM,XC2Q_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XSM_XC2Q_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XSM_XC2Q_PLUS_NC = ( XSM(Y) - XSM(X) ) * XC2Q_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XCM_XC2Q_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XCM,XC2Q_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XCM_XC2Q_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XCM_XC2Q_UNPLUS_NC = XCM(Y) * XC2Q_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XCM_XC2Q_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XCM,XC2Q_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XCM_XC2Q_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XCM_XC2Q_PLUS_NC = ( XCM(Y) - XCM(X) ) * XC2Q_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XBM_XC2Q_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XBM,XC2Q_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XBM_XC2Q_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XBM_XC2Q_UNPLUS_NC = XBM(Y) * XC2Q_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XBM_XC2Q_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XBM,XC2Q_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XBM_XC2Q_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XBM_XC2Q_PLUS_NC = ( XBM(Y) - XBM(X) ) * XC2Q_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XTM_XC2Q_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XTM,XC2Q_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XTM_XC2Q_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XTM_XC2Q_UNPLUS_NC = XTM(Y) * XC2Q_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XTM_XC2Q_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XTM,XC2Q_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XTM_XC2Q_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XTM_XC2Q_PLUS_NC = ( XTM(Y) - XTM(X) ) * XC2Q_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XGL_XC3G_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XGL,XC3G_NC
**
*     Output Variables
*
      DOUBLE PRECISION XGL_XC3G_NC
*
      X = XBJ
      Z = X / Y
      XGL_XC3G_NC = XGL(Y) * XC3G_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XDP_XC3Q_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XDP,XC3Q_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XDP_XC3Q_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XDP_XC3Q_UNPLUS_NC = XDP(Y) * XC3Q_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XDP_XC3Q_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XDP,XC3Q_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XDP_XC3Q_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XDP_XC3Q_PLUS_NC = ( XDP(Y) - XDP(X) ) * XC3Q_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XUP_XC3Q_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XUP,XC3Q_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XUP_XC3Q_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XUP_XC3Q_UNPLUS_NC = XUP(Y) * XC3Q_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XUP_XC3Q_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XUP,XC3Q_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XUP_XC3Q_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XUP_XC3Q_PLUS_NC = ( XUP(Y) - XUP(X) ) * XC3Q_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XSP_XC3Q_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XSP,XC3Q_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XSP_XC3Q_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XSP_XC3Q_UNPLUS_NC = XSP(Y) * XC3Q_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XSP_XC3Q_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XSP,XC3Q_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XSP_XC3Q_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XSP_XC3Q_PLUS_NC = ( XSP(Y) - XSP(X) ) * XC3Q_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XCP_XC3Q_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XCP,XC3Q_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XCP_XC3Q_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XCP_XC3Q_UNPLUS_NC = XCP(Y) * XC3Q_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XCP_XC3Q_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XCP,XC3Q_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XCP_XC3Q_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XCP_XC3Q_PLUS_NC = ( XCP(Y) - XCP(X) ) * XC3Q_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XBP_XC3Q_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XBP,XC3Q_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XBP_XC3Q_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XBP_XC3Q_UNPLUS_NC = XBP(Y) * XC3Q_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XBP_XC3Q_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XBP,XC3Q_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XBP_XC3Q_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XBP_XC3Q_PLUS_NC = ( XBP(Y) - XBP(X) ) * XC3Q_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XTP_XC3Q_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XTP,XC3Q_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XTP_XC3Q_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XTP_XC3Q_UNPLUS_NC = XTP(Y) * XC3Q_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XTP_XC3Q_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XTP,XC3Q_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XTP_XC3Q_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XTP_XC3Q_PLUS_NC = ( XTP(Y) - XTP(X) ) * XC3Q_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XDM_XC3Q_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XDM,XC3Q_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XDM_XC3Q_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XDM_XC3Q_UNPLUS_NC = XDM(Y) * XC3Q_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XDM_XC3Q_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XDM,XC3Q_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XDM_XC3Q_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XDM_XC3Q_PLUS_NC = ( XDM(Y) - XDM(X) ) * XC3Q_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XUM_XC3Q_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XUM,XC3Q_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XUM_XC3Q_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XUM_XC3Q_UNPLUS_NC = XUM(Y) * XC3Q_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XUM_XC3Q_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XUM,XC3Q_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XUM_XC3Q_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XUM_XC3Q_PLUS_NC = ( XUM(Y) - XUM(X) ) * XC3Q_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XSM_XC3Q_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XSM,XC3Q_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XSM_XC3Q_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XSM_XC3Q_UNPLUS_NC = XSM(Y) * XC3Q_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XSM_XC3Q_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XSM,XC3Q_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XSM_XC3Q_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XSM_XC3Q_PLUS_NC = ( XSM(Y) - XSM(X) ) * XC3Q_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XCM_XC3Q_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XCM,XC3Q_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XCM_XC3Q_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XCM_XC3Q_UNPLUS_NC = XCM(Y) * XC3Q_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XCM_XC3Q_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XCM,XC3Q_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XCM_XC3Q_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XCM_XC3Q_PLUS_NC = ( XCM(Y) - XCM(X) ) * XC3Q_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XBM_XC3Q_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XBM,XC3Q_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XBM_XC3Q_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XBM_XC3Q_UNPLUS_NC = XBM(Y) * XC3Q_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XBM_XC3Q_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XBM,XC3Q_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XBM_XC3Q_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XBM_XC3Q_PLUS_NC = ( XBM(Y) - XBM(X) ) * XC3Q_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XTM_XC3Q_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XTM,XC3Q_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XTM_XC3Q_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XTM_XC3Q_UNPLUS_NC = XTM(Y) * XC3Q_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XTM_XC3Q_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XTM,XC3Q_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XTM_XC3Q_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XTM_XC3Q_PLUS_NC = ( XTM(Y) - XTM(X) ) * XC3Q_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XGL_XCLG_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XGL,XCLG_NC
**
*     Output Variables
*
      DOUBLE PRECISION XGL_XCLG_NC
*
      X = XBJ
      Z = X / Y
      XGL_XCLG_NC = XGL(Y) * XCLG_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XPS_XCLQ_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
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
      DOUBLE PRECISION XCLQ_UNPLUS_PS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XPS_XCLQ_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XPS_XCLQ_UNPLUS_NC = XSG(Y) * XCLQ_UNPLUS_PS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XDP_XCLQ_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XDP,XCLQ_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XDP_XCLQ_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XDP_XCLQ_UNPLUS_NC = XDP(Y) * XCLQ_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XDP_XCLQ_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XDP,XCLQ_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XDP_XCLQ_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XDP_XCLQ_PLUS_NC = ( XDP(Y) - XDP(X) ) * XCLQ_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XUP_XCLQ_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XUP,XCLQ_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XUP_XCLQ_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XUP_XCLQ_UNPLUS_NC = XUP(Y) * XCLQ_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XUP_XCLQ_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XUP,XCLQ_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XUP_XCLQ_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XUP_XCLQ_PLUS_NC = ( XUP(Y) - XUP(X) ) * XCLQ_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XSP_XCLQ_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XSP,XCLQ_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XSP_XCLQ_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XSP_XCLQ_UNPLUS_NC = XSP(Y) * XCLQ_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XSP_XCLQ_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XSP,XCLQ_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XSP_XCLQ_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XSP_XCLQ_PLUS_NC = ( XSP(Y) - XSP(X) ) * XCLQ_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XCP_XCLQ_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XCP,XCLQ_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XCP_XCLQ_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XCP_XCLQ_UNPLUS_NC = XCP(Y) * XCLQ_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XCP_XCLQ_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XCP,XCLQ_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XCP_XCLQ_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XCP_XCLQ_PLUS_NC = ( XCP(Y) - XCP(X) ) * XCLQ_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XBP_XCLQ_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XBP,XCLQ_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XBP_XCLQ_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XBP_XCLQ_UNPLUS_NC = XBP(Y) * XCLQ_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XBP_XCLQ_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XBP,XCLQ_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XBP_XCLQ_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XBP_XCLQ_PLUS_NC = ( XBP(Y) - XBP(X) ) * XCLQ_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XTP_XCLQ_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XTP,XCLQ_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XTP_XCLQ_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XTP_XCLQ_UNPLUS_NC = XTP(Y) * XCLQ_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XTP_XCLQ_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XTP,XCLQ_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XTP_XCLQ_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XTP_XCLQ_PLUS_NC = ( XTP(Y) - XTP(X) ) * XCLQ_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XDM_XCLQ_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XDM,XCLQ_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XDM_XCLQ_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XDM_XCLQ_UNPLUS_NC = XDM(Y) * XCLQ_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XDM_XCLQ_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XDM,XCLQ_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XDM_XCLQ_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XDM_XCLQ_PLUS_NC = ( XDM(Y) - XDM(X) ) * XCLQ_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XUM_XCLQ_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XUM,XCLQ_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XUM_XCLQ_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XUM_XCLQ_UNPLUS_NC = XUM(Y) * XCLQ_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XUM_XCLQ_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XUM,XCLQ_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XUM_XCLQ_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XUM_XCLQ_PLUS_NC = ( XUM(Y) - XUM(X) ) * XCLQ_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XSM_XCLQ_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XSM,XCLQ_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XSM_XCLQ_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XSM_XCLQ_UNPLUS_NC = XSM(Y) * XCLQ_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XSM_XCLQ_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XSM,XCLQ_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XSM_XCLQ_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XSM_XCLQ_PLUS_NC = ( XSM(Y) - XSM(X) ) * XCLQ_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XCM_XCLQ_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XCM,XCLQ_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XCM_XCLQ_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XCM_XCLQ_UNPLUS_NC = XCM(Y) * XCLQ_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XCM_XCLQ_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XCM,XCLQ_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XCM_XCLQ_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XCM_XCLQ_PLUS_NC = ( XCM(Y) - XCM(X) ) * XCLQ_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XBM_XCLQ_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XBM,XCLQ_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XBM_XCLQ_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XBM_XCLQ_UNPLUS_NC = XBM(Y) * XCLQ_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XBM_XCLQ_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
      include "../commons/scale.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XBM,XCLQ_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XBM_XCLQ_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XBM_XCLQ_PLUS_NC = ( XBM(Y) - XBM(X) ) * XCLQ_PLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XTM_XCLQ_UNPLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
      include "../commons/scale.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XTM,XCLQ_UNPLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XTM_XCLQ_UNPLUS_NC
*
      X = XBJ
      Z = X / Y
      XTM_XCLQ_UNPLUS_NC = XTM(Y) * XCLQ_UNPLUS_NC(Z) / Y
*
      RETURN
*
      END
*
************************************************************************
      FUNCTION XTM_XCLQ_PLUS_NC(Y)
*
      IMPLICIT NONE
*
      include "../commons/xbj.h"
      include "../commons/scale.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION Z
      DOUBLE PRECISION X
      DOUBLE PRECISION XTM,XCLQ_PLUS_NC
**
*     Output Variables
*
      DOUBLE PRECISION XTM_XCLQ_PLUS_NC
*
      X = XBJ
      Z = X / Y
      XTM_XCLQ_PLUS_NC = ( XTM(Y) - XTM(X) ) * XCLQ_PLUS_NC(Z) / Y
*
      RETURN
*
      END
