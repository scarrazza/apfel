************************************************************************
*
*     xpdf.f:
*
*     The following functions return x times the chosen PDF.
*
************************************************************************
*
*     Gluon
*
************************************************************************
      FUNCTION XGL(Y)
*
      IMPLICIT NONE
*
      include "../commons/scale.h"
      include "../commons/pdf.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION YPDF(-6:6),XPDF
**
*     Output Variables
*
      DOUBLE PRECISION XGL
*
      IF(PDF.EQ."LHA")THEN
         CALL EVOLVEPDF(Y,Q,YPDF)
         XGL = YPDF(0)
      ELSEIF(PDF.EQ."APF")THEN
         XGL = xPDF(0,Y)
      ENDIF
*
      RETURN
      END
*
************************************************************************
*
*     Singlet
*
************************************************************************
      FUNCTION XSG(Y)
*
      IMPLICIT NONE
*
      include "../commons/scale.h"
      include "../commons/pdf.h"
      include "../commons/vfns.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      INTEGER I
      DOUBLE PRECISION YPDF(-6:6),XPDF
**
*     Output Variables
*
      DOUBLE PRECISION XSG
*
      XSG = 0D0
      IF(PDF.EQ."LHA")THEN
         CALL EVOLVEPDF(Y,Q,YPDF)
         IF(VFNS.EQ."ZMVN")THEN
            DO I=1,6
               XSG = XSG + YPDF(I) + YPDF(-I)
            ENDDO
         ELSEIF(VFNS.EQ."FFNS".OR.VFNS.EQ."FFN0")THEN
            DO I=1,3
               XSG = XSG + YPDF(I) + YPDF(-I)
            ENDDO
         ENDIF
      ELSEIF(PDF.EQ."APF")THEN
         IF(VFNS.EQ."ZMVN")THEN
            DO I=1,6
               XSG = XSG + XPDF(I,Y) + XPDF(-I,Y)
            ENDDO
         ELSEIF(VFNS.EQ."FFNS".OR.VFNS.EQ."FFN0")THEN
            DO I=1,3
               XSG = XSG + XPDF(I,Y) + XPDF(-I,Y)
            ENDDO
         ENDIF
      ENDIF
*
      RETURN
      END
*
************************************************************************
*
*     NEUTRAL CURRENT combinations.
*
************************************************************************
      FUNCTION XDP(Y)
*
      IMPLICIT NONE
*
      include "../commons/scale.h"
      include "../commons/pdf.h"
      include "../commons/nucleon.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION YPDF(-6:6),XPDF
**
*     Output Variables
*
      DOUBLE PRECISION XDP
*
      IF(PDF.EQ."LHA")THEN
         CALL EVOLVEPDF(Y,Q,YPDF)
         IF(INUCL.EQ.0)THEN
            XDP = YPDF(1) + YPDF(-1)
         ELSEIF(INUCL.EQ.1)THEN
            XDP = YPDF(2) + YPDF(-2)
         ELSEIF(INUCL.EQ.2)THEN
            XDP = ( YPDF(1) + YPDF(-1) + YPDF(2) + YPDF(-2) ) / 2D0
         ENDIF
      ELSEIF(PDF.EQ."APF")THEN
         IF(INUCL.EQ.0)THEN
            XDP = XPDF(1,Y) + XPDF(-1,Y)
         ELSEIF(INUCL.EQ.1)THEN
            XDP = XPDF(2,Y) + XPDF(-2,Y)
         ELSEIF(INUCL.EQ.2)THEN
            XDP = (XPDF(1,Y) + XPDF(-1,Y) + XPDF(2,Y) + XPDF(-2,Y))/2D0
         ENDIF
      ENDIF
*
      RETURN
      END
*
************************************************************************
      FUNCTION XDM(Y)
*
      IMPLICIT NONE
*
      include "../commons/scale.h"
      include "../commons/pdf.h"
      include "../commons/nucleon.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION YPDF(-6:6),XPDF
**
*     Output Variables
*
      DOUBLE PRECISION XDM
*
      IF(PDF.EQ."LHA")THEN
         CALL EVOLVEPDF(Y,Q,YPDF)
         IF(INUCL.EQ.0)THEN
            XDM = YPDF(1) - YPDF(-1)
         ELSEIF(INUCL.EQ.1)THEN
            XDM = YPDF(2) - YPDF(-2)
         ELSEIF(INUCL.EQ.2)THEN
            XDM = ( YPDF(1) - YPDF(-1) + YPDF(2) - YPDF(-2) ) / 2D0
         ENDIF
      ELSEIF(PDF.EQ."APF")THEN
         IF(INUCL.EQ.0)THEN
            XDM = XPDF(1,Y) - XPDF(-1,Y)
         ELSEIF(INUCL.EQ.1)THEN
            XDM = XPDF(2,Y) - XPDF(-2,Y)
         ELSEIF(INUCL.EQ.2)THEN
            XDM = (XPDF(1,Y) - XPDF(-1,Y) + XPDF(2,Y) - XPDF(-2,Y))/2D0
         ENDIF
      ENDIF
*
      RETURN
      END
*
************************************************************************
      FUNCTION XUP(Y)
*
      IMPLICIT NONE
*
      include "../commons/scale.h"
      include "../commons/pdf.h"
      include "../commons/nucleon.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION YPDF(-6:6),XPDF
**
*     Output Variables
*
      DOUBLE PRECISION XUP
*
      IF(PDF.EQ."LHA")THEN
         CALL EVOLVEPDF(Y,Q,YPDF)
         IF(INUCL.EQ.0)THEN
            XUP = YPDF(2) + YPDF(-2)
         ELSEIF(INUCL.EQ.1)THEN
            XUP = YPDF(1) + YPDF(-1)
         ELSEIF(INUCL.EQ.2)THEN
            XUP = ( YPDF(2) + YPDF(-2) + YPDF(1) + YPDF(-1) ) / 2D0
         ENDIF
      ELSEIF(PDF.EQ."APF")THEN
         IF(INUCL.EQ.0)THEN
            XUP = XPDF(2,Y) + XPDF(-2,Y)
         ELSEIF(INUCL.EQ.1)THEN
            XUP = XPDF(1,Y) + XPDF(-1,Y)
         ELSEIF(INUCL.EQ.2)THEN
            XUP = (XPDF(2,Y) + XPDF(-2,Y) + XPDF(1,Y) + XPDF(-1,Y))/2D0
         ENDIF
      ENDIF
*
      RETURN
      END
*
************************************************************************
      FUNCTION XUM(Y)
*
      IMPLICIT NONE
*
      include "../commons/scale.h"
      include "../commons/pdf.h"
      include "../commons/nucleon.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION YPDF(-6:6),XPDF
**
*     Output Variables
*
      DOUBLE PRECISION XUM
*
      IF(PDF.EQ."LHA")THEN
         CALL EVOLVEPDF(Y,Q,YPDF)
         IF(INUCL.EQ.0)THEN
            XUM = YPDF(2) - YPDF(-2)
         ELSEIF(INUCL.EQ.1)THEN
            XUM = YPDF(1) - YPDF(-1)
         ELSEIF(INUCL.EQ.2)THEN
            XUM = ( YPDF(2) - YPDF(-2) + YPDF(1) - YPDF(-1) ) / 2D0
         ENDIF
      ELSEIF(PDF.EQ."APF")THEN
         IF(INUCL.EQ.0)THEN
            XUM = XPDF(2,Y) - XPDF(-2,Y)
         ELSEIF(INUCL.EQ.1)THEN
            XUM = XPDF(1,Y) - XPDF(-1,Y)
         ELSEIF(INUCL.EQ.2)THEN
            XUM = (XPDF(2,Y) - XPDF(-2,Y) + XPDF(1,Y) - XPDF(-1,Y))/2D0
         ENDIF
      ENDIF
*
      RETURN
      END
*
************************************************************************
      FUNCTION XSP(Y)
*
      IMPLICIT NONE
*
      include "../commons/scale.h"
      include "../commons/pdf.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION YPDF(-6:6),XPDF
**
*     Output Variables
*
      DOUBLE PRECISION XSP
*
      IF(PDF.EQ."LHA")THEN
         CALL EVOLVEPDF(Y,Q,YPDF)
         XSP = YPDF(3) + YPDF(-3)
      ELSEIF(PDF.EQ."APF")THEN
         XSP = XPDF(3,Y) + XPDF(-3,Y)
      ENDIF
*
      RETURN
      END
*
************************************************************************
      FUNCTION XSM(Y)
*
      IMPLICIT NONE
*
      include "../commons/scale.h"
      include "../commons/pdf.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION YPDF(-6:6),XPDF
**
*     Output Variables
*
      DOUBLE PRECISION XSM
*
      IF(PDF.EQ."LHA")THEN
         CALL EVOLVEPDF(Y,Q,YPDF)
         XSM = YPDF(3) - YPDF(-3)
      ELSEIF(PDF.EQ."APF")THEN
         XSM = XPDF(3,Y) - XPDF(-3,Y)
      ENDIF
*
      RETURN
      END
*
************************************************************************
      FUNCTION XCP(Y)
*
      IMPLICIT NONE
*
      include "../commons/scale.h"
      include "../commons/pdf.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION YPDF(-6:6),XPDF
**
*     Output Variables
*
      DOUBLE PRECISION XCP
*
      IF(PDF.EQ."LHA")THEN
         CALL EVOLVEPDF(Y,Q,YPDF)
         XCP = YPDF(4) + YPDF(-4)
      ELSEIF(PDF.EQ."APF")THEN
         XCP = XPDF(4,Y) + XPDF(-4,Y)
      ENDIF
*
      RETURN
      END
*
************************************************************************
      FUNCTION XCM(Y)
*
      IMPLICIT NONE
*
      include "../commons/scale.h"
      include "../commons/pdf.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION YPDF(-6:6),XPDF
**
*     Output Variables
*
      DOUBLE PRECISION XCM
*
      IF(PDF.EQ."LHA")THEN
         CALL EVOLVEPDF(Y,Q,YPDF)
         XCM = YPDF(4) - YPDF(-4)
      ELSEIF(PDF.EQ."APF")THEN
         XCM = XPDF(4,Y) - XPDF(-4,Y)
      ENDIF
*
      RETURN
      END
*
************************************************************************
      FUNCTION XBP(Y)
*
      IMPLICIT NONE
*
      include "../commons/scale.h"
      include "../commons/pdf.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION YPDF(-6:6),XPDF
**
*     Output Variables
*
      DOUBLE PRECISION XBP
*
      IF(PDF.EQ."LHA")THEN
         CALL EVOLVEPDF(Y,Q,YPDF)
         XBP = YPDF(5) + YPDF(-5)
      ELSEIF(PDF.EQ."APF")THEN
         XBP = XPDF(5,Y) + XPDF(-5,Y)
      ENDIF
*
      RETURN
      END
*
************************************************************************
      FUNCTION XBM(Y)
*
      IMPLICIT NONE
*
      include "../commons/scale.h"
      include "../commons/pdf.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION YPDF(-6:6),XPDF
**
*     Output Variables
*
      DOUBLE PRECISION XBM
*
      IF(PDF.EQ."LHA")THEN
         CALL EVOLVEPDF(Y,Q,YPDF)
         XBM = YPDF(5) - YPDF(-5)
      ELSEIF(PDF.EQ."APF")THEN
         XBM = XPDF(5,Y) - XPDF(-5,Y)
      ENDIF
*
      RETURN
      END
*
************************************************************************
      FUNCTION XTP(Y)
*
      IMPLICIT NONE
*
      include "../commons/scale.h"
      include "../commons/pdf.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION YPDF(-6:6),XPDF
**
*     Output Variables
*
      DOUBLE PRECISION XTP
*
      IF(PDF.EQ."LHA")THEN
         CALL EVOLVEPDF(Y,Q,YPDF)
         XTP = YPDF(6) + YPDF(-6)
      ELSEIF(PDF.EQ."APF")THEN
         XTP = XPDF(6,Y) + XPDF(-6,Y)
      ENDIF
*
      RETURN
      END
*
************************************************************************
      FUNCTION XTM(Y)
*
      IMPLICIT NONE
*
      include "../commons/scale.h"
      include "../commons/pdf.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION YPDF(-6:6),XPDF
**
*     Output Variables
*
      DOUBLE PRECISION XTM
*
      IF(PDF.EQ."LHA")THEN
         CALL EVOLVEPDF(Y,Q,YPDF)
         XTM = YPDF(6) - YPDF(-6)
      ELSEIF(PDF.EQ."APF")THEN
         XTM = XPDF(6,Y) - XPDF(-6,Y)
      ENDIF
*
      RETURN
      END
*
************************************************************************
*
*     CHARGED CURRENT combinations.
*
************************************************************************
      FUNCTION XD(Y)
*
      IMPLICIT NONE
*
      include "../commons/scale.h"
      include "../commons/pdf.h"
      include "../commons/neut.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION YPDF(-6:6),XPDF
**
*     Output Variables
*
      DOUBLE PRECISION XD
*
      IF(PDF.EQ."LHA")THEN
         CALL EVOLVEPDF(Y,Q,YPDF)
         IF(IH.EQ.4)THEN
            IF(NEUT.EQ.1)THEN
               XD = YPDF(1)
            ELSEIF(NEUT.EQ.-1)THEN
               XD = YPDF(-1)
            ENDIF
         ELSEIF(IH.EQ.5)THEN
            IF(NEUT.EQ.1)THEN
               XD = YPDF(-1)
            ELSEIF(NEUT.EQ.-1)THEN
               XD = YPDF(1)
            ENDIF
         ELSEIF(IH.EQ.6)THEN
            IF(NEUT.EQ.1)THEN
               XD = YPDF(1)
            ELSEIF(NEUT.EQ.-1)THEN
               XD = YPDF(-1)
            ENDIF
         ENDIF
      ELSEIF(PDF.EQ."APF")THEN
         IF(IH.EQ.4)THEN
            IF(NEUT.EQ.1)THEN
               XD = XPDF(1,Y)
            ELSEIF(NEUT.EQ.-1)THEN
               XD = XPDF(-1,Y)
            ENDIF
         ELSEIF(IH.EQ.5)THEN
            IF(NEUT.EQ.1)THEN
               XD = XPDF(-1,Y)
            ELSEIF(NEUT.EQ.-1)THEN
               XD = XPDF(1,Y)
            ENDIF
         ELSEIF(IH.EQ.6)THEN
            IF(NEUT.EQ.1)THEN
               XD = XPDF(1,Y)
            ELSEIF(NEUT.EQ.-1)THEN
               XD = XPDF(-1,Y)
            ENDIF
         ENDIF
      ENDIF
*
      RETURN
      END
*
************************************************************************
      FUNCTION XU(Y)
*
      IMPLICIT NONE
*
      include "../commons/scale.h"
      include "../commons/pdf.h"
      include "../commons/neut.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION YPDF(-6:6),XPDF
**
*     Output Variables
*
      DOUBLE PRECISION XU
*
      IF(PDF.EQ."LHA")THEN
         CALL EVOLVEPDF(Y,Q,YPDF)
         IF(IH.EQ.4)THEN
            IF(NEUT.EQ.1)THEN
               XU = YPDF(2)
            ELSEIF(NEUT.EQ.-1)THEN
               XU = YPDF(-2)
            ENDIF
         ELSEIF(IH.EQ.5)THEN
            IF(NEUT.EQ.1)THEN
               XU = YPDF(-2)
            ELSEIF(NEUT.EQ.-1)THEN
               XU = YPDF(2)
            ENDIF
         ELSEIF(IH.EQ.6)THEN
            IF(NEUT.EQ.1)THEN
               XU = YPDF(2)
            ELSEIF(NEUT.EQ.-1)THEN
               XU = YPDF(-2)
            ENDIF
         ENDIF
      ELSEIF(PDF.EQ."APF")THEN
         IF(IH.EQ.4)THEN
            IF(NEUT.EQ.1)THEN
               XU = XPDF(2,Y)
            ELSEIF(NEUT.EQ.-1)THEN
               XU = XPDF(-2,Y)
            ENDIF
         ELSEIF(IH.EQ.5)THEN
            IF(NEUT.EQ.1)THEN
               XU = XPDF(-2,Y)
            ELSEIF(NEUT.EQ.-1)THEN
               XU = XPDF(2,Y)
            ENDIF
         ELSEIF(IH.EQ.6)THEN
            IF(NEUT.EQ.1)THEN
               XU = XPDF(2,Y)
            ELSEIF(NEUT.EQ.-1)THEN
               XU = XPDF(-2,Y)
            ENDIF
         ENDIF
      ENDIF
*
      RETURN
      END
*
************************************************************************
      FUNCTION XS(Y)
*
      IMPLICIT NONE
*
      include "../commons/scale.h"
      include "../commons/pdf.h"
      include "../commons/neut.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION YPDF(-6:6),XPDF
**
*     Output Variables
*
      DOUBLE PRECISION XS
*
      IF(PDF.EQ."LHA")THEN
         CALL EVOLVEPDF(Y,Q,YPDF)
         IF(NEUT.EQ.1)THEN
            XS = YPDF(3)
         ELSEIF(NEUT.EQ.-1)THEN
            XS = YPDF(-3)
         ENDIF
      ELSEIF(PDF.EQ."APF")THEN
         IF(NEUT.EQ.1)THEN
            XS = XPDF(3,Y)
         ELSEIF(NEUT.EQ.-1)THEN
            XS = XPDF(-3,Y)
         ENDIF
      ENDIF
*
      RETURN
      END
*
************************************************************************
      FUNCTION XC(Y)
*
      IMPLICIT NONE
*
      include "../commons/scale.h"
      include "../commons/pdf.h"
      include "../commons/neut.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION YPDF(-6:6),XPDF
**
*     Output Variables
*
      DOUBLE PRECISION XC
*
      IF(PDF.EQ."LHA")THEN
         CALL EVOLVEPDF(Y,Q,YPDF)
         IF(NEUT.EQ.1)THEN
            XC = YPDF(-4)
         ELSEIF(NEUT.EQ.-1)THEN
            XC = YPDF(4)
         ENDIF
      ELSEIF(PDF.EQ."APF")THEN
         IF(NEUT.EQ.1)THEN
            XC = XPDF(-4,Y)
         ELSEIF(NEUT.EQ.-1)THEN
            XC = XPDF(4,Y)
         ENDIF
      ENDIF
*
      RETURN
      END
*
************************************************************************
      FUNCTION XB(Y)
*
      IMPLICIT NONE
*
      include "../commons/scale.h"
      include "../commons/pdf.h"
      include "../commons/neut.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION YPDF(-6:6),XPDF
**
*     Output Variables
*
      DOUBLE PRECISION XB
*
      IF(PDF.EQ."LHA")THEN
         CALL EVOLVEPDF(Y,Q,YPDF)
         IF(NEUT.EQ.1)THEN
            XB = YPDF(5)
         ELSEIF(NEUT.EQ.-1)THEN
            XB = YPDF(-5)
         ENDIF
      ELSEIF(PDF.EQ."APF")THEN
         IF(NEUT.EQ.1)THEN
            XB = XPDF(5,Y)
         ELSEIF(NEUT.EQ.-1)THEN
            XB = XPDF(-5,Y)
         ENDIF
      ENDIF
*
      RETURN
      END
*
************************************************************************
      FUNCTION XT(Y)
*
      IMPLICIT NONE
*
      include "../commons/scale.h"
      include "../commons/pdf.h"
      include "../commons/neut.h"
**
*     Input Variables
*
      DOUBLE PRECISION Y
**
*     Internal Variables
*
      DOUBLE PRECISION YPDF(-6:6),XPDF
**
*     Output Variables
*
      DOUBLE PRECISION XT
*
      IF(PDF.EQ."LHA")THEN
         CALL EVOLVEPDF(Y,Q,YPDF)
         IF(NEUT.EQ.1)THEN
            XT = YPDF(-6)
         ELSEIF(NEUT.EQ.-1)THEN
            XT = YPDF(6)
         ENDIF
      ELSEIF(PDF.EQ."APF")THEN
         IF(NEUT.EQ.1)THEN
            XT = XPDF(-6,Y)
         ELSEIF(NEUT.EQ.-1)THEN
            XT = XPDF(6,Y)
         ENDIF
      ENDIF
*
      RETURN
      END
