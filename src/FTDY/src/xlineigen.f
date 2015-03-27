************************************************************************
*
*     ELIN is used for interpolation of NLO DY coeff function of a 
*     sqrt(log(1/x)) grid with a total of nxDY points
*
************************************************************************
      FUNCTION ELIN(I,X)
*
      IMPLICIT NONE
*
      include "../commons/mxgridsizeDY.h"
      include "../commons/xxDY.h"
      include "../commons/xgridDY.h"
**
*     Input Variables
*
      INTEGER I
      DOUBLE PRECISION X
**
*     Internal Variables
*
      DOUBLE PRECISION Y,YY(nxDY)
      DOUBLE PRECISION HP,HM
**
*     Output Variables
*
      DOUBLE PRECISION ELIN
*
      ELIN = 0D0
      HP   = 0D0
      HM   = 0D0
*
*     Note FastNLO-like x grid assumed
*
      Y     = DSQRT(DLOG10(1D0/X))
      YY(I) = DSQRT(DLOG10(1d0/xxDY(I)))
*
      IF(I.LT.nxDY) YY(I+1) = DSQRT(DLOG10(1d0/xxDY(I+1)))
      IF(I.LT.nxDY) HP  = YY(I+1) - YY(I)
      IF(I.GT.1)    YY(I-1) = DSQRT(DLOG10(1d0/xxDY(I-1)))
      IF(I.GT.1)    HM = YY(I) - YY(I-1)
*
*     if x<x1 => log(1/x) < log(1/x1) => reverse inequalities
*
      IF(I.EQ.1)THEN
         IF(Y.GE.YY(I+1)) ELIN = ( YY(I+1) - Y ) / HP
      ELSEIF(I.EQ.nxDY)THEN
         IF(Y.LE.YY(I-I)) ELIN = ( Y - YY(I-1) ) / HM
      ELSE
         IF(Y.LE.YY(I-1).AND.Y.GE.YY(I))THEN
            ELIN = ( Y - YY(I-1) ) / HM
         ELSEIF(Y.LT.YY(I).AND.Y.GE.YY(I+1))THEN
            ELIN = ( YY(I+1) - Y ) / HP
         ENDIF
      ENDIF
*
      RETURN
      END
