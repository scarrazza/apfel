************************************************************************
*
*     SetRenQRatio.f:
*
*     This subroutine sets the ratio between renormalization and
*     DIS scales.
*
************************************************************************
      subroutine SetRenQRatio(ratioR)
*
      implicit none
*
      include "../commons/krenQ.h"
*
*     Variables
*
      double precision ratioR
*
      krenQ   = ratioR * ratioR
      InKrenQ = "done"
*
      return
      end
