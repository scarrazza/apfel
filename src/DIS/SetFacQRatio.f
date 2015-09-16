************************************************************************
*
*     SetFacQRatio.f:
*
*     This subroutine sets the ratio between factorization and
*     DIS scales.
*
************************************************************************
      subroutine SetFacQRatio(ratioF)
*
      implicit none
*
      include "../commons/kfacQ.h"
*
*     Variables
*
      double precision ratioF
*
      kfacQ   = ratioF * ratioF
      InKfacQ = "done"
*
      return
      end
