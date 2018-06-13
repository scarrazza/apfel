************************************************************************
*
*     SetTargetDIS.f:
*
*     This subroutine sets the DIS process.
*
************************************************************************
      subroutine SetTargetDIS(tar)
*
      implicit none
*
      include "../commons/TargetDIS.h"
*
*     Variables
*
      character*(*) tar
*
      TargetDIS = trim(tar)
      InTargetDIS = "done"
*
      return
      end
