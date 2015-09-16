************************************************************************
*
*     SetProcessDIS.f:
*
*     This subroutine sets the DIS process.
*
************************************************************************
      subroutine SetProcessDIS(pr)
*
      implicit none
*
      include "../commons/ProcessDIS.h"
*
*     Variables
*
      character*2 pr
*
      ProcessDIS = pr
      InProcessDIS = "done"
*
      return
      end
