************************************************************************
*
*     SetProjectileDIS.f:
*
*     This subroutine sets the DIS process.
*
************************************************************************
      subroutine SetProjectileDIS(lept)
*
      implicit none
*
      include "../commons/ProjectileDIS.h"
*
*     Variables
*
      character*(*) lept
*
      ProjectileDIS = trim(lept)
      InProjectileDIS = "done"
*
      return
      end
