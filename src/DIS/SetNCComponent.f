************************************************************************
*
*     SetNCComponent.f:
*
*     This subroutine sets the component of the NC structure functions
*     to be selected.
*
************************************************************************
      subroutine SetNCComponent(cm)
*
      implicit none
*
      include "../commons/NCComponent.h"
*
*     Variables
*
      character*2 cm
*
      NCComponent = cm
      InNCComponent = "done"
*
      return
      end
