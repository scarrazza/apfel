************************************************************************
*
*     SetTensorGluons.f:
*
*     This subroutine enables or disables time-like evolution for the
*     fragmentation functions.
*
************************************************************************
      subroutine SetTensorGluons(ng)
*
      implicit none
*
      include "../commons/TensorGluons.h"
*
      integer ng
      logical tg
*
      nTG            = ng
      InTensorGluons = "done"
*
      return
      end
