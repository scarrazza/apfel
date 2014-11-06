************************************************************************
*
*     SetPerturbativeOrder.f:
*
*     This subroutine sets the perturbative oder of the evolution.
*
************************************************************************
      subroutine SetPerturbativeOrder(pto)
*
      implicit none
*
      include "../commons/ipt.h"
*
*     Variables
*
      integer pto
*
      ipt  = pto
      InPt = "done"
*
      return
      end
