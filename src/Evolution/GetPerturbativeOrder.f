************************************************************************
*
*     GetPerturbativeOrder.f:
*
*     This function returns the perturbative oder of the evolution.
*
************************************************************************
      function GetPerturbativeOrder()
*
      implicit none
*
      include "../commons/ipt.h"
*
*     Variables
*
      integer GetPerturbativeOrder
*
      GetPerturbativeOrder = ipt
*
      return
      end
