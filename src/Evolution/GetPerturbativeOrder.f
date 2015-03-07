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
      if(inPt.ne."done")then
         write(6,*) "GetPerturbativeOrder: Parameter not initialized"
         write(6,*) "Set it by means of 'SetPerturbativeOrder'"
         write(6,*) "  "
         call exit(-10)
      endif
*
      GetPerturbativeOrder = ipt
*
      return
      end
