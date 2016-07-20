************************************************************************
*
*     initHELL.f:
*
*     This routine initializes the HELL code by Bonvini for the small-x
*     resummation.
*
************************************************************************
      subroutine initHELL(la,pt)
*
      implicit none
**
*     Input Variables
*
      integer la,pt
*
      call HELLLogOrder(la)
      call HELLOrder(pt)
      call HELL
*
      return
      end
