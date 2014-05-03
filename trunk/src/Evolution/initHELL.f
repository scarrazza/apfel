************************************************************************
*
*     initHELL.f:
*
*     This routine initializes the HELL code by Bonvini for the small-x
*     resummation.
*
************************************************************************
      subroutine initHELL(LogAcc,ipt,asmc,asmb,asmt)
*
      implicit none
**
*     Input Variables
*
      double precision LogAcc,ipt
      double precision asmc,asmb,asmt
**
*     Internal Variables
*
      call HELLLogOrder(LogAcc)
      call HELLOrder(ipt)
      call HELL(asmc,asmb,asmt)
*
      return
      end
