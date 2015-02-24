************************************************************************
*
*     SetSmallxResummation.f:
*
*     This subroutine enables or disables the small-x resummation allowing
*     to specify the log accuracy ("LL" and "NLL"). The log accuracy is
*     ignored if the resummation is set to .false..
*
************************************************************************
      subroutine SetSmallxResummation(sx,la)
*
      implicit none
*
      include "../commons/Smallx.h"
*
      logical sx
      character*3 la
*
      Smallx = sx
      LogAcc = - 1
      if(la(1:2).eq."LL")then
         LogAcc = 0
      elseif(la(1:3).eq."NLL")then
         LogAcc = 1
      endif
      InSmallx = "done"
*
      if(Smallx)then
         write(6,*) "Small-x resummation currently not available."
         write(6,*) "Disable it to use APFEL."
         write(6,*) " "
         call exit(-10)
      endif
*
      return
      end
