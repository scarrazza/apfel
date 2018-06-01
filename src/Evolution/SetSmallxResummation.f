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
      character*(*) la
*
      Smallx = sx
      LogAcc = - 1
      if(la.eq."LL")then
         LogAcc = 0
      elseif(la.eq."NLL")then
         LogAcc = 1
      endif
      InSmallx = "done"
*
c      if(Smallx)then
c         write(6,*) "Small-x resummation currently not available."
c         write(6,*) "Disable it to use APFEL."
c         write(6,*) " "
c         call exit(-10)
c      endif
*
      return
      end
