************************************************************************
*
*     EnableWelcomeMessage.f:
*
*     This subroutine enables or disables the priting of the Welcome
*     message.
*
************************************************************************
      subroutine EnableWelcomeMessage(wc)
*
      implicit none
*
      include "../commons/Welcome.h"
*
      logical wc
*
      Welcome = wc
      InWelcome = "done"
*
      return
      end
