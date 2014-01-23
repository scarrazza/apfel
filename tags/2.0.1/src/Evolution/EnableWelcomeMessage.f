************************************************************************
*
*     EnableWelcomeMessage.f:
*
*     This subroutine enable or disable the priting of the Welcome
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
