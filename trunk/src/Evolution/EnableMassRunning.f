************************************************************************
*
*     EnableMassRunning.f:
*
*     This subroutine enables or disables the running of the MSbar
*     masses if they have been set.
*
************************************************************************
      subroutine EnableMassRunning(mr)
*
      implicit none
*
      include "../commons/MassRunning.h"
*
      logical mr
*
      MassRunning   = mr
      InMassRunning = "done"
*
      return
      end
