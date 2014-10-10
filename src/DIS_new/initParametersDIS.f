************************************************************************
*
*     initParametersDIS.f:
*
*     It sets all the parameters for the computation of the DIS structure
*     functions if they were not set exernally before by the user.
*
************************************************************************
      subroutine initParametersDIS
*
      implicit none
*
      include "../commons/MassScheme.h"
*
*     Initialize default parameters (those that were not initialized before)
*
      if(InMassScheme.ne."done") call SetMassScheme("ZM-VFNS")   ! "FONLL-A", "FONLL-B", "FONLL-C" or "FFNS"
*
      return
      end
