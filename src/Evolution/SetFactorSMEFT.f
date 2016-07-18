************************************************************************
*
*     SetFactorSMEFT.f:
*
*     This subroutine sets the value of C_G / Lambda^2 required by the
*     SMEFT approach.
*
************************************************************************
      subroutine SetFactorSMEFT(eft)
*
      implicit none
*
      include "../commons/SMEFT.h"
*
*     Variables
*
      double precision eft
*
      factSMEFT = eft
      InSMEFT   = "done"
*
      return
      end
