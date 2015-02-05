************************************************************************
*
*     SetGFermi.f:
*
*     This subroutine sets the Fermi constant.
*
************************************************************************
      subroutine SetGFermi(gf)
*
      implicit none
*
      include "../commons/GFermi.h"
*
*     Variables
*
      double precision gf
*
      GFermi   = gf
      InGFermi = "done"
*
      return
      end
