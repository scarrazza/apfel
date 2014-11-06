************************************************************************
*
*     SetMSbarMasses.f:
*
*     This subroutine sets as a default the heavy quark MSbar masses. 
*
************************************************************************
      subroutine SetMSbarMasses(mc,mb,mt)
*
      implicit none
*
      include "../commons/m2th.h"
      include "../commons/mass_scheme.h"
*
*     Variables
*
      double precision mc,mb,mt
*
      mass_scheme = "MSbar"
      m2th(4)     = mc * mc
      m2th(5)     = mb * mb
      m2th(6)     = mt * mt
      InMasses    = "done"
*
      return
      end
