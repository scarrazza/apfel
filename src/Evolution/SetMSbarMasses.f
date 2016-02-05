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
      m2ph(4)     = mc * mc
      m2ph(5)     = mb * mb
      m2ph(6)     = mt * mt
      InMasses    = "done"
*
      m2q(4) = m2ph(4)
      m2q(5) = m2ph(5)
      m2q(6) = m2ph(6)
*
      return
      end
