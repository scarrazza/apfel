************************************************************************
*
*     SetMassMatchingScales.f:
*
*     This subroutine sets as a default the mass matching scales to
*     be a factor of the heavy quark heavy quark pole masses. 
*
************************************************************************
      subroutine SetMassMatchingScales(kmc,kmb,kmt)
*
      implicit none
*
      include "../commons/m2th.h"
*
*     Variables
*
      double precision kmc,kmb,kmt
*
      k2th(4)     = kmc * kmc
      k2th(5)     = kmb * kmb
      k2th(6)     = kmt * kmt
      InThrRatios = "done"
*
      return
      end
