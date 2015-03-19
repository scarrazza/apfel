************************************************************************
*
*     switchGluonPhoton.f:
*
*     This routine exchanges the gluon PDF with the Photon PDF at the
*     initial scale.
*
************************************************************************
      subroutine switchGluonPhoton
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/f0ph.h"
*
*     Variables
*
      integer alpha
      double precision fb(0:nint_max)
*
      do alpha=0,nin(igrid)
         fb(alpha)      = f0lep(0,alpha)
         f0lep(0,alpha) = f0ph(0,alpha)
         f0ph(0,alpha)  = fb(alpha)
      enddo
*
      return
      end
