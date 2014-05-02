************************************************************************
*
*     initIntegralsQCDRes.f:
*
*     This routine initializes the integrals of small-x resummed splitting
*     functions and interpolation functions on the grid in alphas.
*
************************************************************************
      subroutine initIntegralsQCDRes
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/gridAlpha.h"
**
*     Variables
*
      integer alpha,tau
*
      do tau=0,na
         do alpha=0,nin(igrid)-1
            call RSLintegralsQCDRes(0,alpha,tau)
         enddo
      enddo
*
      return
      end
