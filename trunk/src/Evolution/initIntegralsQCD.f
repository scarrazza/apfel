************************************************************************
*
*     initIntegralsQCD.f:
*
*     This routine initializes the integrals of splitting functions and
*     and interpolation functions for a given number of active flavours 
*     nf in QCD.
*
************************************************************************
      subroutine initIntegralsQCD(nf)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/TimeLike.h"
**
*     Input Variables
*
      integer nf
**
*     Internal Variables
*
      integer alpha
*
      if(TimeLike)then
         do alpha=0,nin(igrid)-1
            call RSLintegralsQCDT(nf,0,alpha)
         enddo
      else
         do alpha=0,nin(igrid)-1
            call RSLintegralsQCD(nf,0,alpha)
         enddo
      endif
*
      return
      end
