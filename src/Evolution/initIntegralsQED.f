************************************************************************
*
*     initIntegralsQED.f:
*
*     This routine initializes the integrals of splitting functions and
*     interpolation functions for a given number of active flavours nf 
*     in QCD.
*
************************************************************************
      subroutine initIntegralsQED(nf,nl)
*
      implicit none
*
      include "../commons/grid.h"
**
*     Input Variables
*
      integer nf,nl
**
*     Internal Variables
*
      integer alpha,beta
*
      if(IsExt(igrid))then
*
*     If this is an external grid, compute the integrals for
*     the entire splitting matrix ...
*
         do alpha=0,nin(igrid)-1
            do beta=alpha,nin(igrid)-1
               call RSLintegralsQED(nf,nl,alpha,beta)
            enddo
         enddo
      else
*
*     ... otherwise only for the first line
*
         do alpha=0,nin(igrid)-1
            call RSLintegralsQED(nf,nl,0,alpha)
         enddo
      endif
*
      return
      end
