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
      include "../commons/ipt.h"
      include "../commons/PDFEvolution.h"
      include "../commons/grid.h"
      include "../commons/gridAlpha.h"
**
*     Variables
*
      integer pt,ptmin,ptmax
      integer alpha,beta,tau
*
*     According to the DGLAP solution initialize the right
*     integrals
*
      if(PDFEvol.eq."truncated")then
         ptmin = 0
      else
         ptmin = min(ipt,1)
      endif
      ptmax = min(ipt,1)
*
*     Initialize integrals 
*
      if(IsExt(igrid))then
         do pt=ptmin,ptmax
            call HELLOrder(pt)
            do tau=0,na
               do alpha=0,nin(igrid)-1
                  do beta=alpha,nin(igrid)-1
                     call RSLintegralsQCDRes(pt,alpha,beta,tau)
                  enddo
               enddo
            enddo
         enddo
      else
         do pt=ptmin,ptmax
            call HELLOrder(pt)
            do tau=0,na
               do alpha=0,nin(igrid)-1
                  call RSLintegralsQCDRes(pt,0,alpha,tau)
               enddo
            enddo
         enddo
      endif
*
      return
      end
