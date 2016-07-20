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
      include "../commons/Smallx.h"
      include "../commons/PDFEvolution.h"
      include "../commons/grid.h"
      include "../commons/gridAlpha.h"
      include "../commons/ipt.h"
**
*     Variables
*
      integer la,lamin,lamax
      integer alpha,beta,tau
*
*     According to the DGLAP solution initialize the right
*     integrals
*
      if(PDFEvol.eq."truncated".or.
     1   PDFEvol.eq."expandalpha")then
         lamin = 0
      else
         lamin = min(LogAcc,1)
      endif
      lamax = min(LogAcc,1)
*
*     Initialize integrals
*
      if(IsExt(igrid))then
         do la=lamin,lamax
            call initHELL(la,ipt)
            do tau=0,na
               do alpha=0,nin(igrid)-1
                  do beta=alpha,nin(igrid)-1
                     call RSLintegralsQCDRes(la,alpha,beta,tau)
                  enddo
               enddo
            enddo
         enddo
      else
         do la=lamin,lamax
            call initHELL(la,ipt)
            do tau=0,na
               do alpha=0,nin(igrid)-1
                  call RSLintegralsQCDRes(la,0,alpha,tau)
               enddo
            enddo
         enddo
      endif
*
      return
      end
