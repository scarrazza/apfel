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
**
*     Variables
*
      integer alpha,beta,tau
      integer la
*
*     Initialize integrals
*
      if(IsExt(igrid))then
         do tau=0,na
            do alpha=0,nin(igrid)-1
               do beta=alpha,nin(igrid)-1
                  call RSLintegralsQCDRes(LogAcc,alpha,beta,tau)
               enddo
            enddo
         enddo
         if((PDFEvol.eq."truncated".or.PDFEvol.eq."expandalpha").and.
     1        LogAcc.gt.0)then
            do la=0,LogAcc-1
               call HELLLogOrder(la)
               do tau=0,na
                  do alpha=0,nin(igrid)-1
                     do beta=alpha,nin(igrid)-1
                        call RSLintegralsQCDRes(la,alpha,beta,tau)
                     enddo
                  enddo
               enddo
            enddo
         endif
      else
         do tau=0,na
            do alpha=0,nin(igrid)-1
               call RSLintegralsQCDRes(LogAcc,0,alpha,tau)
            enddo
         enddo
         if((PDFEvol.eq."truncated".or.PDFEvol.eq."expandalpha").and.
     1        LogAcc.gt.0)then
            do la=0,LogAcc-1
               call HELLLogOrder(la)
               do tau=0,na
                  do alpha=0,nin(igrid)-1
                     call RSLintegralsQCDRes(la,0,alpha,tau)
                  enddo
               enddo
            enddo
         endif
      endif
*
      return
      end
