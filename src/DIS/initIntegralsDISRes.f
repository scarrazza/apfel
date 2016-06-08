************************************************************************
*
*     initIntegralsDISRes.f:
*
*     This routine initializes the integrals of small-x resummed coefficient
*     functions and interpolation functions on the grid in alphas.
*
************************************************************************
      subroutine initIntegralsDISRes
*
      implicit none
*
      include "../commons/Smallx.h"
      include "../commons/grid.h"
      include "../commons/gridAlpha.h"
**
*     Variables
*
      integer alpha,beta,tau
*
*     Initialize integrals 
*
      call HELLLogOrder(LogAcc)
      if(IsExt(igrid))then
         do tau=0,na
            do alpha=0,nin(igrid)-1
               do beta=alpha,nin(igrid)-1
                  call RSLintegralsDISRes(alpha,beta,tau)
               enddo
            enddo
         enddo
      else
         do tau=0,na
            do alpha=0,nin(igrid)-1
               call RSLintegralsDISRes(0,alpha,tau)
            enddo
         enddo
      endif
*
      return
      end
