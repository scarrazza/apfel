************************************************************************
*
*     IncludeIntrinsicCharm.f:
*
*     This routine includes in the precomputed coefficient functions
*     the intrinsic charm contributions to the massive coefficient
*     functions computed in RSLintegralsDIS.f.
*
************************************************************************
      subroutine IncludeIntrinsicCharm
*
      implicit none
*
      include "../commons/ipt.h"
      include "../commons/grid.h"
      include "../commons/coeffhqmellin.h"
      include "../commons/integralsDIS.h"
      include "../commons/MassScheme.h"
      include "../commons/Nf_FF.h"
**
*     Internal Variables
*
      integer ixi
      integer gbound
      integer alpha,beta
      double precision lambda,eta
      double precision w_int,win
*
*     if the ZM-VFNS has been selected no IC contribution
*     has to be included.
*
      if(MassScheme.eq."ZM-VFNS") return
*     
*     If an external grid is found compute the whole operator
*
      gbound = 0
      if(IsExt(igrid)) gbound = nin(igrid) - 1
*
*     LO
*
      do ixi=1,nxir
         lambda = 1d0 / xigrid(ixi*xistep)
         eta    = 2d0 / ( 1d0 + dsqrt( 1d0 + 4d0 * lambda ) )
         do beta=0,gbound
            do alpha=beta,nin(igrid)-1
*
*     Neutral current
*
*     FFNS
               win = w_int(inter_degree(igrid),alpha,
     1              xg(igrid,beta)/eta)
*
               SC2mNC(igrid,ixi,3,0,beta,alpha) = ( 2d0 - eta ) * win
               SCLmNC(igrid,ixi,3,0,beta,alpha) = 4d0 * ( 1d0 - eta )
     1              / ( 2d0 - eta ) * win
*     FFN0
               if(alpha.eq.beta) SC2m0NC(igrid,ixi,3,0,beta,alpha) = 1d0
*
*     Charged current
*
*     (even if the IC contribution is a non-singlet one, put it in the 
*     pure-singlet slot (e.g. SC2mCC(igrid,ixi,2,0,beta,alpha) because
*     for now this slot is never used. If one day there will be the 
*     O(as^2) correction to CC, that contain a pure-singlet piece, I
*     will need to extend the CC arrays from 3 to 4.)
*
               if(alpha.eq.beta)then
*     FFNS
                  SC2mCC(igrid,ixi,2,0,beta,alpha) = ( 1d0 + lambda )
                  SCLmCC(igrid,ixi,2,0,beta,alpha) = lambda
                  SC3mCC(igrid,ixi,2,0,beta,alpha) = 1d0
*     FFN0
                  SC2m0CC(igrid,ixi,2,0,beta,alpha) = 1d0
                  SC3m0CC(igrid,ixi,2,0,beta,alpha) = 1d0
               endif
            enddo
         enddo
      enddo
*
      return
      end
