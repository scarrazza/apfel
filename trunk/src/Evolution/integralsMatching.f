************************************************************************
*
*     integralsMatching.f:
*
*     This function returns the convolution of the matching conditions
*     with the interpolation functions for the the pair of grid indices 
*     (alpha,beta) for singlet and non-singlet combinations:
*
*     kk  =  1   2   3   4   5 
*            V   qq  qg  gq  gg
* 
************************************************************************
      function integralsMatching(alpha,beta,coup,kk)
*
      implicit none
*
      include "../commons/ipt.h"
      include "../commons/grid.h"
      include "../commons/integrals.h"
      include "../commons/TimeLike.h"
**
*     Input Variables
*
      integer alpha,beta,kk
      integer ptstep
      double precision coup
**
*     Internal Variables
*
      integer pt
**
*     Output Variables
*
      double precision integralsMatching
*
      integralsMatching = 0d0
*
*     Return if it attempts to integrate for x > 1
*
      if(beta.ge.nin(igrid).or.alpha.ge.nin(igrid)) return
*
*     Integrals
*
*     For the space-like evolution one can skip the order alphas because it is
*     always equal to zero but this can't be done for the time-like evolution 
*
      ptstep = 2
      if(TimeLike) ptstep = 1
      do pt=0,ipt,ptstep
         integralsMatching = integralsMatching 
     1                     + coup**pt * SM(igrid,kk,pt,alpha,beta)
      enddo
*
      return
      end
