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
      function integralsMatching(nf,alpha,beta,coup,kk,sgn)
*
      implicit none
*
      include "../commons/ColorFactors.h"
      include "../commons/ipt.h"
      include "../commons/grid.h"
      include "../commons/integrals.h"
      include "../commons/TimeLike.h"
      include "../commons/m2th.h"
**
*     Input Variables
*
      integer nf,alpha,beta,kk,sgn
      double precision coup
**
*     Internal Variables
*
      integer pt,ptstep
      integer s
      double precision fact
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
      s      = 1
      if(TimeLike.or.k2th(nf).ne.1d0) ptstep = 1
      do pt=0,ipt,ptstep
         if(pt.ne.0) s = sgn
         integralsMatching = integralsMatching 
     1                     + s * coup**pt
     2                     * SM(igrid,nf,kk,pt,alpha,beta)
      enddo
*
*     In case of backwards evolution, apart from changing the sign of NLO and NNLO terms,
*     an additional term to the NNLO should be added, because:
*
*     ( 1 + a A1 + a^2 A2 )^{-1} = 1 - a A1 - a^2 ( A2 - A1^2 ) + O(a^3)
*                                = 1 - a A1 - a^2 A2 + [a^2 A1^2] + O(a^3)
*
*     but being A1 different from zero only for the singlet with the form:
*
*     A1 = | 0  A1qg |  => A1^2 = A1gg * A1
*          | 0  A1gg |
*
*     with:
*
*     A1gg = - 4 / 3 TR * delta(1-x) ln(k2th)
*
      if(.not.TimeLike.and.
     1     sgn.eq.-1.and.ipt.ge.2.and.k2th(nf).ne.1d0.and.
     2     (kk.eq.3.or.kk.eq.5).and.
     3     alpha.eq.beta)then
         fact = - 4d0 * TR * dlog(k2th(nf)) / 3d0 ! A1gg
         integralsMatching = integralsMatching
     1        + fact * coup**2 * SM(igrid,nf,kk,1,alpha,beta)
      endif
*
      return
      end
