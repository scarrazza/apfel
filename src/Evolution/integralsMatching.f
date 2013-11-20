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
**
*     Input Variables
*
      integer alpha,beta,kk
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
      do pt=0,ipt,2
         integralsMatching = integralsMatching 
     1                     + coup**pt * SM(igrid,kk,pt,alpha,beta)
      enddo
*
      return
      end
