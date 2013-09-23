************************************************************************
*
*     integralsQCD.f:
*
*     This function returns the convolution of the splitting functions
*     with the interpolation functions for a given number of flavours
*     nf and for the the pair of grid indices (alpha,beta) for singlet
*     and non-singlet in QCD according to:
*
*     kk  =  1   2   3   4   5   6   7
*            +   -   V   qq  qg  gq  gg
* 
************************************************************************
      function integralsQCD(alpha,beta,coup,kk)
*
      implicit none
*
      include "../commons/ipt.h"
      include "../commons/grid.h"
      include "../commons/integrals.h"
      include "../commons/wrap.h"
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
      double precision integralsQCD
*
      integralsQCD = 0d0
*
*     Return if it attempts to integrate for x > 1
*
      if(beta.ge.nin(igrid).or.alpha.ge.nin(igrid)) return
*
*     Integrals
*
      do pt=1,ipt+1
         integralsQCD = integralsQCD 
     1                + coup**pt * SP(igrid,wnf,kk,pt-1,alpha,beta)
      enddo
*
      return
      end
