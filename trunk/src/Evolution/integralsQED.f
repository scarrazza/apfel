************************************************************************
*
*     integralsQED.f:
*
*     This function returns the convolution of the splitting functions
*     with the interpolation functions for a given number of flavours
*     nf and for the the pair of grid indices (alpha,beta) for singlet
*     and non-singlet in QED according to:
*
*     kk  =  1   2   3   4   5   6   7
*            +   -   V   qq  qg  gq  gg
* 
************************************************************************
      function integralsQED(alpha,beta,coup,kk)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/PDFEvolution.h"
      include "../commons/integrals.h"
      include "../commons/wrap.h"
**
*     Input Variables
*
      integer alpha,beta,kk
      double precision coup
      double precision beta0qed
**
*     Output Variables
*
      double precision integralsQED
*
      integralsQED = 0d0
*
*     Return if it attempts to integrate for x > 1
*
      if(beta.ge.nin(igrid).or.alpha.ge.nin(igrid)) return
*
*     Integrals
*
      if(PDFEvol.eq."exactmu")then
         integralsQED = coup * SQ(igrid,wnf,kk,alpha,beta)
      else
         integralsQED = SQ(igrid,wnf,kk,alpha,beta) 
     1                / beta0qed(wnf) / coup
      endif
*
      return
      end
