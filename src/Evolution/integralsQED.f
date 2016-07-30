************************************************************************
*
*     integralsQED.f:
*
*     This function returns the convolution of the splitting functions
*     with the interpolation functions for a given number of flavours
*     nf and for the the pair of grid indices (alpha,beta) for singlet
*     and non-singlet in QED according to:
*
*     kk  =  1    2    3    4
*            nspu nsmu nspd nsmd
*
*            5    6    7    8
*            gg   ggm  gS   gD
*
*            9    10   11   12
*            gmg  gmgm gmS  gmD
*
*            13   14   15   16
*            Sg   Sgm  SS   SD
*
*            17   18   19   20
*            Dg   Dgm  DS   DD
*
*            21   22
*            VV   VDV (DVDV = VV, DVV = VDV)
*
*            23   24   25
*            LL   gmL  Lgm
* 
************************************************************************
      function integralsQED(alpha,beta,coupQED,coupQCD,kk)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/PDFEvolution.h"
      include "../commons/integrals.h"
      include "../commons/wrap.h"
      include "../commons/Th.h"
      include "../commons/MaxFlavourPDFs.h"
      include "../commons/MaxFlavourAlpha.h"
      include "../commons/ipt.h"
      include "../commons/NLOQEDCorrections.h"
**
*     Input Variables
*
      integer alpha,beta,kk
      integer bnf,pnf
      integer jpt
      double precision coupQED,coupQCD
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
      jpt = 0
      if(NLOQED.and.ipt.ge.1) jpt = 2
*
*     Define the number of flavours to be used in the
*     beta function and in the splitting functions.
*
      bnf = min(wnf,nfMaxAlpha)
      pnf = min(wnf,nfMaxPDFs)
*
*     Integrals
*
      integralsQED = coupQED * SQ(igrid,pnf,wnl,kk,0,alpha,beta)
      if(jpt.ge.1) integralsQED = integralsQED
     1                          + coupQCD * coupQED
     2                          * SQ(igrid,pnf,wnl,kk,1,alpha,beta)
      if(jpt.ge.2) integralsQED = integralsQED
     1                          + coupQED**2
     2                          * SQ(igrid,pnf,wnl,kk,2,alpha,beta)
*
      return
      end
