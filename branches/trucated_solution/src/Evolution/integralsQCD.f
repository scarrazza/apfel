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
      include "../commons/PDFEvolution.h"
      include "../commons/integrals.h"
      include "../commons/wrap.h"
      include "../commons/Smallx.h"
      include "../commons/gridAlpha.h"
      include "../commons/integralsRes.h"
      include "../commons/EpsTrunc.h"
**
*     Input Variables
*
      integer alpha,beta,tau,kk
      double precision coup
      double precision fbeta,beta0,beta1,beta2,b1,b2
      double precision c1,c2
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
      if(PDFEvol.eq."exactmu")then
         do pt=1,ipt+1
            integralsQCD = integralsQCD 
     1                   + coup**pt * SP(igrid,wnf,kk,pt-1,alpha,beta)
     2                   * EpsTrunc**(pt-1)
         enddo
*
         if(Smallx.and.kk.ge.4)then
*     find "tau" such that ag(tau) <= coup < ag(tau+1)
            do tau=0,na-1
               if(nfg(tau).eq.wnf.and.
     1            ag(tau).ge.coup.and.ag(tau+1).lt.coup) goto 102
            enddo
*     interpolation coefficient
 102        c1 = ( ag(tau+1) - coup ) / ( ag(tau+1) - ag(tau) )
            c2 = ( coup - ag(tau) ) / ( ag(tau+1) - ag(tau) )
*
            integralsQCD = integralsQCD 
     1                   + c1 * SPRes(igrid,kk,alpha,beta,tau)
     2                   + c2 * SPRes(igrid,kk,alpha,beta,tau+1)
         endif
      elseif(PDFEvol.eq."exactalpha")then
         do pt=1,ipt+1
            integralsQCD = integralsQCD 
     1                   + coup**pt * SP(igrid,wnf,kk,pt-1,alpha,beta)
         enddo
         if(Smallx.and.kk.ge.4)then
*     find "tau" such that ag(tau) <= coup < ag(tau+1)
            do tau=0,na-1
               if(nfg(tau).eq.wnf.and.
     1            ag(tau).ge.coup.and.ag(tau+1).lt.coup) goto 103
            enddo
*     interpolation coefficient
 103        c1 = ( ag(tau+1) - coup ) / ( ag(tau+1) - ag(tau) )
            c2 = ( coup - ag(tau) ) / ( ag(tau+1) - ag(tau) )
*
            integralsQCD = integralsQCD 
     1                   + c1 * SPRes(igrid,kk,alpha,beta,tau)
     2                   + c2 * SPRes(igrid,kk,alpha,beta,tau+1)
         endif
         integralsQCD = integralsQCD / fbeta(coup,wnf,ipt)
      elseif(PDFEvol.eq."expandalpha")then
         integralsQCD = SP(igrid,wnf,kk,0,alpha,beta)
         if(ipt.ge.1)then
            b1 = beta1(wnf) / beta0(wnf)
            integralsQCD = integralsQCD 
     1                   + coup * ( SP(igrid,wnf,kk,1,alpha,beta) 
     2                   - b1 * SP(igrid,wnf,kk,0,alpha,beta) )
         endif
         if(ipt.ge.2)then
            b2 = beta2(wnf) / beta0(wnf)
            integralsQCD = integralsQCD 
     1           + coup**2d0 * ( SP(igrid,wnf,kk,2,alpha,beta) 
     2           - b1 * SP(igrid,wnf,kk,1,alpha,beta)
     3           + ( b1**2d0 - b2 ) * SP(igrid,wnf,kk,0,alpha,beta) )
         endif
         integralsQCD = - integralsQCD / beta0(wnf) / coup
      endif
*
      return
      end
