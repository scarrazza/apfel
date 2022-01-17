************************************************************************
*
*     RSLintegralspDIS.f:
*
*     This routine evaluates and dump on the common integralsDIS.h once and 
*     for all for the pair (alpha,beta) the integral of the sum of Regular
*     and Singular part plus the Local part of all the coefficient functions
*     for all the needed orders in polarised DIS.
*     The index sf runs like that:
*
*        sf  Coefficient Function
*     --------------------------
*        1           G1
*        2           GL
*        3           G4

*
*     while the index k runs like that:
*
*        k      combination
*     --------------------------
*        1         gluon  
*        2      pure-singlet
*        3    non-singlet-plus
*        4    non-singlet-minus
*
************************************************************************
      subroutine RSLintegralspDIS(beta,alpha)
*
      implicit none
*
c      include "../commons/ipt.h"
      include "../commons/wrapDIS.h"
      include "../commons/grid.h"
      include "../commons/coeffhqmellin.h"
      include "../commons/integralsDIS.h"
      include "../commons/Nf_FF.h"
      include "../commons/ColorFactors.h"
**
*     Input Variables
*
      integer beta,alpha
**
*     Internal Variables
*
      integer bound,inf,ipt_ZM
      double precision G4L(0:2),GLL(0:2),G1L(0:2)
      double precision fL
      double precision dgauss,a,b,eps(2)
      double precision integrandspDISzm
      double precision integG4(0:2),integGL(0:2),integG1(0:2)
      double precision G4ns1RS,G4ns1L,GLns1RS,G1ns1RS,G1ns1L,G1g1R
      double precision G4NS1PC,G1NS1PC
      external integrandspDISzm
c      data eps / 5d-8, 1d-3 /
c      data eps / 5d-8, 1d-5 /
      data eps / 1d-5, 1d-3 /
c      data eps / 1d-7, 1d-5 /
*
*     Adjustment of the bounds of the integrals
*
      if(alpha.lt.beta) return
*
      bound = alpha-inter_degree(igrid)
      if(alpha.lt.inter_degree(igrid)) bound = 0
      a = max(xg(igrid,beta),xg(igrid,beta)/xg(igrid,alpha+1))
      b = min(1d0,xg(igrid,beta)/xg(igrid,bound))
*
      fL = 0d0
      if(alpha.eq.beta) fL = 1d0
*
*     The ZM scheme is cheap enough to compute all orders anyway. This
*     is useful because in some cases one may need different absolute
*     orderes at the same time, like with F2 and FL.
*
      ipt_ZM = 2
*
*     Initialize Integrals
*
      do inf=3,6
         do k=1,4
            do wipt=0,ipt_ZM
               SC2zm(igrid,inf,k,wipt,beta,alpha) = 0d0
               SCLzm(igrid,inf,k,wipt,beta,alpha) = 0d0
               SC3zm(igrid,inf,k,wipt,beta,alpha) = 0d0
            enddo
         enddo
*
*     Variables needed for wrapping the integrand functions
*
         wnf    = inf
         walpha = alpha
         wbeta  = beta
*
*     Precompute integrals
*
         if(ipt_ZM.ge.1)then
            wipt = 1
*
            sf = 3
            k = 3
            G4ns1RS = dgauss(integrandspDISzm,a,b,eps(wipt))
            G4ns1L  = G4NS1PC(a)
*
            sf = 2
            k = 3
            GLns1RS = dgauss(integrandspDISzm,a,b,eps(wipt))
*
            sf = 1
            k = 1
            G1g1R = dgauss(integrandspDISzm,a,b,eps(wipt))
            k = 3
            G1ns1RS = dgauss(integrandspDISzm,a,b,eps(wipt))
            G1ns1L  = G1NS1PC(a)
         endif
*
         do k=1,4
*
*     LO
*
*     G4
            G4L(0) = 0d0
            if(k.eq.3.or.k.eq.4) G4L(0) = 1d0
            integG4(0) = 0d0
*     GL
            GLL(0) = 0d0
            integGL(0) = 0d0
*     G1
            G1L(0) = 0d0
            if(k.eq.3.or.k.eq.4) G1L(0) = 1d0
            integG1(0) = 0d0
*
*     NLO
*
            if(ipt_ZM.ge.1)then
*     Gluon
               if(k.eq.1)then
*     G4
                  G4L(1) = 0d0
                  integG4(1) = 0d0
*     GL
                  GLL(1) = 0d0
                  integGL(1) = 0d0
*     G1
                  G1L(1) = 0d0
                  integG1(1) = G1g1R
*     Pure-Singlet
               elseif(k.eq.2)then
*     G4
                  G4L(1) = 0d0
                  integG4(1) = 0d0
*     GL
                  GLL(1) = 0d0
                  integGL(1) = 0d0
*     G1
                  G1L(1) = 0d0
                  integG1(1) = 0d0
*     Non-singlet-plus/minus
               elseif(k.eq.3.or.k.eq.4)then
*     G4
                  G4L(1) = G4ns1L
                  integG4(1) = G4ns1RS
*     GL
                  GLL(1) = 0d0
                  integGL(1) = GLns1RS
*     G1
                  G1L(1) = G1ns1L
                  integG1(1) = G1ns1RS
               endif
            endif
*
*     Integrals
*
            do wipt=0,ipt_ZM
               SC2zm(igrid,inf,k,wipt,beta,alpha) = integG1(wipt) 
     1                                            + G1L(wipt) * fL
               SCLzm(igrid,inf,k,wipt,beta,alpha) = integGL(wipt) 
     1                                            + GLL(wipt) * fL
               SC3zm(igrid,inf,k,wipt,beta,alpha) = integG4(wipt) 
     1                                            + G4L(wipt) * fL
            enddo
         enddo
      enddo
*
      return
      end
