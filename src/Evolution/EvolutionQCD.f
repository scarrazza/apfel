************************************************************************
*
*     EvolutionQCD.f:
*
*     This routine returns the singlet and the non-singlet evolution
*     operators on the x-space grid between the scales m20 and mu2 for
*     nf active flavours given the initial evolution operators in QCD.
*
************************************************************************
      subroutine EvolutionQCD(muF20,muF2)
*
      implicit none
*
      include "../commons/Evs.h"
      include "../commons/Nf_FF.h"
      include "../commons/m2th.h"
      include "../commons/grid.h"
      include "../commons/wrap.h"
      include "../commons/MaxFlavourPDFs.h"
      include "../commons/MaxFlavourAlpha.h"
      include "../commons/f0ph.h"
      include "../commons/fph.h"
**
*     Input Variables
*
      double precision muF20,muF2
**
*     Internal Variables
*
      integer inf,jnf
      integer i,alpha
      integer sgn
      integer nfi,nff,nfmax
      double precision mu2i(3:7),mu2f(3:7)
      double precision f0ev(0:13,0:nint_max)
      double precision fev(0:13,0:nint_max)
      double precision fpy(-6:6,0:nint_max)
      double precision Sg(2,0:nint_max),Sg0(2,0:nint_max)
      double precision V(0:nint_max),V0(0:nint_max)
      double precision V3(0:nint_max),V30(0:nint_max)
      double precision V8(0:nint_max),V80(0:nint_max)
      double precision V15(0:nint_max),V150(0:nint_max)
      double precision V24(0:nint_max),V240(0:nint_max)
      double precision V35(0:nint_max),V350(0:nint_max)
      double precision T3(0:nint_max),T30(0:nint_max)
      double precision T8(0:nint_max),T80(0:nint_max)
      double precision T15(0:nint_max),T150(0:nint_max)
      double precision T24(0:nint_max),T240(0:nint_max)
      double precision T35(0:nint_max),T350(0:nint_max)
      double precision tiny
      parameter(tiny=1d-10)
*
*     Initial Conditions
*
      call initPDFs(muF20)
*
*     Rotate initial PDFs into the evolution basis
*
      call PDFphys2evQCD(f0ph,f0ev)
*
*     Define maximun number of flavours
*
      nfmax = max(nfMaxPDFs,nfMaxAlpha)
*
*     Mass scheme
*
*     Fixed Flavour Number Scheme
*
      if(Evs.eq."FF")then
         wnf = Nf_FF
*     If initial and final energies are equal,
*     return immediately the intial conditions
         if(muF2.eq.muF20)then
            do alpha=0,nin(igrid)
               do i=0,13
                  fev(i,alpha) = f0ev(i,alpha)
               enddo
            enddo
            goto 104
         endif
*
         do alpha=0,nin(igrid)
            Sg0(1,alpha) = f0ev(1,alpha)
            Sg0(2,alpha) = f0ev(2,alpha)
            V0(alpha)    = f0ev(3,alpha)
            V30(alpha)   = f0ev(4,alpha)
            V80(alpha)   = f0ev(5,alpha)
            V150(alpha)  = f0ev(6,alpha)
            V240(alpha)  = f0ev(7,alpha)
            V350(alpha)  = f0ev(8,alpha)
            T30(alpha)   = f0ev(9,alpha)
            T80(alpha)   = f0ev(10,alpha)
            T150(alpha)  = f0ev(11,alpha)
            T240(alpha)  = f0ev(12,alpha)
            T350(alpha)  = f0ev(13,alpha)
         enddo
*     Singlet
         call odeintsgQCDf(muF20,muF2,Sg0,Sg)
*     Non-Singlet
         call odeintnsQCDf(3,muF20,muF2,V0,V)
         call odeintnsQCDf(2,muF20,muF2,V30,V3)
         call odeintnsQCDf(2,muF20,muF2,V80,V8)
         call odeintnsQCDf(1,muF20,muF2,T30,T3)
         call odeintnsQCDf(1,muF20,muF2,T80,T8)
         if(Nf_FF.eq.3)then
            do alpha=0,nin(igrid)
               V15(alpha) = V(alpha)
               V24(alpha) = V(alpha)
               V35(alpha) = V(alpha)
               T15(alpha) = Sg(1,alpha)
               T24(alpha) = Sg(1,alpha)
               T35(alpha) = Sg(1,alpha)
            enddo
         endif
         if(Nf_FF.eq.4)then
            call odeintnsQCDf(2,muF20,muF2,V150,V15)
            call odeintnsQCDf(1,muF20,muF2,T150,T15)
            do alpha=0,nin(igrid)
               V24(alpha) = V(alpha)
               V35(alpha) = V(alpha)
               T24(alpha) = Sg(1,alpha)
               T35(alpha) = Sg(1,alpha)
            enddo
         endif
         if(Nf_FF.eq.5)then
            call odeintnsQCDf(2,muF20,muF2,V150,V15)
            call odeintnsQCDf(2,muF20,muF2,V240,V24)
            call odeintnsQCDf(1,muF20,muF2,T150,T15)
            call odeintnsQCDf(1,muF20,muF2,T240,T24)
            do alpha=0,nin(igrid)
               V35(alpha) = V(alpha)
               T35(alpha) = Sg(1,alpha)
            enddo
         endif
         if(Nf_FF.eq.6)then
            call odeintnsQCDf(2,muF20,muF2,V150,V15)
            call odeintnsQCDf(2,muF20,muF2,V240,V24)
            call odeintnsQCDf(2,muF20,muF2,V350,V35)
            call odeintnsQCDf(1,muF20,muF2,T150,T15)
            call odeintnsQCDf(1,muF20,muF2,T240,T24)
            call odeintnsQCDf(1,muF20,muF2,T350,T35)
         endif
*
*     Variable Flavour Number Scheme
*
      elseif(Evs.eq."VF")then
*
*     Find initial and final number of flavours
*
         if(muF2.gt.m2th(6))then
            nff = 6
         elseif(muF2.gt.m2th(5))then
            nff = 5
         elseif(muF2.gt.m2th(4))then
            nff = 4
         else
            nff = 3
         endif
         if(nff.gt.nfmax) nff = nfmax
*
         if(muF20.gt.m2th(6))then
            nfi = 6
         elseif(muF20.gt.m2th(5))then
            nfi = 5
         elseif(muF20.gt.m2th(4))then
            nfi = 4
         else
            nfi = 3
         endif
         if(nfi.gt.nfmax) nfi = nfmax
*     If initial and final energies are equal, return immediately the intial conditions
         if(muF2.eq.muF20)then
            do alpha=0,nin(igrid)
               do i=0,13
                  fev(i,alpha) = f0ev(i,alpha)
               enddo
            enddo
            goto 104
         elseif(muF2.gt.muF20)then
            sgn = 1
         elseif(muF2.lt.muF20)then
            sgn = - 1
         endif
*
         mu2i(nfi) = muF20
         if(sgn.eq.1)then
            do inf=nfi+1,nff
               mu2i(inf) = m2th(inf)
            enddo
            do inf=nfi,nff-1
               mu2f(inf) = m2th(inf+1) - tiny
            enddo
         elseif(sgn.eq.-1)then
            do inf=nfi-1,nff,sgn
               mu2i(inf) = m2th(inf+1) + tiny
            enddo
            do inf=nfi,nff+1,sgn
               mu2f(inf) = m2th(inf)
            enddo
         endif
         mu2f(nff) = muF2
*
         do inf=nfi,nff,sgn
            jnf = min(inf,nfMaxPDFs)
            wnf = inf
*
*     Apply matching conditions for the backward evolution
*
            if(sgn.eq.-1.and.inf.lt.nfi.and.inf.lt.nfMaxPDFs)then
               call MatchPDFs(inf+1,sgn,f0ev)
            endif
*
            do alpha=0,nin(igrid)
               Sg0(1,alpha) = f0ev(1,alpha)
               Sg0(2,alpha) = f0ev(2,alpha)
               V0(alpha)    = f0ev(3,alpha)
               V30(alpha)   = f0ev(4,alpha)
               V80(alpha)   = f0ev(5,alpha)
               V150(alpha)  = f0ev(6,alpha)
               V240(alpha)  = f0ev(7,alpha)
               V350(alpha)  = f0ev(8,alpha)
               T30(alpha)   = f0ev(9,alpha)
               T80(alpha)   = f0ev(10,alpha)
               T150(alpha)  = f0ev(11,alpha)
               T240(alpha)  = f0ev(12,alpha)
               T350(alpha)  = f0ev(13,alpha)
            enddo
*     Singlet
            call odeintsgQCDf(mu2i(inf),mu2f(inf),Sg0,Sg)
*     Non-Singlet
            call odeintnsQCDf(3,mu2i(inf),mu2f(inf),V0,V)
            call odeintnsQCDf(2,mu2i(inf),mu2f(inf),V30,V3)
            call odeintnsQCDf(2,mu2i(inf),mu2f(inf),V80,V8)
            call odeintnsQCDf(1,mu2i(inf),mu2f(inf),T30,T3)
            call odeintnsQCDf(1,mu2i(inf),mu2f(inf),T80,T8)
            if(jnf.eq.3)then
               do alpha=0,nin(igrid)
                  V15(alpha) = V(alpha)
                  V24(alpha) = V(alpha)
                  V35(alpha) = V(alpha)
                  T15(alpha) = Sg(1,alpha)
                  T24(alpha) = Sg(1,alpha)
                  T35(alpha) = Sg(1,alpha)
               enddo
            endif
            if(jnf.eq.4)then
               call odeintnsQCDf(2,mu2i(inf),mu2f(inf),V150,V15)
               call odeintnsQCDf(1,mu2i(inf),mu2f(inf),T150,T15)
               do alpha=0,nin(igrid)
                  V24(alpha) = V(alpha)
                  V35(alpha) = V(alpha)
                  T24(alpha) = Sg(1,alpha)
                  T35(alpha) = Sg(1,alpha)
               enddo
            endif
            if(jnf.eq.5)then
               call odeintnsQCDf(2,mu2i(inf),mu2f(inf),V150,V15)
               call odeintnsQCDf(2,mu2i(inf),mu2f(inf),V240,V24)
               call odeintnsQCDf(1,mu2i(inf),mu2f(inf),T150,T15)
               call odeintnsQCDf(1,mu2i(inf),mu2f(inf),T240,T24)
               do alpha=0,nin(igrid)
                  V35(alpha) = V(alpha)
                  T35(alpha) = Sg(1,alpha)
               enddo
            endif
            if(jnf.eq.6)then
               call odeintnsQCDf(2,mu2i(inf),mu2f(inf),V150,V15)
               call odeintnsQCDf(2,mu2i(inf),mu2f(inf),V240,V24)
               call odeintnsQCDf(2,mu2i(inf),mu2f(inf),V350,V35)
               call odeintnsQCDf(1,mu2i(inf),mu2f(inf),T150,T15)
               call odeintnsQCDf(1,mu2i(inf),mu2f(inf),T240,T24)
               call odeintnsQCDf(1,mu2i(inf),mu2f(inf),T350,T35)
            endif
*
*     Update initial coditions for the next step
*
            do alpha=0,nin(igrid)
               f0ev(1,alpha)  = Sg(1,alpha)
               f0ev(2,alpha)  = Sg(2,alpha)
               f0ev(3,alpha)  = V(alpha)
               f0ev(4,alpha)  = V3(alpha)
               f0ev(5,alpha)  = V8(alpha)
               f0ev(6,alpha)  = V15(alpha)
               f0ev(7,alpha)  = V24(alpha)
               f0ev(8,alpha)  = V35(alpha)
               f0ev(9,alpha)  = T3(alpha)
               f0ev(10,alpha) = T8(alpha)
               f0ev(11,alpha) = T15(alpha)
               f0ev(12,alpha) = T24(alpha)
               f0ev(13,alpha) = T35(alpha)
            enddo
*
*     Apply matching conditions for the forward evolution
*
            if(sgn.eq.1.and.inf.lt.nff.and.inf.lt.nfMaxPDFs)then
               call MatchPDFs(inf+1,sgn,f0ev)
            endif
         enddo
      endif
*
*     Put evolved PDFs in a single vector
*
      do alpha=0,nin(igrid)
         fev(1,alpha)  = Sg(1,alpha)
         fev(2,alpha)  = Sg(2,alpha)
         fev(3,alpha)  = V(alpha)
         fev(4,alpha)  = V3(alpha)
         fev(5,alpha)  = V8(alpha)
         fev(6,alpha)  = V15(alpha)
         fev(7,alpha)  = V24(alpha)
         fev(8,alpha)  = V35(alpha)
         fev(9,alpha)  = T3(alpha)
         fev(10,alpha) = T8(alpha)
         fev(11,alpha) = T15(alpha)
         fev(12,alpha) = T24(alpha)
         fev(13,alpha) = T35(alpha)
      enddo
*
*     Rotate PDFs back into the physical basis
*
 104  call PDFevQCD2phys(fev,fpy)
*
*     Put evolved PDFs in the physical basis in the common used for the interpolation
*
      do alpha=0,nin(igrid)
         do i=-6,6
            fph(igrid,i,alpha) = fpy(i,alpha)
         enddo
         fgamma(igrid,alpha) = f0lep(0,alpha)
      enddo
*
      return
      end
