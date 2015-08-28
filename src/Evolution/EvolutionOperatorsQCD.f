************************************************************************
*
*     EvolutionOperatorsQCD.f:
*
*     This routine returns the singlet and the non-singlet evolution
*     operators on the x-space grid between the scales m20 and mu2 for
*     nf active flavours given the initial evolution operators in QCD.
*
************************************************************************
      subroutine EvolutionOperatorsQCD(muF20,muF2)
*
      implicit none
*
      include "../commons/Evs.h"
      include "../commons/Nf_FF.h"
      include "../commons/ipt.h"
      include "../commons/m2th.h"
      include "../commons/grid.h"
      include "../commons/EvolutionMatrices.h"
      include "../commons/wrap.h"
      include "../commons/MaxFlavourPDFs.h"
      include "../commons/MaxFlavourAlpha.h"
**
*     Input Variables
*
      double precision muF20,muF2
**
*     Internal Variables
*
      integer inf,nfmax
      double precision mu2i(3:7),mu2f(3:7)
      double precision M0sg(2,2,0:nint_max,0:nint_max)
      double precision M0nsp(0:nint_max,0:nint_max)
      double precision M0nsm(0:nint_max,0:nint_max)
      double precision M0nsv(0:nint_max,0:nint_max)
      double precision Msg(2,2,0:nint_max,0:nint_max)
      double precision Mnsp(0:nint_max,0:nint_max)
      double precision Mnsm(0:nint_max,0:nint_max)
      double precision Mnsv(0:nint_max,0:nint_max)
      double precision tiny
      parameter(tiny=1d-14)
*
*     Initial Conditions (Unity matrices)
*
      call IdentityOperatorsQCD(M0sg,M0nsp,M0nsm,M0nsv)
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
         nfi = Nf_FF
         nff = Nf_FF
         wnf = Nf_FF
*
         sgn = 1
*     If initial and final energies are equal return immediately the intial conditions
         if(muF2.eq.muF20)then
            call EqualOperatorsQCDnf(Nf_FF,M0sg,M0nsp,M0nsm,M0nsv,
     1                               MQCDsg,MQCDnsp,MQCDnsm,MQCDnsv)
            return
         endif
*     Singlet
         call odeintsgQCD(muF20,muF2,M0sg,Msg)
*     Non-Singlet
         if(ipt.eq.0)then
            call odeintnsQCD(1,muF20,muF2,M0nsp,Mnsp)
            call EqualOperatorsQCDnf(Nf_FF,Msg,Mnsp,Mnsp,Mnsp,
     1                               MQCDsg,MQCDnsp,MQCDnsm,MQCDnsv)
         elseif(ipt.eq.1)then
            call odeintnsQCD(1,muF20,muF2,M0nsp,Mnsp)
            call odeintnsQCD(2,muF20,muF2,M0nsm,Mnsm)
            call EqualOperatorsQCDnf(Nf_FF,Msg,Mnsp,Mnsm,Mnsm,
     1                               MQCDsg,MQCDnsp,MQCDnsm,MQCDnsv)
         elseif(ipt.eq.2)then
            call odeintnsQCD(1,muF20,muF2,M0nsp,Mnsp)
            call odeintnsQCD(2,muF20,muF2,M0nsm,Mnsm)
            call odeintnsQCD(3,muF20,muF2,M0nsv,Mnsv)
            call EqualOperatorsQCDnf(Nf_FF,Msg,Mnsp,Mnsm,Mnsv,
     1                               MQCDsg,MQCDnsp,MQCDnsm,MQCDnsv)
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
*     If initial and final energies are equal return immediately the intial conditions
         if(muF2.eq.muF20)then
            call EqualOperatorsQCDnf(nfi,M0sg,M0nsp,M0nsm,M0nsv,
     1                               MQCDsg,MQCDnsp,MQCDnsm,MQCDnsv)
            sgn = 1
            return
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
            wnf = inf
*     Singlet
            call odeintsgQCD(mu2i(inf),mu2f(inf),M0sg,Msg)
*     Non-Singlet
            if(ipt.eq.0)then
               call odeintnsQCD(1,mu2i(inf),mu2f(inf),M0nsp,Mnsp)
               call EqualOperatorsQCDnf(inf,Msg,Mnsp,Mnsp,Mnsp,
     1                                  MQCDsg,MQCDnsp,MQCDnsm,MQCDnsv)
            elseif(ipt.eq.1)then
               call odeintnsQCD(1,mu2i(inf),mu2f(inf),M0nsp,Mnsp)
               call odeintnsQCD(2,mu2i(inf),mu2f(inf),M0nsm,Mnsm)
               call EqualOperatorsQCDnf(inf,Msg,Mnsp,Mnsm,Mnsm,
     1                               MQCDsg,MQCDnsp,MQCDnsm,MQCDnsv)
            elseif(ipt.eq.2)then
               call odeintnsQCD(1,mu2i(inf),mu2f(inf),M0nsp,Mnsp)
               call odeintnsQCD(2,mu2i(inf),mu2f(inf),M0nsm,Mnsm)
               call odeintnsQCD(3,mu2i(inf),mu2f(inf),M0nsv,Mnsv)
               call EqualOperatorsQCDnf(inf,Msg,Mnsp,Mnsm,Mnsv,
     1                                  MQCDsg,MQCDnsp,MQCDnsm,MQCDnsv)
            endif
         enddo
      endif
*
      return
      end
