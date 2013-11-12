************************************************************************
*
*     EvolutionOperatorsQED.f:
*
*     This routine returns the single and the non-singlet evolution
*     operators on the x-spage grid between the scales m20 and mu2 for
*     nf active flavours given the initial evolution operators in QED.
*
************************************************************************
      subroutine EvolutionOperatorsQED(muF20,muF2)
*
      implicit none
*
      include "../commons/Evs.h"
      include "../commons/Nf_FF.h"
      include "../commons/m2th.h"
      include "../commons/grid.h"
      include "../commons/EvolutionMatrices.h"
      include "../commons/wrap.h"
      include "../commons/MaxFlavourPDFs.h"
**
*     Input Variables
*
      double precision muF20,muF2
**
*     Internal Variables
*
      integer inf
      double precision mu2i(3:6),mu2f(4:7)
      double precision M0sg(3,3,0:nint_max,0:nint_max)
      double precision M0nsp(0:nint_max,0:nint_max)
      double precision M0nsm(0:nint_max,0:nint_max)
      double precision Msg(3,3,0:nint_max,0:nint_max)
      double precision Mnsp(0:nint_max,0:nint_max)
      double precision Mnsm(0:nint_max,0:nint_max)
*
*     Initial Conditions (Unity matrices)
*
      call IdentityOperatorsQED(M0sg,M0nsp,M0nsm)
*
*     Mass scheme
*
*     Fixed Flavour Number Scheme
*
      if(Evs.eq."FF")then
         nfi = Nf_FF
         nff = Nf_FF
         wnf = Nf_FF
*     If initial and final energies are equal return immediately the intial conditions
         if(muF2.eq.muF20)then
            call EqualOperatorsQEDnf(Nf_FF,M0sg,M0nsp,M0nsm,
     1                                     MQEDsg,MQEDnsp,MQEDnsm)
            return
         endif
*     Singlet
         call odeintsgQED(muF20,muF2,M0sg,Msg)
*     Non-Singlet
         call odeintnsQED(1,muF20,muF2,M0nsp,Mnsp)
         call odeintnsQED(2,muF20,muF2,M0nsm,Mnsm)
*
*     Put arrays into common arrays depending on the number of active flavours
*
         call EqualOperatorsQEDnf(Nf_FF,Msg,Mnsp,Mnsm,
     1                                  MQEDsg,MQEDnsp,MQEDnsm)
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
         if(nff.gt.nfMaxPDFs) nff = nfMaxPDFs
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
         if(nfi.gt.nfMaxPDFs) nfi = nfMaxPDFs
*     If initial and final energies are equal return immediately the intial conditions
         if(muF2.eq.muF20)then
            call EqualOperatorsQEDnf(nfi,M0sg,M0nsp,M0nsm,
     1                                   MQEDsg,MQEDnsp,MQEDnsm)
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
               mu2f(inf) = m2th(inf+1)
            enddo
         elseif(sgn.eq.-1)then
            do inf=nfi-1,nff,sgn
               mu2i(inf) = m2th(inf+1)
            enddo
            do inf=nfi,nff+1,sgn
               mu2f(inf) = m2th(inf)
            enddo
         endif
         mu2f(nff) = muF2
*
         do inf=nfi,nff
            wnf = inf
*     Singlet
            call odeintsgQED(mu2i(inf),mu2f(inf),M0sg,Msg)
*     Non-Singlet
            call odeintnsQED(1,mu2i(inf),mu2f(inf),M0nsp,Mnsp)
            call odeintnsQED(2,mu2i(inf),mu2f(inf),M0nsm,Mnsm)
*
*     Put arrays into common arrays depending on the number of active flavours
*
            call EqualOperatorsQEDnf(inf,Msg,Mnsp,Mnsm,
     1                                   MQEDsg,MQEDnsp,MQEDnsm)
         enddo
      endif
*
      return
      end
