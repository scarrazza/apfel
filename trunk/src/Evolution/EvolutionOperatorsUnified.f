************************************************************************
*
*     EvolutionOperatorsUnified.f:
*
*     This routine returns the singlet and the non-singlet evolution
*     operators on the x-space grid between the scales m20 and mu2 for
*     nf active flavours given the initial evolution operators in the
*     unified evolution basis.
*
************************************************************************
      subroutine EvolutionOperatorsUnified(muF20,muF2)
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
**
*     Input Variables
*
      double precision muF20,muF2
**
*     Internal Variables
*
      integer inf
      double precision mu2i(3:6),mu2f(4:7)
      double precision M0sg1(4,4,0:nint_max,0:nint_max)
      double precision M0sg2(2,2,0:nint_max,0:nint_max)
      double precision M0nspu(0:nint_max,0:nint_max)
      double precision M0nspd(0:nint_max,0:nint_max)
      double precision M0nsmu(0:nint_max,0:nint_max)
      double precision M0nsmd(0:nint_max,0:nint_max)
      double precision Msg1(4,4,0:nint_max,0:nint_max)
      double precision Msg2(2,2,0:nint_max,0:nint_max)
      double precision Mnspu(0:nint_max,0:nint_max)
      double precision Mnspd(0:nint_max,0:nint_max)
      double precision Mnsmu(0:nint_max,0:nint_max)
      double precision Mnsmd(0:nint_max,0:nint_max)
      double precision tiny
      parameter(tiny=1d-10)

      integer alpha,beta,i,j
*
*     Initial Conditions (Unity matrices)
*
      call IdentityOperatorsUnified(M0sg1,M0sg2,M0nspu,M0nspd,
     1                                          M0nsmu,M0nsmd)
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
            call EqualOperatorsUnifiednf(Nf_FF,
     1           M0sg1,M0sg2,M0nspu,M0nspd,M0nsmu,M0nsmd,
     2           MUnisg1,Munisg2,MUninspu,MUninspd,MUninsmu,MUninsmd)
            return
         endif
*     Singlet 1
         call odeintsgUnifiedS1(muF20,muF2,M0sg1,Msg1)
*     Singlet 2
         call odeintsgUnifiedS2(muF20,muF2,M0sg2,Msg2)
*     Non-Singlet plus-up
         call odeintnsUnified(1,muF20,muF2,M0nspu,Mnspu)
*     Non-Singlet plus-down
         call odeintnsUnified(2,muF20,muF2,M0nspd,Mnspd)
*     Non-Singlet minus-up
         call odeintnsUnified(3,muF20,muF2,M0nsmu,Mnsmu)
*     Non-Singlet minus-down
         call odeintnsUnified(4,muF20,muF2,M0nsmd,Mnsmd)
*
         call EqualOperatorsUnifiednf(Nf_FF,
     1        Msg1,Msg2,Mnspu,Mnspd,Mnsmu,Mnsmd,
     2        MUnisg1,MUnisg2,MUninspu,MUninspd,MUninsmu,MUninsmd)
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
            call EqualOperatorsUnifiednf(nfi,
     1           M0sg1,M0sg2,M0nspu,M0nspd,M0nsmu,M0nsmd,
     2           MUnisg1,Munisg2,MUninspu,MUninspd,MUninsmu,MUninsmd)
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
*     Singlet 1
            call odeintsgUnifiedS1(mu2i(inf),mu2f(inf),M0sg1,Msg1)
*     Singlet 2
            call odeintsgUnifiedS2(mu2i(inf),mu2f(inf),M0sg2,Msg2)
*     Non-Singlet plus-up
            call odeintnsUnified(1,mu2i(inf),mu2f(inf),M0nspu,Mnspu)
*     Non-Singlet plus-down
            call odeintnsUnified(2,mu2i(inf),mu2f(inf),M0nspd,Mnspd)
*     Non-Singlet minus-up
            call odeintnsUnified(3,mu2i(inf),mu2f(inf),M0nsmu,Mnsmu)
*     Non-Singlet minus-down
            call odeintnsUnified(4,mu2i(inf),mu2f(inf),M0nsmd,Mnsmd)
*     Match operators
            call EqualOperatorsUnifiednf(inf,
     1           Msg1,Msg2,Mnspu,Mnspd,Mnsmu,Mnsmd,
     2           MUnisg1,Munisg2,MUninspu,MUninspd,MUninsmu,MUninsmd)
         enddo
      endif
*
      return
      end
