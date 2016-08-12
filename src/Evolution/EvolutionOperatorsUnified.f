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
      include "../commons/TauMass.h"
      include "../commons/grid.h"
      include "../commons/EvolutionMatrices.h"
      include "../commons/wrap.h"
      include "../commons/MaxFlavourPDFs.h"
      include "../commons/MaxFlavourAlpha.h"
      include "../commons/LeptEvol.h"
**
*     Input Variables
*
      double precision muF20,muF2
**
*     Internal Variables
*
      integer inf,inl,nfmax
      double precision mu2i(3:7),mu2f(3:7)
      double precision muF20l,muF2l
      double precision M0sg1(5,5,0:nint_max,0:nint_max)
      double precision M0sg2(2,2,0:nint_max,0:nint_max)
      double precision M0nspu(0:nint_max,0:nint_max)
      double precision M0nspd(0:nint_max,0:nint_max)
      double precision M0nsmu(0:nint_max,0:nint_max)
      double precision M0nsmd(0:nint_max,0:nint_max)
      double precision M0nslep(0:nint_max,0:nint_max)
      double precision Msg1(5,5,0:nint_max,0:nint_max)
      double precision Msg2(2,2,0:nint_max,0:nint_max)
      double precision Mnspu(0:nint_max,0:nint_max)
      double precision Mnspd(0:nint_max,0:nint_max)
      double precision Mnsmu(0:nint_max,0:nint_max)
      double precision Mnsmd(0:nint_max,0:nint_max)
      double precision Mnslep(0:nint_max,0:nint_max)
      double precision tiny
      parameter(tiny=1d-10)
*
*     Initial Conditions (Unity matrices)
*
      call IdentityOperatorsUnified(M0sg1,M0sg2,M0nspu,M0nspd,
     1                                          M0nsmu,M0nsmd,M0nslep)
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
*
         nli = 2
         nlf = 2
         nfli(nli) = Nf_FF
         nflf(nli) = Nf_FF
         wnf = Nf_FF
         wnl = nli
         sgn = 1
*     If initial and final energies are equal return immediately the intial conditions
         if(muF2.eq.muF20)then
            call EqualOperatorsUnifiednf(Nf_FF,nli,
     1           M0sg1,M0sg2,M0nspu,M0nspd,M0nsmu,M0nsmd,M0nslep,
     2           MUnisg1,Munisg2,MUninspu,MUninspd,MUninsmu,MUninsmd,
     3           MUninslep)
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
*     Non-Singlet lepton
         if(LeptEvol)then
            call odeintnsUnified(5,muF20,muF2,M0nslep,Mnslep)
         endif
*
         call EqualOperatorsUnifiednf(Nf_FF,nli,
     1        Msg1,Msg2,Mnspu,Mnspd,Mnsmu,Mnsmd,Mnslep,
     2        MUnisg1,MUnisg2,MUninspu,MUninspd,MUninsmu,MUninsmd,
     3        MUninslep)
*
*     Variable Flavour Number Scheme
*
      elseif(Evs.eq."VF")then
*     If initial and final energies are equal return immediately the intial conditions
         if(muF2.ge.muF20)then
            sgn = 1
         elseif(muF2.lt.muF20)then
            sgn = - 1
         endif
*
         nli = 2
         if(muF20.gt.MTau**2) nli = 3
         nlf = 2
         if(muF2.gt.MTau**2)  nlf = 3
*
         if(nli.eq.nlf)then
            muF20l = muF20
            muF2l  = muF2
         else
            muF20l = muF20
            muF2l  = MTau**2
         endif
*
         do inl=nli,nlf,sgn
            wnl = inl
*
*     Find initial and final number of flavours
*
            if(muF2l.gt.m2th(6))then
               nflf(inl) = 6
            elseif(muF2l.gt.m2th(5))then
               nflf(inl) = 5
            elseif(muF2l.gt.m2th(4))then
               nflf(inl) = 4
            else
               nflf(inl) = 3
            endif
            if(nflf(inl).gt.nfmax) nflf(inl) = nfmax
*
            if(muF20l.gt.m2th(6))then
               nfli(inl) = 6
            elseif(muF20l.gt.m2th(5))then
               nfli(inl) = 5
            elseif(muF20l.gt.m2th(4))then
               nfli(inl) = 4
            else
               nfli(inl) = 3
            endif
            if(nfli(inl).gt.nfmax) nfli(inl) = nfmax
*
            mu2i(nfli(inl)) = muF20l
            if(sgn.eq.1)then
               do inf=nfli(inl)+1,nflf(inl)
                  mu2i(inf) = m2th(inf)
               enddo
               do inf=nfli(inl),nflf(inl)-1
                  mu2f(inf) = m2th(inf+1) - tiny
               enddo
            elseif(sgn.eq.-1)then
               do inf=nfli(inl)-1,nflf(inl),sgn
                  mu2i(inf) = m2th(inf+1) + tiny
               enddo
               do inf=nfli(inl),nflf(inl)+1,sgn
                  mu2f(inf) = m2th(inf)
               enddo
            endif
            mu2f(nflf(inl)) = muF2l
*
            do inf=nfli(inl),nflf(inl),sgn
               if(muF2.eq.muF20)then
                  call EqualOperatorsUnifiednf(inf,inl,
     1                 M0sg1,M0sg2,M0nspu,M0nspd,M0nsmu,M0nsmd,M0nslep,
     2                 MUnisg1,Munisg2,
     3                 MUninspu,MUninspd,MUninsmu,MUninsmd,MUninslep)
                  return
               endif
*
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
*     Non-Singlet leptons
               if(LeptEvol)then
                  call odeintnsUnified(4,mu2i(inf),mu2f(inf),M0nslep,
     1                                                       Mnslep)
               endif
*     Match operators
               call EqualOperatorsUnifiednf(inf,inl,
     1              Msg1,Msg2,Mnspu,Mnspd,Mnsmu,Mnsmd,Mnslep,
     2              MUnisg1,Munisg2,MUninspu,MUninspd,MUninsmu,MUninsmd,
     3              MUninslep)
            enddo
            muF20l = MTau**2
            muF2l  = muF2
         enddo
      endif
*
      return
      end
