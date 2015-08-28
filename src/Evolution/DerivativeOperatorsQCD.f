************************************************************************
*
*     DerivativeOperatorsQCD.f:
*
*     This routine returns the singlet and the non-singlet derivative
*     operators on the x-space grid between at the scale muF2.
*
************************************************************************
      subroutine DerivativeOperatorsQCD(muF2)
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
      include "../commons/MaxFlavourAlpha.h"
**
*     Input Variables
*
      double precision muF2
**
*     Internal Variables
*
      integer inf,nfmax
      double precision coup,a_QCD
*
*     Define maximun number of flavours
*
      nfmax = max(nfMaxPDFs,nfMaxAlpha)
*
*     Mass scheme
*
*     Fixed Flavour Number Scheme
      if(Evs.eq."FF")then
         inf = Nf_FF
*     Variable Flavour Number Scheme
      elseif(Evs.eq."VF")then
*     Find number of flavours
         if(muF2.gt.m2th(6))then
            inf = 6
         elseif(muF2.gt.m2th(5))then
            inf = 5
         elseif(muF2.gt.m2th(4))then
            inf = 4
         else
            inf = 3
         endif
         if(inf.gt.nfmax) inf = nfmax
      endif
*
      coup = a_QCD(muF2)
*
*     Compute operators
*
      wnf = inf
*     Singlet
      call DeriveSgQCD(coup,dMQCDsg)
*     Non-Singlet
      call DeriveNsQCD(1,coup,dMQCDnsp)
      call DeriveNsQCD(2,coup,dMQCDnsm)
      call DeriveNsQCD(3,coup,dMQCDnsv)
*
      return
      end
