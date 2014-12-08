************************************************************************
*
*     DerivativeOperatorsQED.f:
*
*     This routine returns the singlet and the non-singlet derivative
*     operators on the x-space grid between at the scale muF2.
*
************************************************************************
      subroutine DerivativeOperatorsQED(muF2)
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
      double precision muF2
**
*     Internal Variables
*
      integer inf
      double precision coup,a_QED
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
         if(inf.gt.nfMaxPDFs) inf = nfMaxPDFs
      endif
*
      coup = a_QED(muF2)
*
*     Compute operators
*
      wnf = inf
*     Singlet
      call DeriveSgQED(coup,dMQEDsg)
*     Non-Singlet
      call DeriveNsQED(1,coup,dMQEDnsp)
      call DeriveNsQED(2,coup,dMQEDnsm)
*
      return
      end
