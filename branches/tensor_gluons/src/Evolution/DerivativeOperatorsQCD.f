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
      include "../commons/consts.h"
**
*     Input Variables
*
      double precision muF2
**
*     Internal Variables
*
      integer inf
      double precision coup,alphasPDF!,a_QCD
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
c      coup = a_QCD(muF2)
      coup = alphasPDF(dsqrt(muF2)) / 4d0 / pi
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
