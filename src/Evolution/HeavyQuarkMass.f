************************************************************************
*
*     HeavyQuarkMass.f:
*
*     This function returns the value of i-th (i=4,5,6) heavy quark mass
*     at the given scale using the parameters of the evolution.
*     Be careful because for the MSbar masses, if kren.ne.1, the matching
*     is not done at the heavy quark thresholds.
*
************************************************************************
      function HeavyQuarkMass(i,Q)
*
      implicit none
*
      include "../commons/consts.h"
      include "../commons/kren.h"
      include "../commons/m2th.h"
      include "../commons/mass_scheme.h"
      include "../commons/MassRunning.h"
**
*     Input Variables
*
      integer i
      double precision Q
**
*     Internal Variables
*
      double precision Q2
      double precision MSbarmass
**
*     Output Variables
*
      double precision HeavyQuarkMass
*
*     Check the consistency of the input
*
      if(i.lt.4.and.i.gt.6)then
         write(6,*) "In HeavyQuarkMass.f:"
         write(6,*) "Invalid heavy quark index i =",i
         call exit(-10)
      endif
*
      if(mass_scheme.eq."MSbar")then
         if(MassRunning)then
            Q2 = Q * Q / kren
            HeavyQuarkMass = MSbarmass(i,Q2)
         else
            HeavyQuarkMass = dsqrt(m2ph(i))
         endif
      elseif(mass_scheme.eq."Pole")then
         HeavyQuarkMass = dsqrt(m2ph(i))
      endif
*
      return
      end
