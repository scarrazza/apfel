************************************************************************
*
*     SetMassScheme.f:
*
*     This subroutine sets the mass sche to be used for the computation
*     of the DIS observables.
*
************************************************************************
      subroutine SetMassScheme(ms)
*
      implicit none
*
      include "../commons/MassScheme.h"
*
*     Variables
*
      character*7 ms
*
      if(ms(1:4).eq."FFNS")then
         MassScheme = ms(1:4)
      else
         MassScheme = ms
      endif
      InMassScheme = "done"
*
      return
      end
