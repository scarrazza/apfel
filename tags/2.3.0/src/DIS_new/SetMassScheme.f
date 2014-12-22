************************************************************************
*
*     SetMassScheme.f:
*
*     This subroutine sets the mass scheme to be used for the computation
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
         if(ms(5:5).eq."3".or.
     1      ms(5:5).eq."4".or.
     2      ms(5:5).eq."5".or.
     3      ms(5:5).eq."6")then
            MassScheme = ms(1:5)
         else
            MassScheme = ms(1:4)
         endif
      elseif(ms(1:4).eq."FFN0")then
         if(ms(5:5).eq."3".or.
     1      ms(5:5).eq."4".or.
     2      ms(5:5).eq."5".or.
     3      ms(5:5).eq."6")then
            MassScheme = ms(1:5)
         else
            MassScheme = ms(1:4)
         endif
      else
         MassScheme = ms
      endif
      InMassScheme = "done"
*
      return
      end
