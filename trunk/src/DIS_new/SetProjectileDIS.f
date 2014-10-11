************************************************************************
*
*     SetProjectileDIS.f:
*
*     This subroutine sets the DIS process.
*
************************************************************************
      subroutine SetProjectileDIS(lept)
*
      implicit none
*
      include "../commons/ProjectileDIS.h"
*
*     Variables
*
      character*12 lept
*
      if(lept(1:8).eq."electron")then
         ProjectileDIS = lept(1:8)
      elseif(lept(1:8).eq."positron")then
         ProjectileDIS = lept(1:8)
      elseif(lept(1:8).eq."neutrino")then
         ProjectileDIS = lept(1:8)
      elseif(lept(1:8).eq."antineutrino")then
         ProjectileDIS = lept
      endif
      InProjectileDIS = "done"
*
      return
      end
