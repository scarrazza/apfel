************************************************************************
*
*     SelectCharge.f:
*
*     This subroutine selects one particular charge in the NC structure
*     functions (needed fo selecting up down and strange components).
*
************************************************************************
      subroutine SelectCharge(selch)
*
      implicit none
*
      include "../commons/SelectedCharge.h"
*
*     Variables
*
      character*7 selch
*
      if(selch(1:4).eq."down")then
         SelectedCharge = selch(1:4)
      elseif(selch(1:2).eq."up")then
         SelectedCharge = selch(1:2)
      elseif(selch(1:7).eq."strange")then
         SelectedCharge = selch(1:7)
      elseif(selch(1:5).eq."charm")then
         SelectedCharge = selch(1:5)
      elseif(selch(1:6).eq."bottom")then
         SelectedCharge = selch(1:6)
      elseif(selch(1:3).eq."top")then
         SelectedCharge = selch(1:3)
      elseif(selch(1:3).eq."all")then
         SelectedCharge = selch(1:3)
      endif
      InSelectedCharge = "done"
*
      return
      end
