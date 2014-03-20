************************************************************************
*
*     SetAlphaEvolution.f:
*
*     This subroutine sets the solution of the coupling equations.
*
************************************************************************
      subroutine SetAlphaEvolution(ae)
*
      implicit none
*
      include "../commons/AlphaEvolution.h"
*
*     Variables
*
      character*8 ae
*
      if(ae(1:5).eq."exact")then
         AlphaEvol = ae(1:5)
      elseif(ae(1:8).eq."expanded")then
         AlphaEvol = ae(1:8)
      elseif(ae(1:6).eq."lambda")then
         AlphaEvol = ae(1:6)
      endif
      InAlphaEvol = "done"
*
      return
      end
