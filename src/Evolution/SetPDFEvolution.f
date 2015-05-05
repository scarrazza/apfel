************************************************************************
*
*     SetPDFEvolution.f:
*
*     This subroutine sets the solution of the coupling equations.
*
************************************************************************
      subroutine SetPDFEvolution(pe)
*
      implicit none
*
      include "../commons/PDFEvolution.h"
*
*     Variables
*
      character*11 pe
*
      if(pe(1:7).eq."exactmu")then
         PDFEvol = pe(1:7)
      elseif(pe(1:10).eq."exactalpha")then
         PDFEvol = pe(1:10)
      elseif(pe(1:11).eq."expandalpha")then
         PDFEvol = pe(1:11)
      elseif(pe(1:9).eq."truncated")then
         PDFEvol = pe(1:9)
      endif
      InPDFEvol = "done"
*
      return
      end
