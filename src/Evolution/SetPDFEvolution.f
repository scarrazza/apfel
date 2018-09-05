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
      character*(*) pe
*
      PDFEvol = trim(pe)
      InPDFEvol = "done"
*
      return
      end
