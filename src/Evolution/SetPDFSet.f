************************************************************************
*
*     SetPDFSet.f:
*
*     This subroutine sets the name of the PDF set to be used at the
*     initial scale.
*
************************************************************************
      subroutine SetPDFSet(name)
*
      implicit none
*
      include "../commons/pdfset.h"
*
*     Variables
*
      character*(*) name
*
      pdfset = name
      pdfsetlen = len(name)
      InPDFs = "done"
*     
      return
      end
