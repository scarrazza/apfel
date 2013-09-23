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
      integer ln
      character*100 name
*
      ln = index(name,char(0)) - 1
      pdfset = name(1:ln)
      InPDFs = "done"
*
      return
      end
