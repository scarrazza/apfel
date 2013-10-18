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
c      ln = index(name,char(0)) - 1
      if(name.eq."ToyLH")then
         ln = 5
      elseif(name.eq."private")then
         ln = 7
      else
         ln = index(name,"LHgrid") + 5
      endif
      pdfset = name(1:ln)
      InPDFs = "done"
*
      return
      end
