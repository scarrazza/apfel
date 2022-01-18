************************************************************************
*
*     SetRenFacRatioPDF.f:
*
*     This subroutine sets the ratio between renormalization and
*     factorization scales in the PDF evolution.
*
************************************************************************
      subroutine SetRenFacRatioPDF(ratio)
*
      implicit none
*
      include "../commons/krenpdf.h"
*
*     Variables
*
      double precision ratio
*
      kren      = ratio * ratio
      InKrenPDF = "done"
*
      return
      end
