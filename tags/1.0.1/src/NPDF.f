************************************************************************
*
*     NPDF.f:
*
*     This function returns the N-th Mellin moment of the i-th PDF 
*     in the physical basis at the final scale.
*
************************************************************************
      function NPDF(i,N)
*
      implicit none
*
      include "../commons/grid.h"
**
*     Input Variables
*
      integer i,N
**
*     Input Variables
*
      integer j,M
      double precision xPDFwrap
      double precision dgauss,a,b,eps
      external xPDFwrap

      common / PDFindex / j,M
**
*     Output Variables
*
      double precision NPDF
*
      j   = i
      M   = N
*
      a   = xmin(1)
      b   = xmax
      eps = 1d-7
      NPDF = dgauss(xPDFwrap,a,b,eps)
*
      return
      end
************************************************************************
*
*     Wrapping of the function xPDF
*
************************************************************************
      function xPDFwrap(x)
*
      implicit none
**
*     Input variables
*
      double precision x
**
*     Internal variables
*
      integer j,M
      double precision xPDF

      common / PDFindex / j,M
**
*     Output variables
*
      double precision xPDFwrap
*
      xPDFwrap = x**(M-2) * xPDF(j,x)
*
      return
      end
