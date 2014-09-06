************************************************************************
*
*     PDFevUni2phys.f:
*
*     This routine converts PDFs from the unified evolution basis to the 
*     physical basis:
*
*     Evolution basis:
*     0   1   2   3   4   5   6   7   8   9   10  11  12  13
*     g   gm  Sig Dsg Tu1 Tu2 Td1 Td2 V   DV  Vu1 Vu2 Vd1 Vd2
*
*     Physical basis:
*     -7  -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6 
*     gm  tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
*
************************************************************************
      subroutine PDFevUni2phys(pdfin,pdfout)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/transUni.h"
**
*     Input Variables
*
      double precision pdfin(0:13,0:nint_max)
**
*     Internal Variables
*
      integer a
      integer i,j
**
*     Output Variables
*
      double precision pdfout(-7:6,0:nint_max)
*
*     Rotate PDFs
*
      do a=0,nin(igrid)
         do i=0,13
            pdfout(i-7,a) = 0d0
            do j=0,13
               pdfout(i-7,a) = pdfout(i-7,a)
     1                       + Tev2phUni(i,j) * pdfin(j,a)
            enddo
         enddo
      enddo
*
      return
      end
