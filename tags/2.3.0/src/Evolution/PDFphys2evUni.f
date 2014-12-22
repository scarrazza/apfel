************************************************************************
*
*     PDFphys2evUni.f:
*
*     This routine converts PDFs from the physical basis to the unified
*     evolution basis:
*
*     Physical basis:
*     -7  -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6 
*     gm  tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
*
*     Evolution basis:
*     0   1   2   3   4   5   6   7   8   9   10  11  12  13
*     g   gm  Sig Dsg Tu1 Tu2 Td1 Td2 V   DV  Vu1 Vu2 Vd1 Vd2
*
************************************************************************
      subroutine PDFphys2evUni(gammain,pdfin,pdfout)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/transUni.h"
**
*     Input Variables
*
      double precision gammain(0:nint_max)
      double precision pdfin(-6:6,0:nint_max)
**
*     Internal Variables
*
      integer a
      integer i,j
      double precision pdfi(-7:6)
**
*     Output Variables
*
      double precision pdfout(0:13,0:nint_max)
*
*     Rotate PDFs
*
      do a=0,nin(igrid)
         pdfi(-7) = gammain(a)
         do i=-6,6
            pdfi(i) = pdfin(i,a)
         enddo
         do i=0,13
            pdfout(i,a) = 0d0
            do j=0,13
               pdfout(i,a) = pdfout(i,a) + Tph2evUni(i,j) * pdfi(j-7)
            enddo
         enddo
      enddo
*
      return
      end
