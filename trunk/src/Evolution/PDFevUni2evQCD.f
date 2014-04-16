************************************************************************
*
*     PDFevUni2evQCD.f:
*
*     This routine converts PDFs from the unified evolution basis to the 
*     QCD evolution basis:
*
*     Unified evolution basis:
*     0   1   2   3   4   5   6   7   8   9   10  11  12  13
*     g   gm  Sig Dsg Tu1 Tu2 Td1 Td2 V   DV  Vu1 Vu2 Vd1 Vd2
*
*     QCD Evolution basis:
*     0   1   2   3   4   5   6   7   8   9  10  11  12  13
*     gm  Sg   g   V  V3  V8 V15 V24 V35  T3  T8 T15 T24 T35
*
************************************************************************
      subroutine PDFevUni2evQCD(pdfin,pdfout)
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
      double precision pdfout(0:13,0:nint_max)
*
*     Rotate PDFs
*
      do a=0,nin(igrid)
         do i=0,13
            pdfout(i,a) = 0d0
            do j=0,13
               pdfout(i,a) = pdfout(i,a)
     1                     + TevUni2evQCD(i,j) * pdfin(j,a)
            enddo
         enddo
      enddo
*
      return
      end
