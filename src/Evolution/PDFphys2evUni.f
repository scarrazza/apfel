************************************************************************
*
*     PDFphys2evUni.f:
*
*     This routine converts PDFs from the physical basis to the unified
*     evolution basis:
*
*     Physical basis quark:
*     -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6 
*     tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
*     Physical basis leptons:
*     -3    -2    -1     0     1     2     3
*     tau+  mu+   e+     gamma e-    mu-   tau-
*
*     Evolution basis quark:
*     0   1   2   3   4   5   6   7   8   9   10  11  12  13
*     g   gm  Sig Dsg Tu1 Tu2 Td1 Td2 V   DV  Vu1 Vu2 Vd1 Vd2
*     Evolution basis leptons:
*     1    2     3     4     5     6
*     Sgl  Tl3   Tl8   Vl    Vl3   Vl8
*
************************************************************************
      subroutine PDFphys2evUni(leptonin,pdfin,leptonout,pdfout)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/transUni.h"
**
*     Input Variables
*
      double precision leptonin(-3:3,0:nint_max)
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
      double precision leptonout(6,0:nint_max)
      double precision pdfout(0:13,0:nint_max)
*
*     Rotate PDFs
*
      do a=0,nin(igrid)
         pdfi(-7) = leptonin(0,a)
         do i=-6,6
            pdfi(i) = pdfin(i,a)
         enddo
         do i=0,13
            pdfout(i,a) = 0d0
            do j=0,13
               pdfout(i,a) = pdfout(i,a) + Tph2evUni(i,j) * pdfi(j-7)
            enddo
         enddo
         leptonout(1,a) = ( leptonin(1,a) + leptonin(-1,a) )
     1                  + ( leptonin(2,a) + leptonin(-2,a) )
     2                  + ( leptonin(3,a) + leptonin(-3,a) )
         leptonout(2,a) = ( leptonin(1,a) + leptonin(-1,a) )
     1                  - ( leptonin(2,a) + leptonin(-2,a) )
         leptonout(3,a) = ( leptonin(1,a) + leptonin(-1,a) )
     1                  + ( leptonin(2,a) + leptonin(-2,a) )
     2            - 2d0 * ( leptonin(3,a) + leptonin(-3,a) )
         leptonout(4,a) = ( leptonin(1,a) - leptonin(-1,a) )
     1                  + ( leptonin(2,a) - leptonin(-2,a) )
     2                  + ( leptonin(3,a) - leptonin(-3,a) )
         leptonout(5,a) = ( leptonin(1,a) - leptonin(-1,a) )
     1                  - ( leptonin(2,a) - leptonin(-2,a) )
         leptonout(6,a) = ( leptonin(1,a) - leptonin(-1,a) )
     1                  + ( leptonin(2,a) - leptonin(-2,a) )
     2            - 2d0 * ( leptonin(3,a) - leptonin(-3,a) )
      enddo
*
      return
      end
