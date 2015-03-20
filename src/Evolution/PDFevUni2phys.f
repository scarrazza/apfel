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
      subroutine PDFevUni2phys(leptonin,pdfin,leptonout,pdfout)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/transUni.h"
**
*     Input Variables
*
      double precision pdfin(0:13,0:nint_max)
      double precision leptonin(6,0:nint_max)
**
*     Internal Variables
*
      integer a
      integer i,j
      double precision pdftmp(-7:6)
**
*     Output Variables
*
      double precision pdfout(-6:6,0:nint_max)
      double precision leptonout(-3:3,0:nint_max)
*
*     Rotate PDFs
*
      do a=0,nin(igrid)
         do i=0,13
            pdftmp(i-7) = 0d0
            do j=0,13
               pdftmp(i-7) = pdftmp(i-7)
     1                     + Tev2phUni(i,j) * pdfin(j,a)
            enddo
         enddo
*
         do i=-6,6
            pdfout(i,a) = pdftmp(i)
         enddo
*
         leptonout(-3,a) = ( ( leptonin(1,a) - leptonin(4,a) ) 
     1                   -   ( leptonin(3,a) - leptonin(6,a) ) ) / 6d0

         leptonout(-2,a) = ( 2d0 * ( leptonin(1,a) - leptonin(4,a) ) 
     1                   -   3d0 * ( leptonin(2,a) - leptonin(5,a) )
     2                   +         ( leptonin(3,a) - leptonin(6,a) ) )
     3                   / 12d0
         leptonout(-1,a) = ( 2d0 * ( leptonin(1,a) - leptonin(4,a) ) 
     1                   +   3d0 * ( leptonin(2,a) - leptonin(5,a) )
     2                   +         ( leptonin(3,a) - leptonin(6,a) ) ) 
     3                   / 12d0
         leptonout(0,a)  = pdftmp(-7)
         leptonout(1,a)  = ( 2d0 * ( leptonin(1,a) + leptonin(4,a) ) 
     1                   +   3d0 * ( leptonin(2,a) + leptonin(5,a) )
     2                   +         ( leptonin(3,a) + leptonin(6,a) ) )
     3                   / 12d0
         leptonout(2,a)  = ( 2d0 * ( leptonin(1,a) + leptonin(4,a) ) 
     1                   -   3d0 * ( leptonin(2,a) + leptonin(5,a) )
     2                   +         ( leptonin(3,a) + leptonin(6,a) ) )
     3                   / 12d0
         leptonout(3,a)  = ( ( leptonin(1,a) + leptonin(4,a) )
     1                   -   ( leptonin(3,a) + leptonin(6,a) ) ) / 6d0 
      enddo
*
      return
      end
