************************************************************************
*
*     PDFphys2evQCD.f:
*
*     This routine converts PDFs from the physical basis to the QCD
*     evolution basis:
*
*     Physical basis:
*      -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6 
*      tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
*
*     Evolution basis:
*       1   2   3   4   5   6   7   8   9  10  11  12  13
*      Sg   g   V  V3  V8 V15 V24 V35  T3  T8 T15 T24 T35
*
************************************************************************
      subroutine PDFphys2evQCD(pdfin,pdfout)
*
      implicit none
*
      include "../commons/grid.h"
**
*     Input Variables
*
      double precision pdfin(-6:6,0:nint_max)
**
*     Internal Variables
*
      integer a
      integer ipdf
**
*     Output Variables
*
      double precision pdfout(0:13,0:nint_max)
*
      do a=0,nin(igrid)
         pdfout(1,a) = 0d0
         do ipdf=1,6
            pdfout(1,a) = pdfout(1,a) + (pdfin(ipdf,a) + pdfin(-ipdf,a))
         enddo
*
         pdfout(2,a)  = pdfin(0,a)
*
         pdfout(3,a) = 0d0
         do ipdf=1,6
            pdfout(3,a) = pdfout(3,a) + (pdfin(ipdf,a) - pdfin(-ipdf,a))
         enddo
*
         pdfout(4,a)  = ( pdfin(2,a) - pdfin(-2,a) ) 
     1                - ( pdfin(1,a) - pdfin(-1,a) )
         pdfout(5,a)  = ( pdfin(2,a) - pdfin(-2,a) )
     1                + ( pdfin(1,a) - pdfin(-1,a) )
     2                - 2d0 * ( pdfin(3,a) - pdfin(-3,a) )
         pdfout(6,a)  = ( pdfin(2,a) - pdfin(-2,a) )
     1                + ( pdfin(1,a) - pdfin(-1,a) )
     2                + ( pdfin(3,a) - pdfin(-3,a) ) 
     3                - 3d0 * ( pdfin(4,a) - pdfin(-4,a) )
         pdfout(7,a)  = ( pdfin(2,a) - pdfin(-2,a) )
     1                + ( pdfin(1,a) - pdfin(-1,a) )
     2                + ( pdfin(3,a) - pdfin(-3,a) )
     3                + ( pdfin(4,a) - pdfin(-4,a) )
     4                - 4d0 * ( pdfin(5,a) - pdfin(-5,a) )
         pdfout(8,a)  = ( pdfin(2,a) - pdfin(-2,a) )
     1                + ( pdfin(1,a) - pdfin(-1,a) )
     2                + ( pdfin(3,a) - pdfin(-3,a) )
     3                + ( pdfin(4,a) - pdfin(-4,a) )
     4                + ( pdfin(5,a) - pdfin(-5,a) ) 
     5                - 5d0 * ( pdfin(6,a) - pdfin(-6,a) )
         pdfout(9,a)  = ( pdfin(2,a) + pdfin(-2,a) )
     1                - ( pdfin(1,a) + pdfin(-1,a) )
         pdfout(10,a) = ( pdfin(2,a) + pdfin(-2,a) )
     1                + ( pdfin(1,a) + pdfin(-1,a) )
     2                - 2d0 * ( pdfin(3,a) + pdfin(-3,a) )
         pdfout(11,a) = ( pdfin(2,a) + pdfin(-2,a) )
     1                + ( pdfin(1,a) + pdfin(-1,a) )
     2                + ( pdfin(3,a) + pdfin(-3,a) ) 
     3                - 3d0 * ( pdfin(4,a) + pdfin(-4,a) )
         pdfout(12,a) = ( pdfin(2,a) + pdfin(-2,a) )
     1                + ( pdfin(1,a) + pdfin(-1,a) )
     2                + ( pdfin(3,a) + pdfin(-3,a) )
     3                + ( pdfin(4,a) + pdfin(-4,a) )
     4                - 4d0 * ( pdfin(5,a) + pdfin(-5,a) )
         pdfout(13,a) = ( pdfin(2,a) + pdfin(-2,a) )
     1                + ( pdfin(1,a) + pdfin(-1,a) )
     2                + ( pdfin(3,a) + pdfin(-3,a) )
     3                + ( pdfin(4,a) + pdfin(-4,a) )
     4                + ( pdfin(5,a) + pdfin(-5,a) ) 
     5                - 5d0 * ( pdfin(6,a) + pdfin(-6,a) )
      enddo
*
      return
      end
