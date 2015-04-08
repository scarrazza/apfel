************************************************************************
*
*     PDFphys2ev.f:
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
      subroutine PDFphys2ev(pdfin,pdfout)
*
      implicit none
**
*     Input Variables
*
      double precision pdfin(-6:6)
**
*     Internal Variables
*
      integer ipdf
**
*     Output Variables
*
      double precision pdfout(13)
*
      pdfout(1) = 0d0
      do ipdf=1,6
         pdfout(1) = pdfout(1) + ( pdfin(ipdf) + pdfin(-ipdf) )
      enddo
*
      pdfout(2)  = pdfin(0)
*
      pdfout(3) = 0d0
      do ipdf=1,6
         pdfout(3) = pdfout(3) + ( pdfin(ipdf) - pdfin(-ipdf) )
      enddo
*
      pdfout(4)  = ( pdfin(2) - pdfin(-2) ) 
     1           - ( pdfin(1) - pdfin(-1) )
      pdfout(5)  = ( pdfin(2) - pdfin(-2) )
     1           + ( pdfin(1) - pdfin(-1) )
     2           - 2d0 * ( pdfin(3) - pdfin(-3) )
      pdfout(6)  = ( pdfin(2) - pdfin(-2) )
     1           + ( pdfin(1) - pdfin(-1) )
     2           + ( pdfin(3) - pdfin(-3) ) 
     3           - 3d0 * ( pdfin(4) - pdfin(-4) )
      pdfout(7)  = ( pdfin(2) - pdfin(-2) )
     1           + ( pdfin(1) - pdfin(-1) )
     2           + ( pdfin(3) - pdfin(-3) )
     3           + ( pdfin(4) - pdfin(-4) )
     4           - 4d0 * ( pdfin(5) - pdfin(-5) )
      pdfout(8)  = ( pdfin(2) - pdfin(-2) )
     1           + ( pdfin(1) - pdfin(-1) )
     2           + ( pdfin(3) - pdfin(-3) )
     3           + ( pdfin(4) - pdfin(-4) )
     4           + ( pdfin(5) - pdfin(-5) ) 
     5           - 5d0 * ( pdfin(6) - pdfin(-6) )
      pdfout(9)  = ( pdfin(2) + pdfin(-2) )
     1           - ( pdfin(1) + pdfin(-1) )
      pdfout(10) = ( pdfin(2) + pdfin(-2) )
     1           + ( pdfin(1) + pdfin(-1) )
     2           - 2d0 * ( pdfin(3) + pdfin(-3) )
      pdfout(11) = ( pdfin(2) + pdfin(-2) )
     1           + ( pdfin(1) + pdfin(-1) )
     2           + ( pdfin(3) + pdfin(-3) ) 
     3           - 3d0 * ( pdfin(4) + pdfin(-4) )
      pdfout(12) = ( pdfin(2) + pdfin(-2) )
     1           + ( pdfin(1) + pdfin(-1) )
     2           + ( pdfin(3) + pdfin(-3) )
     3           + ( pdfin(4) + pdfin(-4) )
     4           - 4d0 * ( pdfin(5) + pdfin(-5) )
      pdfout(13) = ( pdfin(2) + pdfin(-2) )
     1           + ( pdfin(1) + pdfin(-1) )
     2           + ( pdfin(3) + pdfin(-3) )
     3           + ( pdfin(4) + pdfin(-4) )
     4           + ( pdfin(5) + pdfin(-5) ) 
     5           - 5d0 * ( pdfin(6) + pdfin(-6) )
*
      return
      end
