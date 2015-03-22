************************************************************************
*
*     PDFphys2evQED.f:
*
*     This routine converts PDFs from the physical basis to the QED
*     evolution basis:
*
*     Evolution basis:
*       1   2   3   4   5   6   7   8   9  10  11  12  13
*      gam  Sg  D   Duc Dds Dsb Dct u-  d-  s-  c-  b-  t-
*
*     Physical basis:
*      -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6 
*      tb  bb  cb  sb  ub  db  gam  d   u   s   c   b   t
*
************************************************************************
      subroutine PDFphys2evQED(pdfin,pdfout)
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
      integer i,j
      double precision trans(-6:6,-6:6)
**
*     Output Variables
*
      double precision pdfout(0:13,0:nint_max)
*
*     Define rotation matrix
*
      trans(-6,-6) = 0d0 
      trans(-6,-5) = 0d0 
      trans(-6,-4) = 0d0 
      trans(-6,-3) = 0d0 
      trans(-6,-2) = 0d0 
      trans(-6,-1) = 0d0 
      trans(-6,0)  = 1d0 
      trans(-6,1)  = 0d0 
      trans(-6,2)  = 0d0 
      trans(-6,3)  = 0d0 
      trans(-6,4)  = 0d0 
      trans(-6,5)  = 0d0 
      trans(-6,6)  = 0d0 
                         
      trans(-5,-6) = 1d0 
      trans(-5,-5) = 1d0 
      trans(-5,-4) = 1d0 
      trans(-5,-3) = 1d0 
      trans(-5,-2) = 1d0 
      trans(-5,-1) = 1d0 
      trans(-5,0)  = 0d0 
      trans(-5,1)  = 1d0 
      trans(-5,2)  = 1d0 
      trans(-5,3)  = 1d0 
      trans(-5,4)  = 1d0 
      trans(-5,5)  = 1d0 
      trans(-5,6)  = 1d0 
                         
      trans(-4,-6) = 1d0 
      trans(-4,-5) = -1d0
      trans(-4,-4) = 1d0 
      trans(-4,-3) = -1d0
      trans(-4,-2) = 1d0 
      trans(-4,-1) = -1d0
      trans(-4,0)  = 0d0 
      trans(-4,1)  = -1d0
      trans(-4,2)  = 1d0 
      trans(-4,3)  = -1d0
      trans(-4,4)  = 1d0 
      trans(-4,5)  = -1d0
      trans(-4,6)  = 1d0 
                         
      trans(-3,-6) = 0d0 
      trans(-3,-5) = 0d0 
      trans(-3,-4) = -1d0
      trans(-3,-3) = 0d0 
      trans(-3,-2) = 1d0 
      trans(-3,-1) = 0d0 
      trans(-3,0)  = 0d0 
      trans(-3,1)  = 0d0 
      trans(-3,2)  = 1d0 
      trans(-3,3)  = 0d0 
      trans(-3,4)  = -1d0
      trans(-3,5)  = 0d0 
      trans(-3,6)  = 0d0 
                         
      trans(-2,-6) = 0d0 
      trans(-2,-5) = 0d0 
      trans(-2,-4) = 0d0 
      trans(-2,-3) = -1d0
      trans(-2,-2) = 0d0 
      trans(-2,-1) = 1d0 
      trans(-2,0)  = 0d0 
      trans(-2,1)  = 1d0 
      trans(-2,2)  = 0d0 
      trans(-2,3)  = -1d0
      trans(-2,4)  = 0d0 
      trans(-2,5)  = 0d0 
      trans(-2,6)  = 0d0 
                         
      trans(-1,-6) = 0d0 
      trans(-1,-5) = -1d0
      trans(-1,-4) = 0d0 
      trans(-1,-3) = 1d0 
      trans(-1,-2) = 0d0 
      trans(-1,-1) = 0d0 
      trans(-1,0)  = 0d0 
      trans(-1,1)  = 0d0 
      trans(-1,2)  = 0d0 
      trans(-1,3)  = 1d0 
      trans(-1,4)  = 0d0 
      trans(-1,5)  = -1d0
      trans(-1,6)  = 0d0 
                         
      trans(0,-6)  = -1d0
      trans(0,-5)  = 0d0 
      trans(0,-4)  = 1d0 
      trans(0,-3)  = 0d0 
      trans(0,-2)  = 0d0 
      trans(0,-1)  = 0d0 
      trans(0,0)   = 0d0 
      trans(0,1)   = 0d0 
      trans(0,2)   = 0d0 
      trans(0,3)   = 0d0 
      trans(0,4)   = 1d0 
      trans(0,5)   = 0d0 
      trans(0,6)   = -1d0
                         
      trans(1,-6)  = 0d0 
      trans(1,-5)  = 0d0 
      trans(1,-4)  = 0d0 
      trans(1,-3)  = 0d0 
      trans(1,-2)  = -1d0
      trans(1,-1)  = 0d0 
      trans(1,0)   = 0d0 
      trans(1,1)   = 0d0 
      trans(1,2)   = 1d0 
      trans(1,3)   = 0d0 
      trans(1,4)   = 0d0 
      trans(1,5)   = 0d0 
      trans(1,6)   = 0d0 
                         
      trans(2,-6)  = 0d0 
      trans(2,-5)  = 0d0 
      trans(2,-4)  = 0d0 
      trans(2,-3)  = 0d0 
      trans(2,-2)  = 0d0 
      trans(2,-1)  = -1d0
      trans(2,0)   = 0d0 
      trans(2,1)   = 1d0 
      trans(2,2)   = 0d0 
      trans(2,3)   = 0d0 
      trans(2,4)   = 0d0 
      trans(2,5)   = 0d0 
      trans(2,6)   = 0d0 
                         
      trans(3,-6)  = 0d0 
      trans(3,-5)  = 0d0 
      trans(3,-4)  = 0d0 
      trans(3,-3)  = -1d0
      trans(3,-2)  = 0d0 
      trans(3,-1)  = 0d0 
      trans(3,0)   = 0d0 
      trans(3,1)   = 0d0 
      trans(3,2)   = 0d0 
      trans(3,3)   = 1d0 
      trans(3,4)   = 0d0 
      trans(3,5)   = 0d0 
      trans(3,6)   = 0d0 
                         
      trans(4,-6)  = 0d0 
      trans(4,-5)  = 0d0 
      trans(4,-4)  = -1d0
      trans(4,-3)  = 0d0 
      trans(4,-2)  = 0d0 
      trans(4,-1)  = 0d0 
      trans(4,0)   = 0d0 
      trans(4,1)   = 0d0 
      trans(4,2)   = 0d0 
      trans(4,3)   = 0d0 
      trans(4,4)   = 1d0 
      trans(4,5)   = 0d0 
      trans(4,6)   = 0d0 
                         
      trans(5,-6)  = 0d0 
      trans(5,-5)  = -1d0
      trans(5,-4)  = 0d0 
      trans(5,-3)  = 0d0 
      trans(5,-2)  = 0d0 
      trans(5,-1)  = 0d0 
      trans(5,0)   = 0d0 
      trans(5,1)   = 0d0 
      trans(5,2)   = 0d0 
      trans(5,3)   = 0d0 
      trans(5,4)   = 0d0 
      trans(5,5)   = 1d0 
      trans(5,6)   = 0d0 
                         
      trans(6,-6)  = -1d0
      trans(6,-5)  = 0d0 
      trans(6,-4)  = 0d0 
      trans(6,-3)  = 0d0 
      trans(6,-2)  = 0d0 
      trans(6,-1)  = 0d0 
      trans(6,0)   = 0d0 
      trans(6,1)   = 0d0 
      trans(6,2)   = 0d0 
      trans(6,3)   = 0d0 
      trans(6,4)   = 0d0 
      trans(6,5)   = 0d0 
      trans(6,6)   = 1d0 
*
*     Rotate PDFs
*
      do a=0,nin(igrid)
         pdfout(0,a) = 0d0
         do i=-6,6
            pdfout(i+7,a) = 0d0
            do j=-6,6
               pdfout(i+7,a) = pdfout(i+7,a) + trans(i,j) * pdfin(j,a)
            enddo
         enddo
      enddo
*
      return
      end
