************************************************************************
*
*     PDFevQED2phys.f:
*
*     This routine converts PDFs from the QDC evolution basis to the 
*     physical basis:
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
      subroutine PDFevQED2phys(pdfin,pdfout)
*
      implicit none
*
      include "../commons/grid.h"
**
*     Input Variables
*
      double precision pdfin(0:13,0:nint_max)
**
*     Internal Variables
*
      integer a
      integer i,j
      double precision trans(13,13)
**
*     Output Variables
*
      double precision pdfout(-6:6,0:nint_max)
*
*     Define rotation matrix
*
      trans(1,1)  = 0d0
      trans(1,2)  = 1d0
      trans(1,3)  = 1d0
      trans(1,4)  = -2d0
      trans(1,5)  = 0d0
      trans(1,6)  = 0d0
      trans(1,7)  = -4d0
      trans(1,8)  = 0d0
      trans(1,9)  = 0d0
      trans(1,10) = 0d0
      trans(1,11) = 0d0
      trans(1,12) = 0d0
      trans(1,13) = -6d0

      trans(2,1)  = 0d0
      trans(2,2)  = 1d0
      trans(2,3)  = -1d0
      trans(2,4)  = 0d0
      trans(2,5)  = -2d0
      trans(2,6)  = -4d0
      trans(2,7)  = 0d0
      trans(2,8)  = 0d0
      trans(2,9)  = 0d0
      trans(2,10) = 0d0
      trans(2,11) = 0d0
      trans(2,12) = -6d0
      trans(2,13) = 0d0

      trans(3,1)  = 0d0
      trans(3,2)  = 1d0
      trans(3,3)  = 1d0
      trans(3,4)  = -2d0
      trans(3,5)  = 0d0
      trans(3,6)  = 0d0
      trans(3,7)  = 2d0
      trans(3,8)  = 0d0
      trans(3,9)  = 0d0
      trans(3,10) = 0d0
      trans(3,11) = -6d0
      trans(3,12) = 0d0
      trans(3,13) = 0d0

      trans(4,1)  = 0d0
      trans(4,2)  = 1d0
      trans(4,3)  = -1d0
      trans(4,4)  = 0d0
      trans(4,5)  = -2d0
      trans(4,6)  = 2d0
      trans(4,7)  = 0d0
      trans(4,8)  = 0d0
      trans(4,9)  = 0d0
      trans(4,10) = -6d0
      trans(4,11) = 0d0
      trans(4,12) = 0d0
      trans(4,13) = 0d0

      trans(5,1)  = 0d0
      trans(5,2)  = 1d0
      trans(5,3)  = 1d0
      trans(5,4)  = 4d0
      trans(5,5)  = 0d0
      trans(5,6)  = 0d0
      trans(5,7)  = 2d0
      trans(5,8)  = -6d0
      trans(5,9)  = 0d0
      trans(5,10) = 0d0
      trans(5,11) = 0d0
      trans(5,12) = 0d0
      trans(5,13) = 0d0

      trans(6,1)  = 0d0
      trans(6,2)  = 1d0
      trans(6,3)  = -1d0
      trans(6,4)  = 0d0
      trans(6,5)  = 4d0
      trans(6,6)  = 2d0
      trans(6,7)  = 0d0
      trans(6,8)  = 0d0
      trans(6,9)  = -6d0
      trans(6,10) = 0d0
      trans(6,11) = 0d0
      trans(6,12) = 0d0
      trans(6,13) = 0d0

      trans(7,1)  = 12d0
      trans(7,2)  = 0d0
      trans(7,3)  = 0d0
      trans(7,4)  = 0d0
      trans(7,5)  = 0d0
      trans(7,6)  = 0d0
      trans(7,7)  = 0d0
      trans(7,8)  = 0d0
      trans(7,9)  = 0d0
      trans(7,10) = 0d0
      trans(7,11) = 0d0
      trans(7,12) = 0d0
      trans(7,13) = 0d0

      trans(8,1)  = 0d0
      trans(8,2)  = 1d0
      trans(8,3)  = -1d0
      trans(8,4)  = 0d0
      trans(8,5)  = 4d0
      trans(8,6)  = 2d0
      trans(8,7)  = 0d0
      trans(8,8)  = 0d0
      trans(8,9)  = 6d0
      trans(8,10) = 0d0
      trans(8,11) = 0d0
      trans(8,12) = 0d0
      trans(8,13) = 0d0

      trans(9,1)  = 0d0
      trans(9,2)  = 1d0
      trans(9,3)  = 1d0
      trans(9,4)  = 4d0
      trans(9,5)  = 0d0
      trans(9,6)  = 0d0
      trans(9,7)  = 2d0
      trans(9,8)  = 6d0
      trans(9,9)  = 0d0
      trans(9,10) = 0d0
      trans(9,11) = 0d0
      trans(9,12) = 0d0
      trans(9,13) = 0d0

      trans(10,1)  = 0d0
      trans(10,2)  = 1d0
      trans(10,3)  = -1d0
      trans(10,4)  = 0d0
      trans(10,5)  = -2d0
      trans(10,6)  = 2d0
      trans(10,7)  = 0d0
      trans(10,8)  = 0d0
      trans(10,9)  = 0d0
      trans(10,10) = 6d0
      trans(10,11) = 0d0
      trans(10,12) = 0d0
      trans(10,13) = 0d0

      trans(11,1)  = 0d0
      trans(11,2)  = 1d0
      trans(11,3)  = 1d0
      trans(11,4)  = -2d0
      trans(11,5)  = 0d0
      trans(11,6)  = 0d0
      trans(11,7)  = 2d0
      trans(11,8)  = 0d0
      trans(11,9)  = 0d0
      trans(11,10) = 0d0
      trans(11,11) = 6d0
      trans(11,12) = 0d0
      trans(11,13) = 0d0

      trans(12,1)  = 0d0
      trans(12,2)  = 1d0
      trans(12,3)  = -1d0
      trans(12,4)  = 0d0
      trans(12,5)  = -2d0
      trans(12,6)  = -4d0
      trans(12,7)  = 0d0
      trans(12,8)  = 0d0
      trans(12,9)  = 0d0
      trans(12,10) = 0d0
      trans(12,11) = 0d0
      trans(12,12) = 6d0
      trans(12,13) = 0d0

      trans(13,1)  = 0d0
      trans(13,2)  = 1d0
      trans(13,3)  = 1d0
      trans(13,4)  = -2d0
      trans(13,5)  = 0d0
      trans(13,6)  = 0d0
      trans(13,7)  = -4d0
      trans(13,8)  = 0d0
      trans(13,9)  = 0d0
      trans(13,10) = 0d0
      trans(13,11) = 0d0
      trans(13,12) = 0d0
      trans(13,13) = 6d0
*
*     Rotate PDFs
*
      do a=0,nin(igrid)
         do i=1,13
            pdfout(i-7,a) = 0d0
            do j=1,13
               pdfout(i-7,a) = pdfout(i-7,a) 
     1                       + trans(i,j) * pdfin(j,a) / 12d0
            enddo
         enddo
      enddo
*
      return
      end
