************************************************************************
*
*     PDFevQCD2phys.f:
*
*     This routine converts PDFs from the QDC evolution basis to the 
*     physical basis:
*
*     Evolution basis:
*       1   2   3   4   5   6   7   8   9  10  11  12  13
*      Sg   g   V  V3  V8 V15 V24 V35  T3  T8 T15 T24 T35
*
*     Physical basis:
*      -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6 
*      tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
*
************************************************************************
      subroutine PDFevQCD2phys(pdfin,pdfout)
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
      double precision evln2lha(13,13)
**
*     Output Variables
*
      double precision pdfout(-6:6,0:nint_max)
*
*     Define rotation matrix
*
      do i=1,13
         do j=1,13
            evln2lha(i,j) = 0d0
         enddo
      enddo
      do i=1,13
         evln2lha(i,1) = 10.d0
      enddo
      evln2lha(7,1)  = 0.d0
      evln2lha(7,2)  = 120.d0
      do i=1,6 
         evln2lha(i,3) = -10.d0
      enddo
      do i=8,13 
         evln2lha(i,3) = 10.d0
      enddo
      evln2lha(5,4)  = -30.d0
      evln2lha(6,4)  =  30.d0
      evln2lha(8,4)  = -30.d0
      evln2lha(9,4)  =  30.d0
      evln2lha(4,5)  =  20.d0
      evln2lha(5,5)  = -10.d0
      evln2lha(6,5)  = -10.d0
      evln2lha(8,5)  =  10.d0
      evln2lha(9,5)  =  10.d0
      evln2lha(10,5) = -20.d0
      evln2lha(3,6)  =  15.d0
      evln2lha(4,6)  = - 5.d0
      evln2lha(5,6)  = - 5.d0
      evln2lha(6,6)  = - 5.d0
      evln2lha(8,6)  =   5.d0
      evln2lha(9,6)  =   5.d0
      evln2lha(10,6) =   5.d0
      evln2lha(11,6) = -15.d0
      evln2lha(2,7)  =  12.d0
      evln2lha(3,7)  = - 3.d0
      evln2lha(4,7)  = - 3.d0
      evln2lha(5,7)  = - 3.d0
      evln2lha(6,7)  = - 3.d0
      evln2lha(8,7)  =   3.d0
      evln2lha(9,7)  =   3.d0
      evln2lha(10,7) =   3.d0
      evln2lha(11,7) =   3.d0
      evln2lha(12,7) = -12.d0
      evln2lha(1,8)  =  10.d0
      evln2lha(2,8)  = - 2.d0
      evln2lha(3,8)  = - 2.d0
      evln2lha(4,8)  = - 2.d0
      evln2lha(5,8)  = - 2.d0
      evln2lha(6,8)  = - 2.d0
      evln2lha(8,8)  =   2.d0
      evln2lha(9,8)  =   2.d0
      evln2lha(10,8) =   2.d0
      evln2lha(11,8) =   2.d0
      evln2lha(12,8) =   2.d0
      evln2lha(13,8) = -10.d0
      evln2lha(5,9)  =  30.d0
      evln2lha(6,9)  = -30.d0
      evln2lha(8,9)  = -30.d0
      evln2lha(9,9)  =  30.d0
      evln2lha(4,10) = -20.d0
      evln2lha(5,10) =  10.d0
      evln2lha(6,10) =  10.d0
      evln2lha(8,10) =  10.d0
      evln2lha(9,10) =  10.d0
      evln2lha(10,10)= -20.d0
      evln2lha(3,11) = -15.d0
      evln2lha(4,11) =   5.d0
      evln2lha(5,11) =   5.d0
      evln2lha(6,11) =   5.d0
      evln2lha(8,11) =   5.d0
      evln2lha(9,11) =   5.d0
      evln2lha(10,11)=   5.d0
      evln2lha(11,11)= -15.d0
      evln2lha(2,12) = -12.d0
      evln2lha(3,12) =   3.d0
      evln2lha(4,12) =   3.d0
      evln2lha(5,12) =   3.d0
      evln2lha(6,12) =   3.d0
      evln2lha(8,12) =   3.d0
      evln2lha(9,12) =   3.d0
      evln2lha(10,12)=   3.d0
      evln2lha(11,12)=   3.d0
      evln2lha(12,12)= -12.d0
      evln2lha(1,13) = -10.d0
      evln2lha(2,13) =   2.d0
      evln2lha(3,13) =   2.d0
      evln2lha(4,13) =   2.d0
      evln2lha(5,13) =   2.d0
      evln2lha(6,13) =   2.d0
      evln2lha(8,13) =   2.d0
      evln2lha(9,13) =   2.d0
      evln2lha(10,13)=   2.d0
      evln2lha(11,13)=   2.d0
      evln2lha(12,13)=   2.d0
      evln2lha(13,13)= -10.d0
      do i=1,13
         do j=1,13
            evln2lha(i,j) = evln2lha(i,j) / 120.d0
         enddo
      enddo
*
*     Rotate PDFs
*
      do a=0,nin(igrid)
         do i=1,13
            pdfout(i-7,a) = 0d0
            do j=1,13
               pdfout(i-7,a) = pdfout(i-7,a) + evln2lha(i,j)*pdfin(j,a)
            enddo
         enddo
      enddo
*
      return
      end
