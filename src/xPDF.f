************************************************************************
*
*     xPDF.f:
*
*     This function returns the value of the i-th PDF in the physical
*     basis at the final scale and for the bjorken variable x using 
*     the interpolation.
*
************************************************************************
      function xPDF(i,x)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/fph.h"
**
*     Input Variables
*
      integer i
      double precision x
**
*     Input Variables
*
      integer n
      integer alpha
      double precision w_int
**
*     Output Variables
*
      double precision xPDF
*
*     Check consistency of the input variables
*
      if(i.lt.-6.or.i.gt.6)then
         write(6,*) "In xPDF.f:"
         write(6,*) "Invalid PDF index, i =",i
         call exit(-10)
      endif
*
      if(x.lt.xmin(1).or.x.gt.xmax)then
         write(6,*) "In xPDF.f:"
         write(6,*) "Invalid value of x =",x
         call exit(-10)
      endif
*
*     Select the grid
*
      do igrid=1,ngrid
         if(x.ge.xmin(igrid).and.x.lt.xmin(igrid+1))then
            goto 101
         endif
      enddo
*
*     Interpolation
*
 101  xPDF = 0d0
      n = inter_degree(igrid)
      do alpha=0,nin(igrid)
         xPDF = xPDF + w_int(n,alpha,x) * fph(igrid,i,alpha)
      enddo
      if(dabs(xPDF).le.1d-14) xPDF = 0d0
*
      return
      end
