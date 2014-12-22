************************************************************************
*
*     dxPDF.f:
*
*     This function returns the value of the derivative of the i-th PDF
*     in the physical basis at for the bjorken variable x using 
*     the interpolation.
*
************************************************************************
      function dxPDF(i,x)
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
*     Internal Variables
*
      integer n
      integer alpha
      double precision w_int
**
*     Output Variables
*
      double precision dxPDF
*
*     Check consistency of the input variables
*
      if(i.lt.-6.or.i.gt.6)then
         write(6,*) "In dxPDF.f:"
         write(6,*) "Invalid PDF index, i =",i
         call exit(-10)
      endif
*
      if(x.lt.xmin(1).or.x.gt.xmax)then
         write(6,*) "In dxPDF.f:"
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
 101  dxPDF = 0d0
      n = inter_degree(igrid)
      do alpha=0,nin(igrid)
         dxPDF = dxPDF + w_int(n,alpha,x) * dfph(igrid,i,alpha)
      enddo
*
      return
      end
