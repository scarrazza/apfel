************************************************************************
*
*     dxgamma.f:
*
*     This function returns the value of the derivative of the photon 
*     PDF at the for the bjorken variable x using the interpolation.
*
************************************************************************
      function dxgamma(x)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/fph.h"
**
*     Input Variables
*
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
      double precision dxgamma
*
*     Check consistency of the input variable
*
      if(x.lt.xmin(1).or.x.gt.xmax)then
         write(6,*) "In dxgamma.f:"
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
 101  dxgamma = 0d0
      n = inter_degree(igrid)
      do alpha=0,nin(igrid)
         dxgamma = dxgamma + w_int(n,alpha,x) * dfgamma(igrid,alpha)
      enddo
*
      return
      end
