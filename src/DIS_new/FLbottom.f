************************************************************************
*
*     FLbottom.f:
*
*     This function returns the value of the bottom structure function
*     FLb.
*
************************************************************************
      function FLbottom(x)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/StructureFunctions.h"
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
      double precision FLbottom
*
      if(x.lt.xmin(1).or.x.gt.xmax)then
         write(6,*) "In FLbottom.f:"
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
 101  FLbottom = 0d0
      n = inter_degree(igrid)
      do alpha=0,nin(igrid)
         FLbottom = FLbottom + w_int(n,alpha,x) * FL(5,igrid,alpha)
      enddo
      if(dabs(FLbottom).le.1d-14) FLbottom = 0d0
*
      return
      end
