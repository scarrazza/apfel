************************************************************************
*
*     FLtop.f:
*
*     This function returns the value of the inclusive structure function
*     FL.
*
************************************************************************
      function FLtop(x)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/StructureFunctions.h"
      include "../commons/ProcessDIS.h"
      include "../commons/MassScheme.h"
      include "../commons/m2th.h"
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
      double precision FLtop
*
      if(x.lt.xmin(1).or.x.gt.xmax)then
         write(6,*) "In FLtop.f:"
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
 101  FLtop = 0d0
      n = inter_degree(igrid)
      do alpha=0,nin(igrid)
         FLtop = FLtop + w_int(n,alpha,x) * FL(6,igrid,alpha)
      enddo
      if(dabs(FLtop).le.1d-14) FLtop = 0d0
*
      return
      end
