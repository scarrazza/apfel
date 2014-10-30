************************************************************************
*
*     F3charm.f:
*
*     This function returns the value of the charm structure function
*     F3c.
*
************************************************************************
      function F3charm(x)
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
      double precision y
**
*     Output Variables
*
      double precision F3charm
*
*     Rescale the Bjorken's x for the CC structur function in the FFNS
*
      if(ProcessDIS.eq."CC".and.MassScheme(1:4).eq."FFNS")then
         y = x * ( 1d0 + m2th(4) / Q2DIS )
         if(y.gt.1d0) y = 0.9999999999999999d0
      else
         y = x
      endif
*
      if(y.lt.xmin(1).or.y.gt.xmax)then
         write(6,*) "In F3charm.f:"
         write(6,*) "Invalid value of x =",y
         call exit(-10)
      endif
*
*     Select the grid
*
      do igrid=1,ngrid
         if(y.ge.xmin(igrid).and.y.lt.xmin(igrid+1))then
            goto 101
         endif
      enddo
*
*     Interpolation
*
 101  F3charm = 0d0
      n = inter_degree(igrid)
      do alpha=0,nin(igrid)
         F3charm = F3charm + w_int(n,alpha,y) * F3(4,igrid,alpha)
      enddo
      if(dabs(F3charm).le.1d-14) F3charm = 0d0
*
      return
      end
