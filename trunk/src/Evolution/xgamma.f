************************************************************************
*
*     xgamma.f:
*
*     This function returns the value of the photon PDF at the final 
*     scale and for the bjorken variable x using the interpolation.
*
************************************************************************
      function xgamma(x)
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
      double precision xgamma
*
*     Check consistency of the input variable
*
      if(x.lt.xmin(1).or.x.gt.xmax)then
         write(6,*) "In xgamma.f:"
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
 101  xgamma = 0d0
      n = inter_degree(igrid)
      do alpha=0,nin(igrid)
         xgamma = xgamma + w_int(n,alpha,x) * fgamma(igrid,alpha)
      enddo
      if(dabs(xgamma).le.1d-14) xgamma = 0d0
*
      return
      end
*
************************************************************************
*
*     Interpolation on the joint x-space grid.
*
************************************************************************
      function xgammaj(x)
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
      integer alpha!,jgrid
      double precision w_int_gen
**
*     Output Variables
*
      double precision xgammaj
*
*     Check consistency of the input variable
*
      if(x.lt.xmin(1).or.x.gt.xmax)then
         write(6,*) "In xgamma.f:"
         write(6,*) "Invalid value of x =",x
         call exit(-10)
      endif
*
*     Select interpolation degree
*
c      do igrid=1,ngrid
c         if(x.ge.xmin(jgrid).and.x.lt.xmin(jgrid+1))then
c            goto 101
c         endif
c      enddo
c 101  n = inter_degree(jgrid)
      n = inter_degree(0)
*
*     Interpolation
*
      xgammaj = 0d0
      do alpha=0,nin(0)
         xgammaj = xgammaj + w_int_gen(n,alpha,x) * fgamma(0,alpha)
      enddo
      if(dabs(xgammaj).le.1d-14) xgammaj = 0d0
*
      return
      end
