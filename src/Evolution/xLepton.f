************************************************************************
*
*     xLepton.f:
*
*     This function returns the value of the i-th lepton in the physical
*     basis at the final scale and for the bjorken variable x using 
*     the interpolation.
*
************************************************************************
      function xLepton(i,x)
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
      double precision xLepton
*
*     Check consistency of the input variables
*
      if(i.lt.-3.or.i.gt.3)then
         write(6,*) "In xLepton.f:"
         write(6,*) "Invalid Lepton index, i =",i
         call exit(-10)
      endif
*
      if(x.lt.xmin(1).or.x.gt.xmax)then
         write(6,*) "In xLepton.f:"
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
 101  xLepton = 0d0
      n = inter_degree(igrid)
      do alpha=0,nin(igrid)
         xLepton = xLepton + w_int(n,alpha,x) * flepton(igrid,i,alpha)
      enddo
      if(dabs(xLepton).le.1d-12) xLepton = 0d0
*
      return
      end
*
************************************************************************
*
*     Interpolation on the joint x-space grid.
*
************************************************************************
      function xLeptonj(i,x)
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
      integer alpha!,jgrid
      double precision w_int_gen
**
*     Output Variables
*
      double precision xLeptonj
*
*     Check consistency of the input variables
*
      if(i.lt.-3.or.i.gt.3)then
         write(6,*) "In xLepton.f:"
         write(6,*) "Invalid Lepton index, i =",i
         call exit(-10)
      endif
*
      if(x.lt.xmin(1).or.x.gt.xmax)then
         write(6,*) "In xLepton.f:"
         write(6,*) "Invalid value of x =",x
         call exit(-10)
      endif
*
*     Interpolation
*
      xLeptonj = 0d0
      n = inter_degree(0)
      do alpha=0,nin(0)
         xLeptonj = xLeptonj + w_int_gen(n,alpha,x) * flepton(0,i,alpha)
      enddo
      if(dabs(xLeptonj).le.1d-12) xLeptonj = 0d0
*
      return
      end
