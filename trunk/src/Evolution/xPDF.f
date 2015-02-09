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
*     Internal Variables
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
*
************************************************************************
*
*     Interpolation on the joint x-space grid.
*
************************************************************************
      function xPDFj(i,x)
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
      double precision xPDFj
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
*     Select the interpolation degree
*
c      do jgrid=1,ngrid
c         if(x.ge.xmin(jgrid).and.x.lt.xmin(jgrid+1))then
c            goto 101
c         endif
c      enddo
c 101  n = inter_degree(jgrid)
      n = inter_degree(0)
*
*     Interpolation
*
      xPDFj = 0d0
      do alpha=0,nin(0)
         xPDFj = xPDFj + w_int_gen(n,alpha,x) * fph(0,i,alpha)
      enddo
      if(dabs(xPDFj).le.1d-14) xPDFj = 0d0
*
      return
      end
*
************************************************************************
*
*     The following routine returns all PDFs in the physical basis at 
*     the final scale and for the bjorken variable x using the
*     interpolation on the joint grid. No photon included.
*
************************************************************************
      subroutine xPDFall(x,xf)
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
      integer ipdf
      double precision w_int_gen,wgt(0:nint_max)
**
*     Output Variables
*
      double precision xf(-6:6)
*
*     Check consistency of the input variables
*
      if(x.lt.xmin(1).or.x.gt.xmax)then
         write(6,*) "In xPDF.f:"
         write(6,*) "Invalid value of x =",x
         call exit(-10)
      endif
*
*     Interpolation
*
      n = inter_degree(0)
      do alpha=0,nin(0)
         wgt(alpha) = w_int_gen(n,alpha,x)
      enddo
*
      do ipdf=-6,6
         xf(ipdf) = 0d0
         do alpha=0,nin(0)
            xf(ipdf) = xf(ipdf) + wgt(alpha) * fph(0,ipdf,alpha)
         enddo
         if(dabs(xf(ipdf)).le.1d-14) xf(ipdf) = 0d0
      enddo
*
      return
      end
