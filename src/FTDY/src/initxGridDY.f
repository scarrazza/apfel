************************************************************************
*
*     initxGridDY.f:
*
*     Initialize th x-space grid for Drell-Yan.
*
************************************************************************
      subroutine initxGridDY
*
      implicit none
*
      include "../commons/mxgridsizeDY.h"
      include "../commons/xgridDY.h"
      include "../commons/xxDY.h"
*     
      integer ix
      double precision hx,hxlim
*
      if(nxDY.gt.mxgridsizeDY)then
         write(6,*) "ERROR: in initxGrid.f:"
         write(6,*) "Number of points exceeds maximum allowed"
         write(6,*) "nxDY =",nxDY,", mxgridsizeDY =",mxgridsizeDY
         call exit(-10)
      endif
*
*     Grid uniformly spaced in sqrt( log_10 (1/x) ).
*
      hxlim = - dsqrt( dlog10( xmaxDY / xminDY ) )
*
      do ix=1,nxDY
         hx = hxlim * ( 1d0 - ( dble( ix - 1d0 ) ) / ( nxDY - 1d0 ) )
         xxDY(ix) = 10d0**( - hx * hx )
      enddo
*
      if(xxDY(nxDY).ne.xmaxDY)then
         write(6,*) "ERROR: in initxGrid.f:"
         write(6,*) "Wrong value of xmaxDY:",xxDY(nxDY),xmaxDY
         call exit(-10)
      endif
*     
      return
      end
