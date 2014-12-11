************************************************************************
*
*     Dummy subroutine for the external set to ensure the compilation
*
************************************************************************
      subroutine ExternalSetAPFEL(x,xf)
*
      implicit none
**
*     Input Variables
*
      double precision x
**
*     Internal Variables
*
      integer ipdf
**
*     Output Variables
*
      double precision xf(-6:7)
*
*     Set to zero
*
      do ipdf=-6,7
         xf(ipdf) = 0d0
      enddo
*
      return
      end
