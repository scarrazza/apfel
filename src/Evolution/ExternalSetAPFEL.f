************************************************************************
*
*     Dummy subroutines for the external set to ensure the compilation
*
************************************************************************
      subroutine ExternalSetAPFEL(x,Q,xf)
*
      implicit none
**
*     Input Variables
*
      double precision x,Q
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
*
************************************************************************
      subroutine ExternalSetAPFEL1(x,Q,xf)
*
      implicit none
**
*     Input Variables
*
      double precision x,Q
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
*
************************************************************************
      subroutine ExternalSetAPFELLept(x,Q,xl,xf)
*
      implicit none
**
*     Input Variables
*
      double precision x,Q
**
*     Internal Variables
*
      integer ipdf,ilep
**
*     Output Variables
*
      double precision xl(-3:3),xf(-6:7)
*
*     Set to zero
*
      do ipdf=-6,7
         xf(ipdf) = 0d0
      enddo
      do ilep=-3,3
         xl(ilep) = 0d0
      enddo
*
      return
      end
