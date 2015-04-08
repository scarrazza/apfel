c$$$************************************************************************
c$$$*
c$$$*     Dummy subroutines for the external set to ensure the compilation
c$$$*
c$$$************************************************************************
c$$$      subroutine ExternalSetAPFEL(x,Q,xf)
c$$$*
c$$$      implicit none
c$$$**
c$$$*     Input Variables
c$$$*
c$$$      double precision x,Q
c$$$**
c$$$*     Internal Variables
c$$$*
c$$$      integer ipdf
c$$$**
c$$$*     Output Variables
c$$$*
c$$$      double precision xf(-6:7)
c$$$*
c$$$*     Set to zero
c$$$*
c$$$      do ipdf=-6,7
c$$$         xf(ipdf) = 0d0
c$$$      enddo
c$$$*
c$$$      return
c$$$      end
c$$$*
c$$$************************************************************************
c$$$      subroutine ExternalSetAPFEL1(x,Q,xf)
c$$$*
c$$$      implicit none
c$$$**
c$$$*     Input Variables
c$$$*
c$$$      double precision x,Q
c$$$**
c$$$*     Internal Variables
c$$$*
c$$$      integer ipdf
c$$$**
c$$$*     Output Variables
c$$$*
c$$$      double precision xf(-6:7)
c$$$*
c$$$*     Set to zero
c$$$*
c$$$      do ipdf=-6,7
c$$$         xf(ipdf) = 0d0
c$$$      enddo
c$$$*
c$$$      return
c$$$      end
c$$$*
c$$$************************************************************************
c$$$      subroutine ExternalSetAPFELLept(x,Q,irep,xl,xf)
c$$$*
c$$$      implicit none
c$$$**
c$$$*     Input Variables
c$$$*
c$$$      integer irep
c$$$      double precision x,Q
c$$$**
c$$$*     Internal Variables
c$$$*
c$$$      integer ipdf,ilep
c$$$**
c$$$*     Output Variables
c$$$*
c$$$      double precision xl(-3:3),xf(-6:7)
c$$$*
c$$$*     Set to zero
c$$$*
c$$$      do ipdf=-6,7
c$$$         xf(ipdf) = 0d0
c$$$      enddo
c$$$      do ilep=-3,3
c$$$         xl(ilep) = 0d0
c$$$      enddo
c$$$*
c$$$      return
c$$$      end
