************************************************************************
*
*     wcDYqg.f:
*
*     Set of functions that return the qg part of the NLO correction
*     to the double differential DY distributions.
* 
************************************************************************
*
*     Single Integral in x1 and x2
*
************************************************************************
      function SingleIntegrand_x2_QG(jx,x0)
*
      implicit none
*
      include "../commons/colfact.h"
      include "../commons/ngau.h"
*
      integer jx,igau
      double precision SingleIntegrand_x2_QG,I2
      double precision x2,x10,x20,x0(2)
      double precision E2
      double precision elin
      double precision qgterm
*
      x10 = x0(1)
      x20 = x0(2)
*
      I2 = 0d0
      do igau = 1,ngau
         x2 = Y2(igau)
         E2 = elin(jx,x2)
*
*     The factor 2 is due to the AS convention
*
         QGTERM = 2d0 * TR * E2 * 
     &        ( ( x20**2+(x2-x20)**2)*
     &        log(2d0*(x2-x20)*(1d0-x10)/x10/(x2+x20) )  +
     &        2d0*x20*(x2-x20) )/x2**3d0
         
         I2 = I2 + W2(igau) * QGTERM 
         
      enddo
*
      SingleIntegrand_x2_QG = I2
*     
      return
      end         
*
************************************************************************
*
*     Double integral with subtraction only in x1
*
************************************************************************
      function DoubleIntegrand_QG(ix,jx,x0)
*
      implicit none
*
      include "../commons/colfact.h"
      include "../commons/ngau.h"
*
      integer igau,jgau,ix,jx
*
      double precision doubleintegrand_QG,I3
      double precision x1,x2,x10,x20,x0(2)
      double precision E1,E10,E2
      double precision elin
      double precision QGterm,QGTERMPLUS
      double precision tmp
      double precision HC_Y_BIS,GC_Y_BIS
*
      x10 = x0(1)
      x20 = x0(2)
      E10 = elin(ix,x10) ! vanishes if Ix.ne..ixP(1),ixP(1)+1
*
      I3 = 0d0
      do igau = 1,ngau
         tmp = 0d0
         x1 = Y1(igau)
         E1 = elin(ix,x1)
         
         do jgau = 1,ngau
            x2 = Y2(jgau)
            E2 = elin(jx,x2)
            
            QGTERM = 2d0 * TR * E1 * E2 *  hc_Y_BIS(x2,x1,x20,x10)
            QGTERMPLUS = 2d0 * TR * E2
     &           * (  GC_Y_BIS(x2,x1,x20,x10) * E1-
     &           GC_Y_BIS(x2,x10,x20,x10)* E10  )
     $           /(x1-x10)
            
            tmp = tmp + w1(jgau)* (QGTERM + QGTERMPLUS)
         enddo
         I3 = I3 + w2(igau)*tmp
      enddo
*
      DoubleIntegrand_QG = I3
*
      return
      end 
*
************************************************************************
*
*     Subtraction term
*
************************************************************************
      function DoubleIntegrand_sub_QG(jx,x0)
*
      implicit none
*
      include "../commons/colfact.h"
      include "../commons/ngaus.h"
*
      integer igau,jgau,jx
      double precision doubleintegrand_sub_QG,I3
      double precision x1,x2,x10,x20,x0(2)
      double precision E2
      double precision elin
      double precision QGTERMPLUS
      double precision tmp
      double precision GC_Y_BIS
*
      x10 = x0(1)
      x20 = x0(2)
*
      I3 = 0d0
      tmp  = 0d0
      do igau = 1,ngaus 
         tmp = 0d0
         x1 = Y1S(igau)
         do jgau = 1,ngaus
            x2 = Y2S(jgau)
            E2 = elin(jx,x2)
            QGTERMPLUS = 2d0 * TR * E2
     &           * GC_Y_BIS(x2,x10,x20,x10)/(x1-x10)
            
            tmp = tmp + w1S(jgau)* QGTERMPLUS
         enddo
         I3 = I3 + w2S(igau)*tmp
      enddo
*
      DoubleIntegrand_sub_QG = I3
*
      return
      end 
*
************************************************************************
*
*     gq part
*
************************************************************************
*
*     Single integral bewteen x01 and 1
*
************************************************************************
      function SingleIntegrand_x1_GQ(ix,x0)
*
      implicit none
*
      include "../commons/colfact.h"
      include "../commons/ngau.h"
*
      integer igau,ix
      double precision SingleIntegrand_x1_GQ,I1
      double precision x1,x10,x20,x0(2)
      double precision E1,elin
      double precision gqterm
*
      x10 = x0(1)
      x20 = x0(2)
*    
      I1 = 0d0
      do igau = 1,ngau
         x1 = Y1(igau)
         E1 = elin(ix,x1)
*
*     The factor 2 is due to the AS convention
*
         GQTERM = 2d0 * TR * E1 * 
     &        ( ( x10**2+(x1-x10)**2)*
     &        log(2d0*(x1-x10)*(1d0-x20)/x20/(x1+x10) ) +
     &        2d0*x10*(x1-x10)
     &        )/x1**3
         
         I1 = I1 + W1(igau) * GQTERM 
      enddo
*     
      SingleIntegrand_x1_GQ = I1
*     
      return
      end
*
************************************************************************
*
*     Double integral with subtraction only in x
*
************************************************************************
      function DoubleIntegrand_GQ(ix,jx,x0)
*
      implicit none
*
      include "../commons/colfact.h"
      include "../commons/ngau.h"
*
      integer igau,jgau,ix,jx

      double precision doubleintegrand_GQ,I3
      double precision x1,x2,x10,x20,x0(2)
      double precision E1,E2,E20
      double precision elin
      double precision GQterm,GQTERMPLUS
      double precision tmp
      double precision HC_Y_BIS,GC_Y_BIS
*
      x10 = x0(1)
      x20 = x0(2)
      E20 = elin(jx,x20) ! vanishes if jx.ne.ixP(2),ixP(2)+1
*
      I3 = 0d0
      do igau = 1,ngau
         tmp = 0d0
         x1 = Y1(igau)
         E1 = elin(ix,x1)         
         do jgau = 1,ngau
            x2 = Y2(jgau)
            E2 = elin(jx,x2)
            GQTERM = 2d0 * TR * E1 * E2 *  HC_Y_BIS(x1,x2,x10,x20)
            GQTERMPLUS = 2d0 * TR * E1 
     &           * (  GC_Y_BIS(x1,x2,x10,x20) * E2-
     &           GC_Y_BIS(x1,x20,x10,x20)* E20  )
     $           /(x2-x20)               
            tmp = tmp + w1(jgau)* (GQTERM + GQTERMPLUS)
         enddo         
         I3 = I3 + w2(igau)*tmp
      enddo
*
      DoubleIntegrand_GQ = I3
*
      return
      end
*
************************************************************************
*
*     Subtraction term
*
************************************************************************      
      function DoubleIntegrand_sub_GQ(ix,x0)
*
      implicit none
*
      include "../commons/colfact.h"
      include "../commons/ngaus.h"
*
      integer igau,jgau,ix
      double precision doubleintegrand_sub_GQ,I3
      double precision x1,x2,x10,x20,x0(2)
      double precision E1
      double precision elin
      double precision GQTERMPLUS
      double precision tmp
      double precision GC_Y_BIS
*
      x10 = x0(1)
      x20 = x0(2)
*
      I3 = 0d0
*
      do igau = 1,ngaus
         tmp = 0d0
         x1 = Y1S(igau)
         E1 = elin(ix,x1)     
         do jgau = 1,ngaus
            x2 = Y2S(jgau)                           
            GQTERMPLUS = 2d0 * TR * E1 
     &           * GC_Y_BIS(x1,x20,x10,x20)/(x2-x20)
            tmp = tmp + w1S(jgau)* GQTERMPLUS
         enddo
         I3 = I3 + w2S(igau)*tmp
      enddo
*
      DoubleIntegrand_sub_GQ = I3
*
      return
      end
*
************************************************************************
      double precision function hc_Y_BIS(x1,x2,x10,x20)
*
      implicit none
*
      double precision x1,x2,x10,x20
*
      hc_Y_BIS = 2d0*x10*x20*(x10*x20+x1*x2)*
     &     (x1*x10*x2**2+x10*x20**2*x1+2d0*x10**2*x2*x20)/
     &     x1**2/x2**2/(x1*x20+x2*x10)**3
*
      end function
*
************************************************************************
      double precision function GC_Y_BIS(x1,x2,x10,x20)
*
      implicit none
*
      double precision x1,x2,x10,x20
*
      GC_Y_BIS=2d0*x20*(x10**2*x20**2+(x10*x20-x1*x2)**2)
     &     *(x10*x20+x1*x2)
     &     /x1**3/x2**2/(x1*x20+x2*x10)/(x2+x20)
*
      end function
