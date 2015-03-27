************************************************************************
*
*     wcDYqq.f:
*
*     Set of functions that return the qq part of the NLO correction
*     to the double differential DY distributions.
* 
************************************************************************
      function ZeroIntegrand_QQ(x0)
*
      implicit none
*
      include "../commons/colfact.h"
      include "../commons/consts.h"
*
      double precision zerointegrand_qq
      double precision qqbterm,x0(2),x1,x2
*
      real*8 ddilog
      external ddilog
*
      x1 = x0(1)
      x2 = x0(2)
*
*     The factor 2 is due to the AS convention
*
      qqbterm = 2d0 * cf*
     &     (pi**2/3d0-8d0+2d0*(ddilog(x1)+ddilog(x2))+
     &     log(1-x1)**2+log(1-x2)**2+
     &     2d0*log(x1/(1d0-x1))*log(x2/(1d0-x2))) 
*
      zerointegrand_qq = qqbterm
*
      return
      end 
*
************************************************************************
*
*     Single Integral in x2
*
************************************************************************
      function SingleIntegrand_x2_QQ(jx,x0)
*
      implicit none
*
      include "../commons/colfact.h"
      include "../commons/ngau.h"
*      
      integer igau,jx
      double precision SingleIntegrand_x2_QQ,I2
      double precision x2,x10,x20,x0(2)
      double precision E2,E20
      double precision elin
      double precision qqbterm,qqbtermplus
*
      x10 = x0(1)
      x20 = x0(2)
*
*     notice that e20 is zero if jx.ne.ixp(2) or .ixp(2)+1
*     in these cases the subtraction disappears.
*
      E20 = elin(jx,x20)
      I2  = 0d0
      do igau = 1,ngau
         x2 = Y2(igau)
         E2 = elin(jx,x2)
*
         QQBTERM = E2 * (1d0/x2 - x20/x2**2d0
     &        - (x20**2d0+x2**2d0)/(x2**2d0*(x2-x20))
     &        * dlog(x20/x2))
*
         QQBTERMPLUS = 
     &        log(1d0-x20/x2)/(x2-x20) * 
     &        ( (1d0+x20**2/x2**2)* E2 - 2d0*E20 ) +
     &        1d0/(x2-x20)* 
     &        ((1d0+x20**2/x2**2)*
     &        log((2d0*x20*(1d0-x10))/(x10*(x2+x20)))* E2
     &        -    2d0*log((1d0-x10)/x10 )* E20 )
*
         I2 = I2 + W2(igau) * (QQBTERM + QQBTERMPLUS)
      enddo      
*
*     The factor 2 is due to the AS convention
*
      SingleIntegrand_x2_QQ = 2d0*CF * I2
*     
      return
      end
*
************************************************************************
*
*     Single Integral in x1
*
************************************************************************
      function SingleIntegrand_x1_QQ(ix,x0)
*
      implicit none
*
      include "../commons/colfact.h"
      include "../commons/ngau.h"
*
      integer igau,ix
      double precision SingleIntegrand_x1_QQ,I1
      double precision x1,x10,x20,x0(2)
      double precision E1,E10
      double precision elin
      double precision qqbterm,qqbtermplus
*
      x10 = x0(1)
      x20 = x0(2)
*
*     notice that e10 is zero if jx.ne.ixp(1) or .ixp(1)+1
*     in these cases the subtraction disappears.
*
      E10 = elin(ix,x10)
      I1 = 0d0
      do igau = 1,ngau
         x1 = Y1(igau)
         E1 = elin(ix,x1)
*
         QQBTERM = E1 * (1d0/x1 - x10/x1**2d0
     &        - (x10**2+x1**2)/(x1**2*(x1-x10))
     &        * log(x10/x1) )
*
         QQBTERMPLUS =  
     &        log(1d0-x10/x1)/(x1-x10) * 
     &        ( (1d0+x10**2/x1**2)*E1 - 2d0 * E10)
     &        + 1d0/(x1-x10)* 
     &        ((1d0+x10**2/x1**2)
     &        *log((2d0*x10*(1d0-x20))/(x20*(x1+x10)))* E1
     &        - 2d0*log((1d0-x20)/x20 )* E10 ) 
*
         I1 = I1 + W1(igau) * (QQBTERM + QQBTERMPLUS)
      enddo      
*
*     The factor 2 is due to the AS convention
*
      SingleIntegrand_x1_QQ = 2d0 * CF * I1
*
      return
      end
*
************************************************************************
*
*     Single integral subtraction
*
************************************************************************
      function SingleIntegrand_sub_x1(x0)
*
      implicit none
*
      include "../commons/colfact.h"
      include "../commons/ngaus.h"
*
      integer igau
      double precision SingleIntegrand_sub_x1,I1
      double precision x1,x10,x20,x0(2)
      double precision qqbtermplus
*
      x10 = x0(1)
      x20 = x0(2)
*
*     this is the same for both rapidity and xF distributions
*
      I1 = 0d0
      do igau = 1,ngaus
         x1 = Y1S(igau)
*
         QQBTERMPLUS = 
     &        2d0 * log(1d0-x10/x1)/(x1-x10)  
     &        + 2d0/(x1-x10)*  log((1d0-x20)/x20)
*
         I1 = I1 + W1S(igau) * QQBTERMPLUS
      enddo
*
*     The factor 2 is due to the AS convention
*
      SingleIntegrand_sub_x1 = 2d0 * CF * I1
*
      return
      end
*
************************************************************************
      function SingleIntegrand_sub_x2(x0)
*
      implicit none
*
      include "../commons/colfact.h"
      include "../commons/ngaus.h"
*      
      integer igau
      double precision SingleIntegrand_sub_x2,I2
      double precision x2,x10,x20,x0(2)
      double precision qqbtermplus
*
      x10 = x0(1)
      x20 = x0(2)
*
*     this is the same for both rapidity and xF distributions
*
      I2 = 0d0
      do igau = 1,ngaus
         x2 = Y2S(igau)
*
         QQBTERMPLUS = 
     &        2d0 * log(1d0-x20/x2)/(x2-x20) 
     &        +  2d0 * log((1d0-x10)/x10)/(x2-x20) 
*
         I2 = I2 + W2S(igau) * QQBTERMPLUS
      enddo
*
*     The factor 2 is due to the AS convention
*
      SingleIntegrand_sub_x2 = 2d0*CF * I2
*
      return
      end
*
************************************************************************
*
*     Double integral in x1 and x2
*
************************************************************************
      function DoubleIntegrand_QQ(ix,jx,x0)
*
      implicit none
*
      include "../commons/colfact.h"
      include "../commons/ngau.h"
*      
      integer igau,jgau,ix,jx
      double precision doubleintegrand_QQ,I3
      double precision x1,x2,x10,x20,x0(2)
      double precision E1,E10,E2,E20
      double precision elin
      double precision qqbterm,qqbtermplus
      double precision tmp
      double precision ga_bis,ha_bis
*
      x10 = x0(1)
      x20 = x0(2)
      E10 = elin(ix,x10)
      E20 = elin(jx,x20)
*
      I3 = 0d0
      tmp  = 0d0
      do igau = 1,ngau
         x1 = Y1(igau)
         E1 = elin(ix,x1)
         tmp = 0d0
         do jgau = 1,ngau        
            x2 = Y2(jgau)
            E2 = elin(jx,x2)
            QQBTERM = E1 * E2 * ha_bis(x1,x2,x10,x20) ! HA
            QQBTERMPLUS =       !GA
     &           1d0/(x1-x10)/(x2-x20)*
*     NO SUB
     &           ( ga_bis(x1,x2,x10,x20) * E1 * E2 
*     SUB x1: vanishes if ix.ne.ixP(1) or ixP(1)+1
     &           - ga_bis(x10,x2,x10,x20)  * E10 * E2 
*     SUB x2: vanishes if jx.ne.ixP(2) or ixP(2)+1
     &           - ga_bis(x1,x20,x10,x20)  * E1  * E20 
*     DOUBLE SUB: vanishes in the previous 2 cases
     &           + ga_bis(x10,x20,x10,x20) * E10 * E20 ) 
            tmp = tmp + w2(jgau)* (QQBTERM + QQBTERMPLUS)
         enddo
         I3 = I3 + w1(igau)*tmp
      enddo
*
*     The factor 2 is due to the AS convention
* 
      DoubleIntegrand_QQ = 2d0 * CF * I3
*
      return
      end
*
************************************************************************
*
*     Double integral in x1 and x2 - Subtraction terms 
*     (dont depend on y or xf distributions)
*
************************************************************************
      function DoubleIntegrand_sub1_QQ(x0)
*
      implicit none
*
      include "../commons/colfact.h"
      include "../commons/ngaus.h"
*      
      integer igau,jgau
      double precision doubleintegrand_sub1_QQ,I3
      double precision x1,x2,x10,x20,x0(2)
      double precision qqbtermplus
      double precision tmp
*
      x10 = x0(1)
      x20 = x0(2)
*
      I3 = 0d0
      do igau = 1,ngaus
         x1 = Y1S(igau)
         tmp = 0d0
         do jgau = 1,ngaus        
            x2 = Y2S(jgau)
            QQBTERMPLUS = 2d0/(x1-x10)/(x2-x20) 
            tmp = tmp + w2s(jgau)* QQBTERMPLUS
         enddo
         I3 = I3 + w1s(igau)*tmp
      enddo
*
*     The factor 2 is due to the AS convention
* 
      DoubleIntegrand_sub1_QQ = 2d0 * CF * I3
*
      return
      end
*
************************************************************************
      function DoubleIntegrand_sub2_QQ(ix,x0)
*
      implicit none
*
      include "../commons/colfact.h"
      include "../commons/ngaus.h"
*
      integer igau,jgau,ix
      double precision doubleintegrand_sub2_QQ,I3
      double precision x1,x2,x10,x20,x0(2)
      double precision E1,E10
      double precision elin
      double precision qqbtermplus
      double precision tmp
*
      x10 = x0(1)
      x20 = x0(2)
      E10 = elin(ix,x10)
*
      I3 = 0d0
      do igau = 1,ngaus
         x1 = Y1S(igau)
         E1 = elin(ix,x1)
         tmp = 0d0
         do jgau = 1,ngaus        
            x2 = Y2S(jgau)
            QQBTERMPLUS =       !GA
     &           1d0/(x1-x10)/(x2-x20)*
     &           ( (1d0+x10**2/x1**2)* E1 
     &           - 2d0* E10 ) 
            tmp = tmp + w2s(jgau)* QQBTERMPLUS
         enddo
         I3 = I3 + w1s(igau)*tmp
      enddo      
*
*     The factor 2 is due to the AS convention
* 
      DoubleIntegrand_sub2_QQ = 2d0 * CF * I3
*
      return
      end
*
************************************************************************
      function DoubleIntegrand_sub3_QQ(jx,x0)
*
      implicit none
*
      include "../commons/colfact.h"
      include "../commons/ngaus.h"
*
      integer igau,jgau,jx
      double precision doubleintegrand_sub3_QQ,I3
      double precision x1,x2,x10,x20,x0(2)
      double precision E2,E20
      double precision elin
      double precision qqbtermplus
      double precision tmp
*
      x10 = x0(1)
      x20 = x0(2)
      E20 = elin(jx,x20)
*
      I3 = 0d0
      do igau = 1,ngaus
         x1 = Y1S(igau)
         tmp = 0d0
         do jgau = 1,ngaus        
            x2 = Y2S(jgau)
            E2 = elin(jx,x2)
            QQBTERMPLUS =                 
     &           1d0/(x1-x10)/(x2-x20)*
     &           ( (1d0+x20**2/x2**2)* E2 
     &           - 2d0 * E20 )             
            tmp = tmp + w2s(jgau)* QQBTERMPLUS
         enddo
         I3 = I3 + w1s(igau)*tmp
      enddo
*
*     The factor 2 is due to the AS convention
* 
      DoubleIntegrand_sub3_QQ = 2d0 * CF * I3
*
      return
      end
*
************************************************************************
*
*     HA
*
************************************************************************
      double precision function ha_bis(x1,x2,x10,x20)
*
      implicit none
*
      double precision x1,x2,x10,x20
*
      ha_bis = -4d0 * ( x10*x20*(x10*x20 + x1*x2) ) / 
     1     ( x1*x2 * ( x1*x20 + x2*x10 )**2d0 )
*
      end function
*
************************************************************************
*
*     GA
*
************************************************************************
      double precision function ga_bis(x1,x2,x10,x20)
*
      implicit none
*
      double precision x1,x2,x10,x20
*
      ga_bis = 2d0*(x1*x2+x10*x20)*(x10**2*x20**2 + x1**2*x2**2) / 
     1     ( x1**2d0 * x2**2d0 * (x1 + x10) * ( x2 + x20 ) )
*
      end function
