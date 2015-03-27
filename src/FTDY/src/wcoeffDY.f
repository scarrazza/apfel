************************************************************************
*
*     wcoeffDY.f:
*
*     Compute LO and NLO coefficient functions for DY processes according 
*     to a fastNLO-like formalism. The interpolating basis is a linear one.
*
************************************************************************
      subroutine wcoeffDY(idat)
*
      implicit none
*
      include "../commons/mxdata.h"
      include "../commons/mxgridsizeDY.h"
      include "../commons/cDY.h"
      include "../commons/xxDY.h"
      include "../commons/xgridDY.h"
      include "../commons/kinematics.h"
      include "../commons/ngau.h"
      include "../commons/ngaus.h"
**
*     Input Variables
*
      integer idat
**
*     Internal Variables
*
      integer i
      integer ix,jx
      integer ixp(2)
      integer ipt
      double precision x0(2)
      double precision dist(2),mdist(2)
      double precision I0Q,I1Q,I2Q,I3Q
      double precision I1G,I2G,I3QG,I3GQ
      double precision elin,E20,E10
      double precision ZeroIntegrand_QQ,SUB1,SUB2
      double precision SingleIntegrand_sub_X1,SingleIntegrand_sub_X2
      double precision SingleIntegrand_X1_QQ,SingleIntegrand_X2_QQ
      double precision DoubleIntegrand_QQ,DoubleIntegrand_sub1_QQ
      double precision DoubleIntegrand_sub2_QQ,DoubleIntegrand_sub3_QQ
      double precision SingleIntegrand_X1_GQ,SingleIntegrand_X2_QG
      double precision DoubleIntegrand_QG,DoubleIntegrand_sub_QG
      double precision DoubleIntegrand_GQ,DoubleIntegrand_sub_GQ
      double precision cDY_NS_LO(0:1,0:1)
*
      common/cixp/ixp
*
*     Initialization of the DY coefficient functions
*
      do ipt=0,1
         do jx=1,nxDY
            do ix=1,nxDY
               cDY_NS(ipt,jx,ix,idat) = 0d0
               cDY_QG(ipt,jx,ix,idat) = 0d0
               cDY_GQ(ipt,jx,ix,idat) = 0d0
            enddo
         enddo
      enddo
*
      x0(1) = x1dat(idat)
      x0(2) = x2dat(idat)
*
*     find ixp(1) and ixp(2) such that xxDY(ixp) < x0(i) < xxDY(ixp+1)
*
      do i=1,2
         ixp(i)   = 1
         mdist(i) = x0(i) - xxDY(1)
         do ix=1,nxDY-1
            dist(i) = x0(i) - xxDY(ix)
            if(dist(i).ge.(-1d-9).and.dist(i).le.mdist(i))then
               ixp(i) = ix
               mdist(i) = x0(i) - xxDY(ix)
            endif
         enddo
      enddo
*
      ixp1DY(idat) = ixp(1)
      ixp2DY(idat) = ixp(2)
*
*     LO Coefficient functions
*
      call cDY_LO(ixp(1),ixp(2),x0,cDY_NS_LO)
*
      cDY_NS(0,ixp(2),ixp(1),idat)     = cDY_NS_LO(0,0)
      cDY_NS(0,ixp(2),ixp(1)+1,idat)   = cDY_NS_LO(0,1)
      cDY_NS(0,ixp(2)+1,ixp(1),idat)   = cDY_NS_LO(1,0)
      cDY_NS(0,ixp(2)+1,ixp(1)+1,idat) = cDY_NS_LO(1,1)
*
*     NLO Coefficient functions
*
      do ix=ixp(1),nxDY
         do jx=ixp(2),nxDY
            I0Q  = 0d0
            I1Q  = 0d0
            I1G  = 0d0
            I2Q  = 0d0
            I2G  = 0d0
            I3Q  = 0d0
            I3QG = 0d0
            I3GQ = 0d0
*     
            if(jx.le.ixp(2)+1.and.ix.le.ixp(1)+1)then
               E10 = elin(ix,x0(1))
               E20 = elin(jx,x0(2))    
*     
               call gauleg(x0(1),xxDY(ix+1),Y1,W1,ngau) 
               call gauleg(x0(2),xxDY(jx+1),Y2,W2,ngau)
*     
               I0Q  = E10 * E20 * ZeroIntegrand_QQ(x0)                 
               I1Q  = E10 * SingleIntegrand_X2_QQ(jx,x0)   
               I1G  = E10 * SingleIntegrand_X2_QG(jx,x0)    
               I2Q  = E20 * SingleIntegrand_X1_QQ(ix,x0)
               I2G  = E20 * SingleIntegrand_X1_GQ(ix,x0)             
               I3Q  = DoubleIntegrand_QQ(ix,jx,x0)      
               I3QG = DoubleIntegrand_QG(ix,jx,x0)
               I3GQ = DoubleIntegrand_GQ(ix,jx,x0)
*     
               call gauleg(xxDY(ix+1),1d0,Y1S,W1S,ngaus) 
               call gauleg(xxDY(jx+1),1d0,Y2S,W2S,ngaus)
*     
               SUB1 = E10 * E20 * SingleIntegrand_sub_X2(x0)
               I1Q  = I1Q - SUB1
               SUB2 = E20 * E10 * SingleIntegrand_sub_X1(x0)  
               I2Q  = I2Q - SUB2
               I3Q  = I3Q + E10 * E20 * DoubleIntegrand_sub1_QQ(x0)
*     
               call gauleg(x0(1),xxDY(ix+1),Y1S,W1S,ngaus)
               call gauleg(xxDY(jx+1),1d0,Y2S,W2S,ngaus)
*     
               I3Q  = I3Q  - E20 * DoubleIntegrand_sub2_QQ(ix,x0)
               I3GQ = I3GQ - E20 * DoubleIntegrand_sub_GQ(ix,x0)
*     
               call gauleg(xxDY(ix+1),1d0,Y1S,W1S,ngaus) 
               call gauleg(x0(2),xxDY(jx+1),Y2S,W2S,ngaus)
*     
               I3Q  = I3Q  - E10 * DoubleIntegrand_sub3_QQ(jx,x0)
               I3QG = I3QG - E10 * DoubleIntegrand_sub_QG(jx,x0)
            elseif(jx.ge.ixp(2)+2.and.ix.le.ixp(1)+1)then 
               E10  = elin(ix,x0(1))
               call gauleg(x0(1),xxDY(ix+1),Y1,W1,ngau)
               if(jx.eq.nxDY)then
                  call gauleg(xxDY(jx-1),xxDY(jx),Y2,W2,ngau)
               else
                  call gauleg(xxDY(jx-1),xxDY(jx+1),Y2,W2,ngau)
               endif
*     
               I1Q  = E10 * SingleIntegrand_X2_QQ(jx,x0)
               I1G  = E10 * SingleIntegrand_X2_QG(jx,x0)
               I3Q  = DoubleIntegrand_QQ(ix,jx,x0)
               I3QG = DoubleIntegrand_QG(ix,jx,x0)
               I3GQ = DoubleIntegrand_GQ(ix,jx,x0)
*     
               call gauleg(xxDY(ix+1),1d0,Y1S,W1S,ngaus)
               if(jx.eq.nxDY)then
                  call gauleg(xxDY(jx-1),xxDY(jx),Y2S,W2S,ngaus)
               else
                  call gauleg(xxDY(jx-1),xxDY(jx+1),Y2S,W2S,ngaus)
               endif
*     
               I3Q  = I3Q  - E10 * DoubleIntegrand_sub3_QQ(jx,x0)
               I3QG = I3QG - E10 * DoubleIntegrand_sub_QG(jx,x0)
            elseif(jx.le.ixp(2)+1.and.ix.ge.ixp(1)+2)then
               E20 = elin(jx,x0(2))
               if(ix.eq.nxDY)then
                  call gauleg(xxDY(ix-1),xxDY(ix),Y1,W1,ngau)
               else
                  call gauleg(xxDY(ix-1),xxDY(ix+1),Y1,W1,ngau)
               endif
               call gauleg(x0(2),xxDY(jx+1),Y2,W2,ngau)
*     
               I2Q  = E20 * SingleIntegrand_X1_QQ(ix,x0)
               I2G  = E20 * SingleIntegrand_X1_GQ(ix,x0)
               I3Q  = DoubleIntegrand_QQ(ix,jx,x0)               
               I3QG = DoubleIntegrand_QG(ix,jx,x0)
               I3GQ = DoubleIntegrand_GQ(ix,jx,x0)
*     
               if(ix.eq.nxDY)then
                  call gauleg(xxDY(ix-1),xxDY(ix),Y1S,W1S,ngaus)
               else
                  call gauleg(xxDY(ix-1),xxDY(ix+1),Y1S,W1S,ngaus)
               endif
               call gauleg(xxDY(jx+1),1d0,Y2S,W2S,ngaus) 
*     
               I3Q  = I3Q  - E20 * DoubleIntegrand_sub2_QQ(ix,x0)
               I3GQ = I3GQ - E20 * DoubleIntegrand_sub_GQ(ix,x0)
            elseif(jx.ge.ixp(2)+2.and.ix.ge.ixp(1)+2)then
               if(ix.eq.nxDY)then
                  call gauleg(xxDY(ix-1),xxDY(ix),Y1,W1,ngau)
               else
                  call gauleg(xxDY(ix-1),xxDY(ix+1),Y1,W1,ngau)
               endif 
*     
               if(jx.eq.nxDY)then
                  call gauleg(xxDY(jx-1),xxDY(jx),Y2,W2,ngau)
               else
                  call gauleg(xxDY(jx-1),xxDY(jx+1),Y2,W2,ngau)
               endif
*     
               I3Q  = DoubleIntegrand_QQ(ix,jx,x0)
               I3QG = DoubleIntegrand_QG(ix,jx,x0)
               I3GQ = DoubleIntegrand_GQ(ix,jx,x0)
            endif
*     
            cDY_NS(1,jx,ix,idat) = I0Q + I1Q + I2Q + I3Q
            cDY_QG(1,jx,ix,idat) = I1G + I3QG
            cDY_GQ(1,jx,ix,idat) = I2G + I3GQ
         enddo
      enddo
*
      return
      end
