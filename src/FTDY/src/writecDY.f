************************************************************************
*
*     writecDY.f:
*
*     It writes to disk the the Drell-Yan hard cross section
*
************************************************************************
      subroutine writecDY(outputfile)
*
      implicit none
*
      include "../commons/mxdata.h"
      include "../commons/mxgridsizeDY.h"
      include "../commons/xgridDY.h"
      include "../commons/xxDY.h"
      include "../commons/cDY.h"
      include "../commons/kinematics.h"
      include "../commons/set.h"
**
*     Input Variables
*
      character*50 outputfile
**
*     Internal Variables
*
      integer ix,jx,ip,ipt
      integer ls
*
      ls = index(outputfile,".hcx") + 4
      open(unit=10,status="unknown",file=outputfile(1:ls))
*
*     Write name of the set
*
      write(10,*) set
*
*     Write grid
*
      write(10,*) nxDY
      do ix=1,nxDY
         write(10,*) xxDY(ix)
      enddo
*
*     Write Kernels
*
      write(10,*) ndata
      do ip=1,ndata
         write(10,*) ip,obs(ip),x1dat(ip),x2dat(ip),q2dat(ip),
     1               ixp1DY(ip),ixp2DY(ip)
         do ipt=0,1
            do ix=ixp1DY(ip),nxDY
               write(10,*) ipt,(cDY_NS(ipt,jx,ix,ip),jx=ixp2DY(ip),nxDY)
            enddo
         enddo
*
*     No need to write the LO for QG and GQ because it's identically zero
*
         do ix=ixp1DY(ip),nxDY
            write(10,*) ipt,(cDY_QG(1,jx,ix,ip),jx=ixp2DY(ip),nxDY)
         enddo
         do ix=ixp1DY(ip),nxDY
            write(10,*) ipt,(cDY_GQ(1,jx,ix,ip),jx=ixp2DY(ip),nxDY)
         enddo
      enddo
*
      close(10)
*
      return
      end
