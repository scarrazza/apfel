************************************************************************
*
*     readcDY.f:
*
*     It reads to disk the the Drell-Yan hard cross section
*
************************************************************************
      subroutine readcDY(inputfile)
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
      character*50 inputfile
**
*     Internal Variables
*
      integer ix,jx,ip,ipt
      integer idum
      integer ls
*
      ls = index(inputfile,".hcx") + 4
      open(unit=10,status="unknown",file=inputfile(1:ls))
*
*     Read name of the set
*
      read(10,*) set
*
*     Read grid
*
      read(10,*) nxDY
      do ix=1,nxDY
         read(10,*) xxDY(ix)
      enddo
*
*     Read Kernels
*
      read(10,*) ndata
      do ip=1,ndata
         read(10,*) idum,obs(ip),x1dat(ip),x2dat(ip),q2dat(ip),
     1              ixp1DY(ip),ixp2DY(ip)
         do ipt=0,1
            do ix=ixp1DY(ip),nxDY
               read(10,*) idum,(cDY_NS(ipt,jx,ix,ip),jx=ixp2DY(ip),nxDY)
            enddo
         enddo
*
*     No need to read the LO for QG and GQ because it's identically zero
*
         do ix=ixp1DY(ip),nxDY
            read(10,*) idum,(cDY_QG(1,jx,ix,ip),jx=ixp2DY(ip),nxDY)
         enddo
         do ix=ixp1DY(ip),nxDY
            read(10,*) idum,(cDY_GQ(1,jx,ix,ip),jx=ixp2DY(ip),nxDY)
         enddo
      enddo
*
      close(10)
*
      return
      end
