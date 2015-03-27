************************************************************************
*
*     ReadDataFile.f:
*
*     Read the kinematics from the input datafile.
*
************************************************************************
      subroutine ReadDataFile(datafile)
*
      implicit none
*
      include "../commons/mxdata.h"
      include "../commons/xgridDY.h"
      include "../commons/kinematics.h"
      include "../commons/set.h"
**
*     Input Variables
*
      character*100 datafile
**
*     Internal Variables
*
      integer i,idum
      integer ls
      double precision ss,stau
*
*     Read kinematical points from the kinematics.in file
*
      ls = index(datafile,".dat") + 4
      open(unit=20,status="unknown",file=datafile(1:ls))
      read(20,*) set,idum,ndata
*
*     Check that the number of data points in kinematics.in is not too big
*
      if(ndata.gt.mxdata)then
         write(6,*) "ERROR: in ReadDataFile.f:"
         write(6,*) "Number of data points too large,"
         write(6,*) "ndata =",ndata,", mxdata =",mxdata
         call exit(-10)
      endif
*
      nxDY   = 100
      xminDY = 1d0
      xmaxDY = 1d0
*
      do i=1,ndata
         read(20,*) idum,obs(i),ydat(i),q2dat(i),ss
*
         stau = dsqrt(q2dat(i)) / ss
         x1dat(i) = stau * dexp(ydat(i))
         x2dat(i) = stau * dexp(-ydat(i))
*
         if(x1dat(i).le.xminDY) xminDY = x1dat(i)
         if(x2dat(i).le.xminDY) xminDY = x2dat(i)
      enddo
      close(20)
*
      return 
      end
