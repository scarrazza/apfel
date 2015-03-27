************************************************************************
*
*     ComputeHardCrossSectionsDY.f:
*
*     Computes and dumps to file the evolution kernels in x-space for 
*     the DY kinematical points given in the input datafile.
*
************************************************************************
      subroutine ComputeHardCrossSectionsDY(datafile,outputfile)
*
      implicit none
*
      include "../commons/mxdata.h"
      include "../commons/kinematics.h"
**
*     Input variables
*
      character*(*) datafile
      character*(*) outputfile
**
*     Internal Variables
*
      integer i
      double precision t1,t2
*
*     Read kinematics of data points
*
      call ReadDataFile(datafile)
*
*     Initialization of the interpolation grid
*
      call initxGridDY
*
      call cpu_time(t1)
*
      write(6,*) "   "
      write(6,*) "Computing DY hard cross sections on the grid ..."
      write(6,*) "   "
*
      do i=1,ndata
*     report Kinematics
         write(6,*) "Point number =",i," /",ndata
         write(6,*) "Obsevable    = ",obs(i)
         write(6,*) "Kinematics:"
         write(6,*) "   Vs =",dsqrt(Q2dat(i)/x1dat(i)/x2dat(i))," GeV"
         write(6,*) "   Q2 =",Q2dat(i)," GeV^2"
         write(6,*) "   y  =",ydat(i)
         write(6,*) "   x1 =",x1dat(i)
         write(6,*) "   x2 =",x2dat(i)
         write(6,*) "   "
*
         call wcoeffDY(i)
      enddo
*
      call writeCDY(outputfile)
*
      call cpu_time(t2)
*
      write(6,*) "DY hard cross section generated in",t2-t1," s"
      write(6,*) "   "
*
      end
