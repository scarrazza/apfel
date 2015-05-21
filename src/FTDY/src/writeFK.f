************************************************************************
*
*     writeFK.f:
*     
*     It writes to file the FK tables.
*
************************************************************************
      subroutine writeFK(idat,ounit,flmap)
*
      implicit none
*
      include "../commons/mxdata.h"
      include "../commons/mxgridsize.h"
      include "../commons/sigmafk.h"
      include "../commons/xgridDY.h"
      include "../commons/kinematics.h"
      include "../commons/set.h"
      include "../commons/cutoff.h"
**
*     Input Variables
*
      integer idat,ounit
**
*     Internal Variebles
*

      integer nIntervals,nx
      integer i
      integer ipdf,jpdf,jx,kx
      integer jpdf1,jpdf2,kx1,kx2
      integer ixp(2)
      integer length(0:13,0:13)
      double precision x0(2)
      double precision dist(0:mxgridsize)
      double precision xGrid,xg(0:mxgridsize)
      character*22 buffer(0:13,0:13)
**
*     Output Variables
*
      integer flmap(0:13,0:13)
*
      nx = nIntervals()
      do jx=0,nx
         xg(jx) = xGrid(jx)
      enddo
*
*     find ixp(1) and ixp(2)
*
      x0(1) = x1dat(idat)
      x0(2) = x2dat(idat)
      do i=1,2
         ixp(i)  = 0
         dist(0) = x0(i) - xg(0)
         do kx=1,nx
            dist(kx) = x0(i) - xg(kx)
            if(dist(kx)*dist(kx-1).lt.0d0) goto 101
         enddo
 101     ixp(i) = kx - 1
         if(ixp(i).lt.1)then
            write(6,*) "The x-space grid is not wide enough."
            write(6,*) "Decrease the lower bound."
            write(6,*) "   "
            call exit(-10)
         endif
      enddo
*
*     Determine flavour map
*
      do jpdf1=0,13
         do jpdf2=0,13
            flmap(jpdf1,jpdf2) = 0
         enddo
      enddo
      do kx1=ixp(1)-1,nx-1
         do kx2=ixp(2)-1,nx-1
            do jpdf1=0,13
               do jpdf2=0,13
                  if(abs(sigmafkdy(kx1,kx2,jpdf1,jpdf2)).gt.cutoff) 
     1                 flmap(jpdf1,jpdf2) = 1
               enddo
            enddo
         enddo
      enddo
*
      do kx1=ixp(1)-1,nx-1
         do kx2=ixp(2)-1,nx-1
            do ipdf=0,13
               do jpdf=0,13
                  if(flmap(ipdf,jpdf).eq.1)then 
                     length(ipdf,jpdf) = 22
                     write(buffer(ipdf,jpdf),*) 
     1                    sigmafkdy(kx1,kx2,ipdf,jpdf)
                  else
                     length(ipdf,jpdf) = 3
                     buffer(ipdf,jpdf) = " 0 "
                  endif
               enddo
            enddo
            write(ounit,*) idat,kx1,kx2,
     1           ((buffer(ipdf,jpdf)(1:length(ipdf,jpdf)),
     2           jpdf=0,13),ipdf=0,13)
         enddo
      enddo
*
      return
      end
