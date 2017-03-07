************************************************************************
*
*     GetThreshold.f:
*
*     This function returns the value of i-th (i=4,5,6) heavy quark
*     threshold.
*
************************************************************************
      function GetThreshold(i)
*
      implicit none
*
      include "../commons/m2th.h"
**
*     Input Variables
*
      integer i
**
*     Output Variables
*
      double precision GetThreshold
*
*     Check the consistency of the input
*
      if(i.lt.4.and.i.gt.6)then
         write(6,*) "In GetThreshold.f:"
         write(6,*) "Invalid heavy quark index i =",i
         call exit(-10)
      endif
*
      GetThreshold = dsqrt(m2th(i))
*
      return
      end
*
************************************************************************
      function HeavyQuarkThreshold(i)
*
      implicit none
*
      include "../commons/m2th.h"
**
*     Input Variables
*
      integer i
**
*     Output Variables
*
      double precision HeavyQuarkThreshold
*
*     Check the consistency of the input
*
      if(i.lt.4.and.i.gt.6)then
         write(6,*) "In GetThreshold.f:"
         write(6,*) "Invalid heavy quark index i =",i
         call exit(-10)
      endif
*
      HeavyQuarkThreshold = dsqrt(m2th(i))
*
      return
      end
