************************************************************************
*
*     SetTargetDIS.f:
*
*     This subroutine sets the DIS process.
*
************************************************************************
      subroutine SetTargetDIS(tar)
*
      implicit none
*
      include "../commons/TargetDIS.h"
*
*     Variables
*
      character*9 tar
*
      if(tar(1:6).eq."proton")then
         TargetDIS = tar(1:6)
      elseif(tar(1:7).eq."neutron")then
         TargetDIS = tar(1:7)
      elseif(tar.eq."isoscalar")then
         TargetDIS = tar
      endif
      InTargetDIS = "done"
*
      return
      end
