************************************************************************
*
*     GetCKM.f:
*
*     This function returns the absolute values of the entries of the 
*     CKM matrix.
*
************************************************************************
      function GetCKM(u,d)
*
      implicit none
*
      include "../commons/CKM.h"
*
*     Variables
*
      integer u,d
      double precision GetCKM
*
      if(inCKM.ne."done")then
         write(6,*) "GetCKM: Parameter not initialized"
         write(6,*) "Set it by means of 'SetCKM'"
         write(6,*) "  "
         call exit(-10)
      endif
*
      if(u.lt.1.or.u.gt.3)then
         write(6,*) "invalid value of the up index"
         write(6,*) "u =",u
         call exit(-10)
      endif
*
      if(d.lt.1.or.d.gt.3)then
         write(6,*) "invalid value of the down index"
         write(6,*) "d =",u
         call exit(-10)
      endif
*
      if(u.eq.1.and.d.eq.1) GetCKM = V_ud
      if(u.eq.2.and.d.eq.1) GetCKM = V_cd
      if(u.eq.3.and.d.eq.1) GetCKM = V_td
      if(u.eq.1.and.d.eq.2) GetCKM = V_us
      if(u.eq.2.and.d.eq.2) GetCKM = V_cs
      if(u.eq.3.and.d.eq.2) GetCKM = V_ts
      if(u.eq.1.and.d.eq.3) GetCKM = V_ub
      if(u.eq.2.and.d.eq.3) GetCKM = V_cb
      if(u.eq.3.and.d.eq.3) GetCKM = V_tb
*
      return
      end
