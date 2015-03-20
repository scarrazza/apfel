************************************************************************
*
*     EqualOperatorsUnifiednf.f
*
*     This routine equals the first six entries with the second six 
*     where the latter have a further index which is the number of active
*     flavours nf.
*
************************************************************************
      subroutine EqualOperatorsUnifiednf(nf,nl,
     1                                   M0sg1,M0sg2,M0nspu,M0nspd,
     1                                   M0nsmu,M0nsmd,M0nslep,
     3                                   Msg1,Msg2,Mnspu,Mnspd,
     4                                   Mnsmu,Mnsmd,Mnslep)
*
      implicit none
*
      include "../commons/grid.h"
**
*     Input Variables
*
      integer nf,nl
      double precision M0sg1(5,5,0:nint_max,0:nint_max)
      double precision M0sg2(2,2,0:nint_max,0:nint_max)
      double precision M0nspu(0:nint_max,0:nint_max)
      double precision M0nspd(0:nint_max,0:nint_max)
      double precision M0nsmu(0:nint_max,0:nint_max)
      double precision M0nsmd(0:nint_max,0:nint_max)
      double precision M0nslep(0:nint_max,0:nint_max)
**
*     Internal Variables
*
      integer i,j
      integer alpha,beta
**
*     Output Variables
*
      double precision Msg1(3:6,2:3,5,5,0:nint_max,0:nint_max)
      double precision Msg2(3:6,2:3,2,2,0:nint_max,0:nint_max)
      double precision Mnspu(3:6,2:3,0:nint_max,0:nint_max)
      double precision Mnspd(3:6,2:3,0:nint_max,0:nint_max)
      double precision Mnsmu(3:6,2:3,0:nint_max,0:nint_max)
      double precision Mnsmd(3:6,2:3,0:nint_max,0:nint_max)
      double precision Mnslep(2:3,0:nint_max,0:nint_max)
*
      do alpha=0,nin(igrid)
         do beta=0,nin(igrid)
            do i=1,5
               do j=1,5
                  Msg1(nf,nl,i,j,alpha,beta) = M0sg1(i,j,alpha,beta)
               enddo
            enddo
*
            do i=1,2
               do j=1,2
                  Msg2(nf,nl,i,j,alpha,beta) = M0sg2(i,j,alpha,beta)
               enddo
            enddo
*
            Mnspu(nf,nl,alpha,beta) = M0nspu(alpha,beta)
            Mnspd(nf,nl,alpha,beta) = M0nspd(alpha,beta)
            Mnsmu(nf,nl,alpha,beta) = M0nsmu(alpha,beta)
            Mnsmd(nf,nl,alpha,beta) = M0nsmd(alpha,beta)
            Mnslep(nl,alpha,beta)   = M0nsmd(alpha,beta)
         enddo
      enddo
*
      return
      end
