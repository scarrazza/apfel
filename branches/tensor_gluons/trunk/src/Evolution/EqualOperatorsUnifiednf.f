************************************************************************
*
*     EqualOperatorsUnifiednf.f
*
*     This routine equals the first six entries with the second six 
*     where the latter have a further index which is the number of active
*     flavours nf.
*
************************************************************************
      subroutine EqualOperatorsUnifiednf(nf,M0sg1,M0sg2,M0nspu,M0nspd,
     1                                                  M0nsmu,M0nsmd,
     2                                      Msg1,Msg2,Mnspu,Mnspd,
     3                                                Mnsmu,Mnsmd)

*
      implicit none
*
      include "../commons/grid.h"
**
*     Input Variables
*
      integer nf
      double precision M0sg1(4,4,0:nint_max,0:nint_max)
      double precision M0sg2(2,2,0:nint_max,0:nint_max)
      double precision M0nspu(0:nint_max,0:nint_max)
      double precision M0nspd(0:nint_max,0:nint_max)
      double precision M0nsmu(0:nint_max,0:nint_max)
      double precision M0nsmd(0:nint_max,0:nint_max)
**
*     Internal Variables
*
      integer i,j
      integer alpha,beta
**
*     Output Variables
*
      double precision Msg1(3:6,4,4,0:nint_max,0:nint_max)
      double precision Msg2(3:6,2,2,0:nint_max,0:nint_max)
      double precision Mnspu(3:6,0:nint_max,0:nint_max)
      double precision Mnspd(3:6,0:nint_max,0:nint_max)
      double precision Mnsmu(3:6,0:nint_max,0:nint_max)
      double precision Mnsmd(3:6,0:nint_max,0:nint_max)
*
      do alpha=0,nin(igrid)
         do beta=0,nin(igrid)
            do i=1,4
               do j=1,4
                  Msg1(nf,i,j,alpha,beta) = M0sg1(i,j,alpha,beta)
               enddo
            enddo
*
            do i=1,2
               do j=1,2
                  Msg2(nf,i,j,alpha,beta) = M0sg2(i,j,alpha,beta)
               enddo
            enddo
*
            Mnspu(nf,alpha,beta) = M0nspu(alpha,beta)
            Mnspd(nf,alpha,beta) = M0nspd(alpha,beta)
            Mnsmu(nf,alpha,beta) = M0nsmu(alpha,beta)
            Mnsmd(nf,alpha,beta) = M0nsmd(alpha,beta)
         enddo
      enddo
*
      return
      end
