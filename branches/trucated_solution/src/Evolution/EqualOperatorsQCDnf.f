************************************************************************
*
*     EqualOperatorsQCDnf.f
*
*     This routine equals the first four entries with the second four 
*     where the latter have a further index which is the number of active
*     flavours nf.
*
************************************************************************
      subroutine EqualOperatorsQCDnf(nf,M0sg,M0nsp,M0nsm,M0nsv,
     1                                  Msg,Mnsp,Mnsm,Mnsv)
*
      implicit none
*
      include "../commons/grid.h"
**
*     Input Variables
*
      integer nf
      double precision M0sg(2,2,0:nint_max,0:nint_max)
      double precision M0nsp(0:nint_max,0:nint_max)
      double precision M0nsm(0:nint_max,0:nint_max)
      double precision M0nsv(0:nint_max,0:nint_max)
**
*     Internal Variables
*
      integer i,j
      integer alpha,beta
**
*     Output Variables
*
      double precision Msg(3:6,2,2,0:nint_max,0:nint_max)
      double precision Mnsp(3:6,0:nint_max,0:nint_max)
      double precision Mnsm(3:6,0:nint_max,0:nint_max)
      double precision Mnsv(3:6,0:nint_max,0:nint_max)
*
      do alpha=0,nin(igrid)
         do beta=0,nin(igrid)
            Mnsp(nf,alpha,beta) = M0nsp(alpha,beta)
            Mnsm(nf,alpha,beta) = M0nsm(alpha,beta)
            Mnsv(nf,alpha,beta) = M0nsv(alpha,beta)
*
            do i=1,2
               do j=1,2
                  Msg(nf,i,j,alpha,beta) = M0sg(i,j,alpha,beta)
               enddo
            enddo
         enddo
      enddo
*
      return
      end
