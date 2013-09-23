************************************************************************
*
*     EqualOperatorsQEDnf.f
*
*     This routine equals the fisrt and the second entries with the
*     third and the forth respectively where the latter have a further
*     index which is the number of active flavours. 
*     The first and the third entries are supposed to be QED non-singlets 
*     while the secons and the forth QED singlets.
*
************************************************************************
      subroutine EqualOperatorsQEDnf(nf,M0sg,M0nsp,M0nsm,Msg,Mnsp,Mnsm)
*
      implicit none
*
      include "../commons/grid.h"
**
*     Input Variables
*
      integer nf
      double precision M0sg(3,3,0:nint_max,0:nint_max)
      double precision M0nsp(0:nint_max,0:nint_max)
      double precision M0nsm(0:nint_max,0:nint_max)
**
*     Internal Variables
*
      integer i,j
      integer alpha,beta
**
*     Output Variables
*
      double precision Msg(3:6,3,3,0:nint_max,0:nint_max)
      double precision Mnsp(3:6,0:nint_max,0:nint_max)
      double precision Mnsm(3:6,0:nint_max,0:nint_max)
*
      do alpha=0,nin(igrid)
         do beta=0,nin(igrid)
            Mnsp(nf,alpha,beta) = M0nsp(alpha,beta)
            Mnsm(nf,alpha,beta) = M0nsm(alpha,beta)
*
            do i=1,3
               do j=1,3
                  Msg(nf,i,j,alpha,beta) = M0sg(i,j,alpha,beta)
               enddo
            enddo
         enddo
      enddo
*
      return
      end
