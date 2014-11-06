************************************************************************
*
*     IdentityOperatorsQED.f
*
*     This routine returns the identity operators for singlet and 
*     non-singlet used as bonduary conditions for the solution of the
*     DGLAP equation in QED.
*
************************************************************************
      subroutine IdentityOperatorsQED(M0sg,M0nsp,M0nsm)
*
      implicit none
*
      include "../commons/grid.h"
**
*     Internal Variables
*
      integer i,j
      integer alpha,beta
**
*     Output Variables
*
      double precision M0sg(3,3,0:nint_max,0:nint_max)
      double precision M0nsp(0:nint_max,0:nint_max)
      double precision M0nsm(0:nint_max,0:nint_max)
*
      do alpha=0,nin(igrid)
         do beta=0,nin(igrid)
            do i=1,3
               do j=1,3
                  if(alpha.eq.beta.and.i.eq.j)then
                     M0sg(i,j,alpha,beta) = 1d0
                  else
                     M0sg(i,j,alpha,beta) = 0d0
                  endif
               enddo
            enddo
*
            if(alpha.eq.beta)then
               M0nsp(alpha,beta) = 1d0
               M0nsm(alpha,beta) = 1d0
            else
               M0nsp(alpha,beta) = 0d0
               M0nsm(alpha,beta) = 0d0
            endif
         enddo
      enddo
*
      return
      end
