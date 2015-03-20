************************************************************************
*
*     IdentityOperatorsUnified.f
*
*     This routine returns the identity operators for singlet and 
*     non-singlet used as bonduary conditions for the solution of the
*     DGLAP equation in the unified basis.
*
************************************************************************
      subroutine IdentityOperatorsUnified(M0sg1,M0sg2,M0nspu,M0nspd,
     1                                                M0nsmu,M0nsmd)
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
      double precision M0sg1(5,5,0:nint_max,0:nint_max)
      double precision M0sg2(2,2,0:nint_max,0:nint_max)
      double precision M0nspu(0:nint_max,0:nint_max)
      double precision M0nspd(0:nint_max,0:nint_max)
      double precision M0nsmu(0:nint_max,0:nint_max)
      double precision M0nsmd(0:nint_max,0:nint_max)
*
      do alpha=0,nin(igrid)
         do beta=0,nin(igrid)
            do i=1,5
               do j=1,5
                  if(alpha.eq.beta.and.i.eq.j)then
                     M0sg1(i,j,alpha,beta) = 1d0
                  else
                     M0sg1(i,j,alpha,beta) = 0d0
                  endif
               enddo
            enddo
            do i=1,2
               do j=1,2
                  if(alpha.eq.beta.and.i.eq.j)then
                     M0sg2(i,j,alpha,beta) = 1d0
                  else
                     M0sg2(i,j,alpha,beta) = 0d0
                  endif
               enddo
            enddo
*
            if(alpha.eq.beta)then
               M0nspu(alpha,beta) = 1d0
               M0nspd(alpha,beta) = 1d0
               M0nsmu(alpha,beta) = 1d0
               M0nsmd(alpha,beta) = 1d0
            else
               M0nspu(alpha,beta) = 0d0
               M0nspd(alpha,beta) = 0d0
               M0nsmu(alpha,beta) = 0d0
               M0nsmd(alpha,beta) = 0d0
            endif
         enddo
      enddo
*
      return
      end
