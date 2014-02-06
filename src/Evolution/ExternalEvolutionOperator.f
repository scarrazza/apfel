************************************************************************
*
*     ExternalEvolutionOperator.f:
*
*     This subroutine computes the "external" evolution operator on a
*     user-given external grid starting from the "internal" evolution
*     operators on the internal x-space grid of APFEL.
*
*     The evolution operator "M" on the external grid "xint" is given
*     as one-dimensional array in order to be passed to the c++ wrapper.
*
************************************************************************
      subroutine ExternalEvolutionOperator(Q0,Q,n,xint,M)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/EvolutionOperator.h"
      include "../commons/EvolutionMatrices.h"
**
*     Input Variables
*
      double precision Q0,Q
**
*     Internal Variables
*
      integer k
      integer i,j
      integer alpha,beta
      double precision eps
      parameter(eps=1d-12)

      include "../commons/fph.h"
      double precision f0(-6:6,0:n),xpd(-6:6)
**
*     Output Variables
*
      integer n
      double precision xint(0:200)
      double precision M(0:14*14*(200+1)*(200+1)-1)
*
*     Disable welcome message
*
      call EnableWelcomeMessage(.false.)
*
*     Enable computation of the Total Evolution operator
*
      call EnableEvolutionOperator(.true.)
*
*     Initializes integrals on the grids
*
      call InitializeAPFEL
*
*     Compute APFEL evolution operators
*
      call EvolveAPFEL(Q0,Q)
*
*     Set grid
*
      n = nin(0)
      do alpha=0,n
         xint(alpha) = xg(0,alpha)
      enddo
*
*     Compute Evolution Operator
*
      do alpha=0,n
        do beta=0,n
            do i=-6,6
               do j=-6,6
                  k = ( 7 + i ) + 14 * ( ( 7 + j ) 
     1                 + 14 * ( alpha + ( n + 1 ) * beta ) )
                  M(k) = 0d0
               enddo
            enddo
         enddo
      enddo
*
      do alpha=0,n
         do igrid=1,ngrid
            if(xint(alpha).ge.xmin(igrid).and.
     1         xint(alpha).lt.xmin(igrid+1))then
               goto 101
            endif
         enddo
 101     do beta=alpha,n
            do i=-nff,nff
               do j=-nfi,nfi
                  k = ( 7 + i ) + 14 * ( ( 7 + j ) 
     1                 + 14 * ( alpha + ( n + 1 ) * beta ) )
                  M(k) = PhQCD(igrid,i,j,alpha,beta) 
                  if(dabs(M(k)).lt.eps) M(k) = 0d0
               enddo
            enddo
         enddo
      enddo
*
      do alpha=0,n
         call toyLHPDFs(xint(alpha),xpd)
         do i=-6,6
            f0(i,alpha) = xpd(i)
         enddo
      enddo
*
      do alpha=0,n
         do i=-6,6
            fph(0,i,alpha) = 0d0
            do beta=0,n
               do j=-6,6
                  k = ( 7 + i ) + 14 * ( ( 7 + j ) 
     1              + 14 * ( alpha + ( n + 1 ) * beta ) )
                  fph(0,i,alpha) = fph(0,i,alpha) + M(k) * f0(j,beta)
               enddo
            enddo
         enddo
      enddo
c*
c      do alpha=0,n
c         do beta=0,n
c            write(58,*) alpha,beta
c            do i=1,13
c               write(58,"(13(2x,f8.5))")
c     1              (M( i + 14 * ( j 
c     2              + 14 * ( alpha + ( n + 1 ) * beta ) )),j=1,13)
c            enddo
c         enddo
c      enddo
*
      return
      end
