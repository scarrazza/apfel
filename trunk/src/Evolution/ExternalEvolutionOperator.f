************************************************************************
*
*     ExternalEvolutionOperator.f:
*
*     This subroutine computes the "external" evolution operator on a
*     user-given external grid starting from the "internal" evolution
*     operators on the internal x-space grid of APFEL.
*
*     The evolution operator "M" on the external grid "xext" is given
*     as one-dimensional array in order to be passed to the c++ wrapper.
*
************************************************************************
      subroutine ExternalEvolutionOperator(Q0,Q,n,xext,M)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/EvolutionOperator.h"
      include "../commons/EvolutionMatrices.h"
**
*     Input Variables
*
      integer n
      double precision xext(0:n)
      double precision Q0,Q
**
*     Internal Variables
*
      integer k
      integer i,j
      integer alpha,beta,alphap,betap
      double precision w_int_herm,xint(0:200)
      double precision inta(0:nint_max,0:n)
      double precision intb(0:n,0:nint_max)
      double precision eps
      parameter(eps=1d-15)
**
*     Output Variables
*
      double precision M(0:14*14*(n+1)*(n+1)-1)
*
*     Set Evolution Operator to zero
*
      do alphap=0,n
         do betap=0,n
            do i=-6,6
               do j=-6,6
                  k = ( 7 + i ) + 14 * ( ( 7 + j ) 
     1              + 14 * ( alphap + ( n + 1 ) * betap ) )
                  M(k) = 0d0
                  if(Q.eq.Q0.and.betap.eq.alphap.and.j.eq.i) M(k) = 1d0
               enddo
            enddo
         enddo
      enddo
      if(Q.eq.Q0) return
*
*     Disable welcome message
*
      call EnableWelcomeMessage(.false.)
*
*     Enable computation of the Total Evolution operator
*     (it also locks internal subgrids)
*
      call EnableEvolutionOperator(.true.)
*
*     Tune internal grids
*
      call SetNumberOfGrids(3)
      call SetGridParameters(1,80,3,xext(0))
      call SetGridParameters(2,50,5,1d-1)
      call SetGridParameters(3,40,5,8d-1)
*
*     Initializes integrals on the grids
*
      call InitializeAPFEL
*
*     Compute APFEL evolution operators
*
      call EvolveAPFEL(Q0,Q)
*
*     Interpolation functions on the joint grid
*
      do alpha=0,nin(0)
         xint(alpha) = xg(0,alpha)
      enddo
      do alphap=0,n
         do alpha=0,nin(0)
            inta(alpha,alphap) = 
     1           w_int_herm(nin(0),xint,alpha,xext(alphap))
            intb(alphap,alpha) = 
     1           w_int_herm(n,xext,alphap,xint(alpha))
         enddo
      enddo
*
*     Fill Evolution Operator
*
      do i=-nff,nff
         do j=-nfi,nfi
            do alphap=0,n
               do betap=alphap,n
                  k = ( 7 + i ) + 14 * ( ( 7 + j ) 
     1                 + 14 * ( alphap + ( n + 1 ) * betap ) )
                  do alpha=0,nin(0)
                     do beta=alpha,nin(0)
                        M(k) = M(k)
     1                       + inta(alpha,alphap) 
     2                       * PhQCD(0,i,j,alpha,beta) 
     3                       * intb(betap,beta)
                     enddo
                  enddo
                  if(dabs(M(k)).lt.eps) M(k) = 0d0
               enddo
            enddo
         enddo
      enddo
*
      return
      end
