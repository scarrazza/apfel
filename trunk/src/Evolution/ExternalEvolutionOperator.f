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
      double precision Q0,Q
      double precision xext(0:n)
**
*     Internal Variables
*
      integer k
      integer i,j
      integer alpha,beta,alphap,betap
      double precision w_int,w_int_ext
      double precision inta(ngrid_max,0:nint_max,0:n)
      double precision intb(ngrid_max,0:n,0:nint_max)
      double precision eps
      parameter(eps=1d-12)

      include "../commons/fph.h"
      double precision f0(-6:6,0:n),xpd(-6:6)
**
*     Output Variables
*
      double precision M(0:14*14*(n+1)*(n+1)-1)
*
*     Disable welcome message
*
      call EnableWelcomeMessage(.false.)
*
*     Enable computation of the Total Evolution operator
*
      call EnableEvolutionOperator(.true.)
*
*     Checks
*
      if(dabs(xext(n)-1d0).gt.eps)then
         write(6,*) "In src/Evolution/ExternalEvolutionOperator.f:"
         write(6,*) "The upper limit of the external grid must be"
         write(6,*) "equal to 1, while n, xext(n) =",n,xext(n)
         call exit(-10)
      endif
*
      if(xext(0).lt.1d-7)then
         write(6,*) "In src/Evolution/ExternalEvolutionOperator.f:"
         write(6,*) "Lower limit of the external grid must be"
         write(6,*) "bigger than 10^-7, while xext(0) =",xext(0)
         call exit(-10)
      endif
*
*     Tuning of the internal grid
*
      if(xext(0).le.1d-5)then
         call SetNumberOfGrids(3)
         call SetGridParameters(1,80,3,xext(0))
         call SetGridParameters(2,50,5,1d-1)
         call SetGridParameters(3,40,5,8d-1)
      elseif(xext(0).le.1d-3)then
         call SetNumberOfGrids(3)
         call SetGridParameters(1,50,3,xext(0))
         call SetGridParameters(2,40,5,2d-1)
         call SetGridParameters(3,30,5,9d-1)
      else
         call SetNumberOfGrids(2)
         call SetGridParameters(1,50,5,xext(0))
         call SetGridParameters(2,30,5,85d-1)
      endif
      call SetQLimits(0.5d0,Q+eps)
*
*     Initializes integrals on the grids
*
      call InitializeAPFEL
*
*     Compute APFEL evolution operators
*
      call EvolveAPFEL(Q0,Q)
*
*     Interpolation functions
*
      do alphap=0,n
         do igrid=1,ngrid
            k = inter_degree(igrid)
            do alpha=0,nin(igrid)
               inta(igrid,alpha,alphap) = w_int(k,alpha,xext(alphap))
               intb(igrid,alphap,alpha) = 
     1              w_int_ext(n,xext,k,alphap,xg(igrid,alpha))
            enddo
         enddo
      enddo
*
*     Compute Evolution Operator
*
      do alphap=0,n
         do betap=0,n
            do i=0,13
               do j=0,13
                  k = i + 14 * ( j + 14 * ( alphap + (n+1) * betap ) )
                  M(k) = 0d0
               enddo
            enddo
         enddo
      enddo
*
      do alphap=0,n
         do igrid=1,ngrid
            if(xext(alphap).ge.xmin(igrid).and.
     1         xext(alphap).lt.xmin(igrid+1))then
               goto 101
            endif
         enddo
 101     do betap=alphap,n
            do alpha=0,nin(igrid)
               do beta=alpha,nin(igrid)
                  do i=-nff,nff
                     do j=-nfi,nfi
                        k = ( 7 + i ) + 14 * ( ( 7 + j ) 
     1                       + 14 * ( alphap + ( n + 1 ) * betap ) )
                        M(k) = M(k)
     1                       + inta(igrid,alpha,alphap) 
     2                       * PhQCD(igrid,i,j,alpha,beta) 
     3                       * intb(igrid,betap,beta)
                        if(M(k).lt.eps) M(k) = 0d0
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
*
      do alphap=0,n
         call toyLHPDFs(xext(alphap),xpd)
         do i=-6,6
            f0(i,alphap) = xpd(i)
         enddo
      enddo
*
      do alphap=0,n
         do i=-6,6
            fph(0,i,alphap) = 0d0
            do betap=0,n
               do j=-6,6
                  k = ( 7 + i ) + 14 * ( ( 7 + j ) 
     1              + 14 * ( alphap + (n+1) * betap ) )
                  fph(0,i,alphap) = fph(0,i,alphap) + M(k) * f0(j,betap)
               enddo
            enddo
         enddo
      enddo
*
      return
      end
