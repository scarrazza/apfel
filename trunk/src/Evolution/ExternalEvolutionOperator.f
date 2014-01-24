************************************************************************
*
*     ExternalEvolutionOperator.f:
*
*     This subroutine computes the "external" evolution operator on a
*     user-given grid starting from the "internal" evolution operators
*     on the internal x-space grid of APFEL.
*
************************************************************************
      subroutine ExternalEvolutionOperator(Q0,Q,n,xext,M)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/EvolutionOperator.h"
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
      parameter(eps=1d-14)
c      include "../commons/fph.h"
c      double precision f0(-6:6,0:200),xpd(-6:6)
**
*     Output Variables
*
      double precision M(-7:6,-7:6,0:n,0:n)
*
*     Disable welcome message
*
      call EnableWelcomeMessage(.false.)
*
*     Enable computation of the Total Evolution operator
*
      call EnableEvolutionOperator(.true.)
*
*     Here I should put some ad hoc setting of the APFEL grid
*     based on the kinematical coverage of "xext" to make the 
*     computation more optimal.
*     

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
c      stop
*
*     Compute Evolution Operator
*
      do alphap=0,n
         do betap=0,n
            do i=-6,6
               do j=-6,6
                  M(i,j,alphap,betap) = 0d0
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
c                  do i=-6,6
c                     do j=-6,6
c                        M(i,j,alphap,betap) = M(i,j,alphap,betap)
c     1                       + inta(igrid,alpha,alphap) 
c     2                       * PhQCD(igrid,i,j,alpha,beta) 
c     3                       * intb(igrid,betap,beta)
c                     enddo
c                  enddo
                  M(0,0,alphap,betap) = M(0,0,alphap,betap)
     1                 + inta(igrid,alpha,alphap) 
     2                 * PhQCD(igrid,0,0,alpha,beta) 
     3                 * intb(igrid,betap,beta)
                  do i=1,6
                     M(i,i,alphap,betap) = M(i,i,alphap,betap)
     1                    + inta(igrid,alpha,alphap) 
     2                    * PhQCD(igrid,i,i,alpha,beta) 
     3                    * intb(igrid,betap,beta)
                     M(i,0,alphap,betap) = M(i,0,alphap,betap)
     1                    + inta(igrid,alpha,alphap) 
     2                    * PhQCD(igrid,i,0,alpha,beta) 
     3                    * intb(igrid,betap,beta)
                     M(0,i,alphap,betap) = M(0,i,alphap,betap)
     1                    + inta(igrid,alpha,alphap) 
     2                    * PhQCD(igrid,0,i,alpha,beta) 
     3                    * intb(igrid,betap,beta)

                     M(-i,-i,alphap,betap) = M(-i,-i,alphap,betap)
     1                    + inta(igrid,alpha,alphap) 
     2                    * PhQCD(igrid,-i,-i,alpha,beta) 
     3                    * intb(igrid,betap,beta)
                     M(-i,0,alphap,betap) = M(-i,0,alphap,betap)
     1                    + inta(igrid,alpha,alphap) 
     2                    * PhQCD(igrid,-i,0,alpha,beta) 
     3                    * intb(igrid,betap,beta)
                     M(0,-i,alphap,betap) = M(0,-i,alphap,betap)
     1                    + inta(igrid,alpha,alphap) 
     2                    * PhQCD(igrid,0,-i,alpha,beta) 
     3                    * intb(igrid,betap,beta)
                  enddo
               enddo
            enddo
         enddo
      enddo
c*
c      do alpha=0,n
c         call toyLHPDFs(xext(alpha),xpd)
c         do i=-6,6
c            f0(i,alpha) = xpd(i)
c         enddo
c      enddo
c*
c      do alpha=0,n
c         do i=-6,6
c            fph(0,i,alpha) = 0d0
c            do beta=0,n
c               do j=-6,6
c                  fph(0,i,alpha) = fph(0,i,alpha)
c     1                 + M(i,j,alpha,beta) * f0(j,beta) 
c               enddo
c            enddo
c         enddo
c      enddo
*
      return
      end
