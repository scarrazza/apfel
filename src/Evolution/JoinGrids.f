************************************************************************
*
*     JoinGrids.f:
*
*     This routine joins the subgrids used in the calculation in one
*     single grid.
*
************************************************************************
      subroutine JoinGrids
*
      implicit none
*
      include "../commons/EvolOp.h"
      include "../commons/grid.h"
      include "../commons/fph.h"
      include "../commons/EvolutionOperator.h"
      include "../commons/EvolutionMatrices.h"
**
*     Internal Variables
*
      integer i,j,ig
      integer alpha,beta,alphap,betap
      integer jgrid,dgrid,istart,density,offset
      double precision eps
      parameter(eps=1d-12)
*
*     Join PDF grids
*
      ig = -1
      do jgrid=1,ngrid
         do alpha=0,nin(jgrid)
            if(xmin(jgrid+1)-xg(jgrid,alpha).lt.eps) goto 101
            ig = ig + 1
            do i=-6,6
               fph(0,i,ig) = fph(jgrid,i,alpha)
            enddo
            do i=-3,3
               flepton(0,i,ig) = flepton(jgrid,i,alpha)
            enddo
            fgamma(0,ig) = fgamma(jgrid,alpha)
         enddo
 101  enddo
*
*     PDFs for x > 1 (Needed by the interpolation)
*
      do alpha=1,inter_degree(0)
         do i=-6,6
            fph(0,i,nin(0)+alpha) = fph(ngrid,i,nin(ngrid)+alpha)
         enddo
         do i=-3,3
            flepton(0,i,nin(0)+alpha) = 
     1           flepton(ngrid,i,nin(ngrid)+alpha)
         enddo
         fgamma(0,nin(0)+alpha) = fgamma(ngrid,nin(ngrid)+alpha)
      enddo
*
*     If the computation of the external evolution operator has
*     been enabled, join the evolution operators of the sub grid.
*
      if(EvolOp)then
*     
*     Set Evolution Operator to zero
*     
         do alpha=0,nin(0)
            do beta=0,nin(0)
               do i=-6,6
                  do j=-6,6
                     Ph2PhQCD(0,i,j,alpha,beta) = 0
                  enddo
               enddo
            enddo
         enddo
*
*     Fill Evolution Operator
*     
         do alpha=0,nin(0)
*     Determine starting grid
            do dgrid=1,ngrid
               if(alpha.ge.TransitionPoint(dgrid).and.
     1            alpha.lt.TransitionPoint(dgrid+1))then
                  goto 102
               endif
            enddo
            if(alpha.eq.nin(0)) dgrid = ngrid
 102        density = 1
            istart  = alpha
            alphap  = alpha - TransitionPoint(dgrid)
            betap   = alphap
            do jgrid=dgrid,ngrid
               do beta=istart,TransitionPoint(jgrid+1),density
                  do i=0,13
                     do j=0,13
                        Ev2EvQCD(0,i,j,alpha,beta) = 
     1                       Ev2EvQCD(dgrid,i,j,alphap,betap)
                        if(abs(Ev2EvQCD(0,i,j,alpha,beta)).lt.eps)
     1                     Ev2EvQCD(0,i,j,alpha,beta) = 0d0 
                     enddo
                  enddo
*
                  do i=-nff,nff
                     do j=0,13
                        Ev2PhQCD(0,i,j,alpha,beta) = 
     1                       Ev2PhQCD(dgrid,i,j,alphap,betap)
                        if(abs(Ev2PhQCD(0,i,j,alpha,beta)).lt.eps)
     1                     Ev2PhQCD(0,i,j,alpha,beta) = 0d0 
                     enddo
                  enddo
*
                  do i=-nff,nff
                     do j=-nfi,nfi
                        Ph2PhQCD(0,i,j,alpha,beta) = 
     1                       Ph2PhQCD(dgrid,i,j,alphap,betap)
                        if(abs(Ph2PhQCD(0,i,j,alpha,beta)).lt.eps)
     1                     Ph2PhQCD(0,i,j,alpha,beta) = 0d0 
                     enddo
                  enddo
*
                  betap = betap + 1
               enddo
               offset = ( density 
     1              - mod(TransitionPoint(jgrid+1)-istart,density) )
     2              * DensityFactor(jgrid+1)
               istart  = offset + TransitionPoint(jgrid+1)
               density = density * DensityFactor(jgrid+1)
            enddo
         enddo
      endif
*
c      open(unit=19,file="JointGrid.dat",status="unknown")
c      write(19,*) nin(0)!+inter_degree(0)
c      do alpha=0,nin(0)!+inter_degree(0)
c         write(19,*) alpha,xg(0,alpha)
c      enddo
c      close(19)
*
      return
      end
