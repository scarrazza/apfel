************************************************************************
*
*     JoinGrids.f:
*
*     This routine joins the subgrids used in the calculation in one
*     single grid.
*     (For the moment it works only for QCD subgrids)
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
      integer i,j
      integer alpha,beta,alphap,betap
      integer jgrid,dgrid,istart,density,offset
      double precision eps
      parameter(eps=1d-14)
*
*     Join grids
*
      nin(0) = -1
      TransitionPoint(1) = 0
      do jgrid=1,ngrid
         do alpha=0,nin(jgrid)
            if(xmin(jgrid+1)-xg(jgrid,alpha).lt.eps)then
               TransitionPoint(jgrid+1) = nin(0) + 1
               goto 101
            endif
            nin(0) = nin(0) + 1
*     Joining x-space grid ...
            xg(0,nin(0)) = xg(jgrid,alpha)
*     Joining PDFs ...
            do i=-6,6
               fph(0,i,nin(0)) = fph(jgrid,i,alpha)
            enddo
            fgamma(0,nin(0)) = fgamma(jgrid,alpha)
         enddo
 101  enddo
      TransitionPoint(ngrid+1) = nin(0)
*
*     (Indicative) interpolation degree of the joint grid
*     just to be used in the check below.
*
      inter_degree(0) = inter_degree(ngrid)
*
*     Grid and PDFs for x > 1 (Needed by the interpolation)
*
      do alpha=1,inter_degree(0)
         xg(0,nin(0)+alpha) = xg(ngrid,nin(ngrid)+alpha)
         do i=-6,6
            fph(0,i,nin(0)+alpha) = fph(ngrid,i,nin(ngrid)+alpha)
         enddo
         fgamma(0,nin(0)+alpha) = fgamma(ngrid,nin(ngrid)+alpha)
      enddo
*
      if(nin(0)+inter_degree(0).gt.nint_max)then
         write(6,*) "In JoinGrids.f:"
         write(6,*) "Number of points of the joint grid too large:"
         write(6,*) "Maximum value allowed =",nint_max
         write(6,*) "You should reduce it"
         write(6,*) " "
         call exit(-10)
      endif
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
