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
      include "../commons/integrals.h"
      include "../commons/ipt.h"
**
*     Internal Variables
*
      integer i,j,ig
      integer alpha,beta,alphap,betap
      integer jgrid,dgrid,istart,density,offset
      integer inf,isp,pt
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
*     Evolution operators
                  do i=0,13
                     do j=0,13
                        Ev2EvQCD(0,i,j,alpha,beta) =
     1                       Ev2EvQCD(dgrid,i,j,alphap,betap)
                        if(abs(Ev2EvQCD(0,i,j,alpha,beta)).lt.eps)
     1                     Ev2EvQCD(0,i,j,alpha,beta) = 0d0 
                     enddo
                  enddo
*
                  do i=-7,6
                     do j=0,13
                        Ev2PhQCD(0,i,j,alpha,beta) =
     1                       Ev2PhQCD(dgrid,i,j,alphap,betap)
                        if(abs(Ev2PhQCD(0,i,j,alpha,beta)).lt.eps)
     1                     Ev2PhQCD(0,i,j,alpha,beta) = 0d0 
                     enddo
                  enddo
*
                  do i=-7,6
                     do j=-7,6
                        Ph2PhQCD(0,i,j,alpha,beta) =
     1                       Ph2PhQCD(dgrid,i,j,alphap,betap)
                        if(abs(Ph2PhQCD(0,i,j,alpha,beta)).lt.eps)
     1                     Ph2PhQCD(0,i,j,alpha,beta) = 0d0 
                     enddo
                  enddo
*     Splitting functions
                  if(IsExt(dgrid))then
                     do inf=nfi,nff,sgn
                        do isp=1,7
                           do pt=0,ipt
                              SP(0,inf,isp,pt,alpha,beta) =
     1                             SP(dgrid,inf,isp,pt,alphap,betap)
                           enddo
                        enddo
                     enddo
                  else
                     do inf=nfi,nff,sgn
                        do isp=1,7
                           do pt=0,ipt
                              SP(0,inf,isp,pt,alpha,beta) =
     1                             SP(dgrid,inf,isp,pt,0,betap-alphap)
                           enddo
                        enddo
                     enddo
                  endif
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
      return
      end
