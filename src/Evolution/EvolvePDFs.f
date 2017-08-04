************************************************************************
*
*     EvolvePDFs.f:
*
*     This routine combines the evolution operators computed by:
*
*     - EvolutionOperatorsQCD
*     - EvolutionOperatorsUnified
*
*     with the initial PDFs producing the final PDFs on the interpolation
*     grid.
*
************************************************************************
      subroutine EvolvePDFs(jgrid)
*
      implicit none
*
      include "../commons/Th.h"
      include "../commons/grid.h"
      include "../commons/EvolutionMatrices.h"
      include "../commons/f0ph.h"
      include "../commons/fph.h"
      include "../commons/EvolOp.h"
**
*     Input Variables
*
      integer jgrid
**
*     Internal Variables
*
      integer inf,inl
      integer i
      integer alpha
      double precision fevQCD(0:13,0:nint_max),fphQCD(-6:6,0:nint_max)
      double precision fevUni(0:13,0:nint_max),flevUni(6,0:nint_max)
      double precision fphUni(-6:6,0:nint_max),flphUni(-3:3,0:nint_max)
*
************************************************************************
*     QCD evolution
************************************************************************
      if(Th.eq."QCD")then
*     Rotate initial PDFs from physical to QCD evolution basis
         call PDFphys2evQCD(f0ph,fevQCD)
*     Evolve PDFs using the QCD evolution operators
         do inf=nfi,nff,sgn
            call EvolveQCD(inf,fevQCD)
         enddo
*     Rotate evolved PDFs from QCD evolution to physical basis
         call PDFevQCD2phys(fevQCD,fphQCD)
*     Compute the total evolution operator in case required
         if(EvolOp) call JoinOperatorsQCD(jgrid)
*     Put Evolved PDF into the common "fph"
         do alpha=0,nin(jgrid)
            do i=-6,6
               fph(jgrid,i,alpha) = fphQCD(i,alpha)
            enddo
            do i=-3,3
               flepton(jgrid,i,alpha) = f0lep(i,alpha)
            enddo
            fgamma(jgrid,alpha) = f0lep(0,alpha)
         enddo
*
************************************************************************
*     Unified QCD x QED Evolution
************************************************************************
      elseif(Th.eq."QUniD")then
*     Rotate initial PDFs from physical to Unified evolution basis
         call PDFphys2evUni(f0lep,f0ph,flevUni,fevUni)
*     Evolve PDFs using the QCD evolution operators
         do inl=nli,nlf,sgn
            do inf=nfli(inl),nflf(inl),sgn
               call EvolveUni(inf,inl,flevUni,fevUni)
            enddo
         enddo
*     Rotate evolved PDFs from Unified evolution to physical basis
         call PDFevUni2phys(flevUni,fevUni,flphUni,fphUni)
*     Compute the total evolution operator in case required
         if(EvolOp) call JoinOperatorsUni(jgrid)
*     Put Evolved PDF into the common "fph"
         do alpha=0,nin(jgrid)
            do i=-6,6
               fph(jgrid,i,alpha) = fphUni(i,alpha)
            enddo
            do i=-3,3
               flepton(jgrid,i,alpha) = flphUni(i,alpha)
            enddo
            fgamma(jgrid,alpha) = flphUni(0,alpha)
         enddo
      endif
*
      return
      end
