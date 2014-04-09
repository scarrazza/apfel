************************************************************************
*
*     EvolvePDFs.f:
*
*     This routine combines the evolution operators computed by:
*
*     - EvolutionOperatorsQCD
*     - EvolutionOperatorsQED
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
      integer inf
      integer i
      integer alpha
      double precision fevQCD(0:13,0:nint_max),fevQCDb(0:13,0:nint_max)
      double precision fphQCD(-6:6,0:nint_max),fphQCDb(-6:6,0:nint_max)
      double precision fevQED(0:13,0:nint_max),fevQEDb(0:13,0:nint_max)
      double precision fphQED(-6:6,0:nint_max),fphQEDb(-6:6,0:nint_max)
      double precision fgm(0:nint_max),fgmb(0:nint_max)
      double precision fgl(0:nint_max),fglb(0:nint_max)
      double precision fevUni(0:13,0:nint_max)
      double precision fphUni(-7:6,0:nint_max)
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
         do alpha=0,nin(igrid)
            do i=-6,6
               fph(jgrid,i,alpha) = fphQCD(i,alpha)
            enddo
            fgamma(jgrid,alpha) = f0bos(alpha)
         enddo
*
************************************************************************
*     QED evolution
************************************************************************
      elseif(Th.eq."QED")then
*     Exchage gluon with photon at the initial scale
         call switchGluonPhoton
*     Rotate initial PDFs from physical to QED evolution basis
         call PDFphys2evQED(f0ph,fevQED)
*     Evolve PDFs using the QED evolution operators
         do inf=nfi,nff,sgn
            call EvolveQED(inf,fevQED)
         enddo
*     Rotate evolved PDFs from QED evolution to physical basis
         call PDFevQED2phys(fevQED,fphQED)
*     Put Evolved PDF into the common "fph"
         do alpha=0,nin(igrid)
            do i=1,6
               fph(jgrid,i,alpha)  = fphQED(i,alpha)
               fph(jgrid,-i,alpha) = fphQED(-i,alpha)
            enddo
            fph(jgrid,0,alpha)  = f0bos(alpha)
            fgamma(jgrid,alpha) = fphQED(0,alpha)
         enddo
*
************************************************************************
*     QCD x QED evolution (QCD first and QED second) in parallel
************************************************************************
      elseif(Th.eq."QCEDP")then
*     Rotate initial PDFs from physical to QCD evolution basis
         call PDFphys2evQCD(f0ph,fevQCD)
         do alpha=0,nin(igrid)
            fgm(alpha) = f0bos(alpha)
         enddo
*     Evolve PDFs
         do inf=nfi,nff,sgn
*     QCD evolution first
            call EvolveQCD(inf,fevQCD)
*     Rotate QCD evolved PDFs from QCD evolution to physical basis
            call PDFevQCD2phys(fevQCD,fphQCD)
            do alpha=0,nin(igrid)
*     Put gluon "fphQCD(0,alpha)" into fgl and replace it with fgm
               fgl(alpha)      = fphQCD(0,alpha)
               fphQCD(0,alpha) = fgm(alpha)
            enddo
*     Rotate PDFs from physical to QED evolution basis
            call PDFphys2evQED(fphQCD,fevQED)
*     QED evolution second
            call EvolveQED(inf,fevQED)
*     Rotate QED evolved PDFs from QED evolution to physical basis
            call PDFevQED2phys(fevQED,fphQED)
            do alpha=0,nin(igrid)
*     Put photon "fphQED(0,alpha)" into fgm and replace it with fgl
               fgm(alpha)      = fphQED(0,alpha)
               fphQED(0,alpha) = fgl(alpha)
            enddo
*     Rotate PDFs from physical to QCD evolution basis
            call PDFphys2evQCD(fphQED,fevQCD)
         enddo
*     Rotate initial PDFs from QCD evolution to physical basis
         call PDFevQCD2phys(fevQCD,fphQCD)
         do alpha=0,nin(igrid)
            do i=-6,6
               fph(jgrid,i,alpha) = fphQCD(i,alpha)
            enddo
            fgamma(jgrid,alpha) = fgm(alpha)
         enddo
*
************************************************************************
*     QCD x QED evolution (QCD first and QED second) in series
************************************************************************
      elseif(Th.eq."QCEDS")then
*     Rotate initial PDFs from physical to QCD evolution basis
         call PDFphys2evQCD(f0ph,fevQCD)
*     Evolve PDFs using the QCD evolution operators
         do inf=nfi,nff,sgn
            call EvolveQCD(inf,fevQCD)
         enddo
*     Rotate QCD evolved PDFs from QCD evolution to physical basis
         call PDFevQCD2phys(fevQCD,fphQCD)
         do alpha=0,nin(igrid)
*     Put gluon "fphQCD(0,alpha)" into fgl and replace it with f0bos
            fgl(alpha)      = fphQCD(0,alpha)
            fphQCD(0,alpha) = f0bos(alpha)
         enddo
*     Rotate PDFs from physical to QED evolution basis
         call PDFphys2evQED(fphQCD,fevQED)
*     Evolve PDFs using the QED evolution operators
         do inf=nfi,nff,sgn
            call EvolveQED(inf,fevQED)
         enddo
*     Rotate evolved PDFs from QED evolution to physical basis
         call PDFevQED2phys(fevQED,fphQED)
*     Put Evolved PDF into the common "fph"
         do alpha=0,nin(igrid)
            do i=1,6
               fph(jgrid,i,alpha)  = fphQED(i,alpha)
               fph(jgrid,-i,alpha) = fphQED(-i,alpha)
            enddo
            fph(jgrid,0,alpha)  = fgl(alpha)
            fgamma(jgrid,alpha) = fphQED(0,alpha)
         enddo
*
************************************************************************
*     QED x QCD evolution (QED first and QCD second) in parallel
************************************************************************
      elseif(Th.eq."QECDP")then
*     Exchage gluon with photon at the initial scale
         call switchGluonPhoton
*     Rotate initial PDFs from physical to QED evolution basis
         call PDFphys2evQED(f0ph,fevQED)
         do alpha=0,nin(igrid)
            fgl(alpha) = f0bos(alpha)
         enddo
*     Evolve PDFs
         do inf=nfi,nff,sgn
*     QED evolution first
            call EvolveQED(inf,fevQED)
*     Rotate QED evolved PDFs from QED evolution to physical basis
            call PDFevQED2phys(fevQED,fphQED)
            do alpha=0,nin(igrid)
*     Put photon "fphQED(0,alpha)" into fgm and replace it with fgl
               fgm(alpha)      = fphQED(0,alpha)
               fphQED(0,alpha) = fgl(alpha)
            enddo
*     Rotate PDFs from physical to QCD evolution basis
            call PDFphys2evQCD(fphQED,fevQCD)
*     QCD evolution second
            call EvolveQCD(inf,fevQCD)
*     Rotate QCD evolved PDFs from QCD evolution to physical basis
            call PDFevQCD2phys(fevQCD,fphQCD)
            do alpha=0,nin(igrid)
*     Put gluon "fphQCD(0,alpha)" into fgl and replace it with fgm
               fgl(alpha)      = fphQCD(0,alpha)
               fphQCD(0,alpha) = fgm(alpha)
            enddo
*     Rotate PDFs from physical to QED evolution basis
            call PDFphys2evQED(fphQCD,fevQED)
         enddo
*     Rotate initial PDFs from QED evolution to physical basis
         call PDFevQED2phys(fevQED,fphQED)
         do alpha=0,nin(igrid)
            do i=1,6
               fph(jgrid,i,alpha)  = fphQED(i,alpha)
               fph(jgrid,-i,alpha) = fphQED(-i,alpha)
            enddo
            fph(jgrid,0,alpha)  = fgl(alpha)
            fgamma(jgrid,alpha) = fphQED(0,alpha)
         enddo
*
************************************************************************
*     QED x QCD evolution (QED first and QCD second) in series
************************************************************************
      elseif(Th.eq."QECDS")then
*     Exchage gluon with photon at the initial scale
         call switchGluonPhoton
*     Rotate initial PDFs from physical to QED evolution basis
         call PDFphys2evQED(f0ph,fevQED)
*     Evolve PDFs using the QED evolution operators
         do inf=nfi,nff,sgn
            call EvolveQED(inf,fevQED)
         enddo
*     Rotate evolved PDFs from QED evolution to physical basis
         call PDFevQED2phys(fevQED,fphQED)
         do alpha=0,nin(igrid)
*     Put gluon "fphQED(0,alpha)" into fgm and replace it with f0bos
            fgm(alpha)      = fphQED(0,alpha)
            fphQED(0,alpha) = f0bos(alpha)
         enddo
*     Rotate PDFs from physical to QCD evolution basis
         call PDFphys2evQCD(fphQED,fevQCD)
*     Evolve PDFs using the QCD evolution operators
         do inf=nfi,nff,sgn
            call EvolveQCD(inf,fevQCD)
         enddo
*     Rotate evolved PDFs from QCD evolution to physical basis
         call PDFevQCD2phys(fevQCD,fphQCD)
*     Put Evolved PDF into the common "fph"
         do alpha=0,nin(igrid)
            do i=-6,6
               fph(jgrid,i,alpha) = fphQCD(i,alpha)
            enddo
            fgamma(jgrid,alpha) = fgm(alpha)
         enddo
*
************************************************************************
*     ( QCD x QED + QED x QCD ) / 2 evolution in parallel
************************************************************************
      elseif(Th.eq."QavDP")then
         call PDFphys2evQCD(f0ph,fevQCD)
         do alpha=0,nin(igrid)
            fgm(alpha) = f0bos(alpha)
         enddo
         call switchGluonPhoton
         call PDFphys2evQED(f0ph,fevQED)
         do alpha=0,nin(igrid)
            fgl(alpha) = f0bos(alpha)
         enddo
*
         do inf=nfi,nff,sgn
*        QCD x QED
            call EvolveQCD(inf,fevQCD)
            call PDFevQCD2phys(fevQCD,fphQCD)
            do alpha=0,nin(igrid)
               fglb(alpha)     = fphQCD(0,alpha)
               fphQCD(0,alpha) = fgm(alpha)
            enddo
            call PDFphys2evQED(fphQCD,fevQEDb)
            call EvolveQED(inf,fevQEDb)
            call PDFevQED2phys(fevQEDb,fphQEDb)
            do alpha=0,nin(igrid)
               fgm(alpha)       = fphQEDb(0,alpha)
               fphQEDb(0,alpha) = fglb(alpha)
            enddo
*        QCD x QED
            call EvolveQED(inf,fevQED)
            call PDFevQED2phys(fevQED,fphQED)
            do alpha=0,nin(igrid)
               fgmb(alpha)     = fphQED(0,alpha)
               fphQED(0,alpha) = fgl(alpha)
            enddo
            call PDFphys2evQCD(fphQED,fevQCDb)
            call EvolveQCD(inf,fevQCDb)
            call PDFevQCD2phys(fevQCDb,fphQCDb)
            do alpha=0,nin(igrid)
               fgl(alpha)       = fphQCDb(0,alpha)
               fphQCDb(0,alpha) = fgmb(alpha)
            enddo
*        Take the average
            do alpha=0,nin(igrid)
               do i=1,6
                  fphQCDb(i,alpha)  = ( fphQCDb(i,alpha) 
     1                              + fphQEDb(i,alpha) ) / 2d0
                  fphQCDb(-i,alpha) = ( fphQCDb(-i,alpha) 
     1                              + fphQEDb(-i,alpha) ) / 2d0
                  fphQEDb(i,alpha)  = fphQCDb(i,alpha)
                  fphQEDb(-i,alpha) = fphQCDb(-i,alpha)
               enddo
               fphQCDb(0,alpha) = ( fgm(alpha) + fphQCDb(0,alpha) )/2d0
               fphQEDb(0,alpha) = ( fgl(alpha) + fphQEDb(0,alpha) )/2d0
               fgm(alpha)       = fphQCDb(0,alpha)
               fgl(alpha)       = fphQEDb(0,alpha)
            enddo
            call PDFphys2evQCD(fphQEDb,fevQCD)
            call PDFphys2evQED(fphQCDb,fevQED)
         enddo
         call PDFevQCD2phys(fevQCD,fphQCD)
         call PDFevQED2phys(fevQED,fphQED)
*     Put Evolved PDF into the common "fph"
         do alpha=0,nin(igrid)
            do i=-6,6
               fph(jgrid,i,alpha) = fphQCD(i,alpha)
            enddo
            fgamma(jgrid,alpha) = fgm(alpha)
         enddo
*
************************************************************************
*     ( QCD x QED + QED x QCD ) / 2 evolution in series
************************************************************************
      elseif(Th.eq."QavDS")then
*     QCD x QED
         call PDFphys2evQCD(f0ph,fevQCD)
         do alpha=0,nin(igrid)
            fgm(alpha) = f0bos(alpha)
         enddo
*
         do inf=nfi,nff,sgn
            call EvolveQCD(inf,fevQCD)
            call PDFevQCD2phys(fevQCD,fphQCD)
            do alpha=0,nin(igrid)
               fglb(alpha)     = fphQCD(0,alpha)
               fphQCD(0,alpha) = fgm(alpha)
            enddo
            call PDFphys2evQED(fphQCD,fevQEDb)
            call EvolveQED(inf,fevQEDb)
            call PDFevQED2phys(fevQEDb,fphQEDb)
            do alpha=0,nin(igrid)
               fgm(alpha)       = fphQEDb(0,alpha)
               fphQEDb(0,alpha) = fglb(alpha)
            enddo
            call PDFphys2evQCD(fphQEDb,fevQCD)
         enddo
         call PDFevQCD2phys(fevQCD,fphQCD)
*     QCD x QED
         call switchGluonPhoton
         call PDFphys2evQED(f0ph,fevQED)
         do alpha=0,nin(igrid)
            fgl(alpha) = f0bos(alpha)
         enddo
         do inf=nfi,nff,sgn
            call EvolveQED(inf,fevQED)
            call PDFevQED2phys(fevQED,fphQED)
            do alpha=0,nin(igrid)
               fgmb(alpha)     = fphQED(0,alpha)
               fphQED(0,alpha) = fgl(alpha)
            enddo
            call PDFphys2evQCD(fphQED,fevQCDb)
            call EvolveQCD(inf,fevQCDb)
            call PDFevQCD2phys(fevQCDb,fphQCDb)
            do alpha=0,nin(igrid)
               fgl(alpha)       = fphQCDb(0,alpha)
               fphQCDb(0,alpha) = fgmb(alpha)
            enddo
            call PDFphys2evQED(fphQCDb,fevQED)
         enddo
         call PDFevQED2phys(fevQED,fphQED)
*     Take the average
         do alpha=0,nin(igrid)
            do i=1,6
               fph(jgrid,i,alpha)  = ( fphQCD(i,alpha) 
     1                             + fphQED(i,alpha) ) / 2d0
               fph(jgrid,-i,alpha) = ( fphQCD(-i,alpha) 
     1                             + fphQED(-i,alpha) ) / 2d0
            enddo
            fph(jgrid,0,alpha)  = ( fgl(alpha) + fphQCD(0,alpha) ) / 2d0
            fgamma(jgrid,alpha) = ( fgm(alpha) + fphQED(0,alpha) ) / 2d0
         enddo
*
************************************************************************
*     Unified QCD x QED Evolution
************************************************************************
      elseif(Th.eq."QUniD")then
*     Rotate initial PDFs from physical to Unified evolution basis
         call PDFphys2evUni(f0bos,f0ph,fevUni)
*     Evolve PDFs using the QCD evolution operators
         do inf=nfi,nff,sgn
            call EvolveUni(inf,fevUni)
         enddo
*     Rotate evolved PDFs from Unified evolution to physical basis
         call PDFevUni2phys(fevUni,fphUni)
*     Put Evolved PDF into the common "fph"
         do alpha=0,nin(igrid)
            do i=-6,6
               fph(jgrid,i,alpha) = fphUni(i,alpha)
            enddo
            fgamma(jgrid,alpha) = fphUni(-7,alpha)
         enddo
      endif
*
      return
      end
