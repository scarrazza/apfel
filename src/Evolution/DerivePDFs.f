************************************************************************
*
*     DerivePDFs.f:
*
*     This routine combines the derivative operators computed by:
*
*     - DerivativeOperatorsQCD
*     - DerivativeOperatorsQED
*
************************************************************************
      subroutine DerivePDFs(jgrid)
*
      implicit none
*
      include "../commons/Th.h"
      include "../commons/grid.h"
      include "../commons/EvolutionMatrices.h"
      include "../commons/f0ph.h"
      include "../commons/fph.h"
**
*     Input Variables
*
      integer jgrid
**
*     Internal Variables
*
      integer i
      integer alpha
      double precision dfevQCD(0:13,0:nint_max)
      double precision dfevQCDb(0:13,0:nint_max)
      double precision dfphQCD(-6:6,0:nint_max)
      double precision dfphQCDb(-6:6,0:nint_max)
      double precision dfevQED(0:13,0:nint_max)
      double precision dfevQEDb(0:13,0:nint_max)
      double precision dfphQED(-6:6,0:nint_max)
      double precision dfphQEDb(-6:6,0:nint_max)
      double precision dfgm(0:nint_max),dfgmb(0:nint_max)
      double precision dfgl(0:nint_max),dfglb(0:nint_max)
*
************************************************************************
*     QCD evolution
************************************************************************
      if(Th.eq."QCD")then
*     Rotate initial PDFs from physical to QCD evolution basis
         call PDFphys2evQCD(f0ph,dfevQCD)
*     Derive PDFs using the QCD evolution operators
         call DeriveQCD(dfevQCD)
*     Rotate evolved PDFs from QCD evolution to physical basis
         call PDFevQCD2phys(dfevQCD,dfphQCD)
*     Put Evolved PDF into the common "fph"
         do alpha=0,nin(igrid)
            do i=-6,6
               dfph(jgrid,i,alpha) = dfphQCD(i,alpha)
            enddo
            dfgamma(jgrid,alpha) = 0d0
         enddo
*
************************************************************************
*     QED evolution
************************************************************************
      elseif(Th.eq."QED")then
*     Exchage gluon with photon at the initial scale
         call switchGluonPhoton
*     Rotate initial PDFs from physical to QED evolution basis
         call PDFphys2evQED(f0ph,dfevQED)
*     Derive PDFs using the QED evolution operators
         call DeriveQED(dfevQED)
*     Rotate evolved PDFs from QED evolution to physical basis
         call PDFevQED2phys(dfevQED,dfphQED)
*     Put Evolved PDF into the common "fph"
         do alpha=0,nin(igrid)
            do i=1,6
               dfph(jgrid,i,alpha)  = dfphQED(i,alpha)
               dfph(jgrid,-i,alpha) = dfphQED(-i,alpha)
            enddo
            dfph(jgrid,0,alpha)  = 0d0
            dfgamma(jgrid,alpha) = dfphQED(0,alpha)
         enddo
*
************************************************************************
*     QCD x QED evolution (QCD first and QED second)
************************************************************************
      elseif(Th.eq."QCEDP".or.Th.eq."QCEDS")then
*     Rotate initial PDFs from physical to QCD evolution basis
         call PDFphys2evQCD(f0ph,dfevQCD)
         do alpha=0,nin(igrid)
            dfgm(alpha) = f0bos(alpha)
         enddo
*     Derive PDFs
*     QCD first
         call DeriveQCD(dfevQCD)
*     Rotate QCD evolved PDFs from QCD evolution to physical basis
         call PDFevQCD2phys(dfevQCD,dfphQCD)
         do alpha=0,nin(igrid)
*     Put gluon "fphQCD(0,alpha)" into fgl and replace it with fgm
            dfgl(alpha)      = dfphQCD(0,alpha)
            dfphQCD(0,alpha) = dfgm(alpha)
         enddo
*     Rotate PDFs from physical to QED evolution basis
         call PDFphys2evQED(dfphQCD,dfevQED)
*     QED second
         call DeriveQED(dfevQED)
*     Rotate QED evolved PDFs from QED evolution to physical basis
         call PDFevQED2phys(dfevQED,dfphQED)
         do alpha=0,nin(igrid)
*     Put photon "fphQED(0,alpha)" into fgm and replace it with fgl
            dfgm(alpha)      = dfphQED(0,alpha)
            dfphQED(0,alpha) = dfgl(alpha)
         enddo
         do alpha=0,nin(igrid)
            do i=-6,6
               dfph(jgrid,i,alpha) = dfphQED(i,alpha)
            enddo
            dfgamma(jgrid,alpha) = dfgm(alpha)
         enddo
*
************************************************************************
*     QED x QCD evolution (QED first and QCD second)
************************************************************************
      elseif(Th.eq."QECDP".or.Th.eq."QECDS")then
*     Exchage gluon with photon at the initial scale
         call switchGluonPhoton
*     Rotate initial PDFs from physical to QED evolution basis
         call PDFphys2evQED(f0ph,dfevQED)
         do alpha=0,nin(igrid)
            dfgl(alpha) = f0bos(alpha)
         enddo
*     Evolve PDFs
*     QEDfirst
         call DeriveQED(dfevQED)
*     Rotate QED evolved PDFs from QED evolution to physical basis
         call PDFevQED2phys(dfevQED,dfphQED)
         do alpha=0,nin(igrid)
*     Put photon "fphQED(0,alpha)" into fgm and replace it with fgl
            dfgm(alpha)      = dfphQED(0,alpha)
            dfphQED(0,alpha) = dfgl(alpha)
         enddo
*     Rotate PDFs from physical to QCD evolution basis
         call PDFphys2evQCD(dfphQED,dfevQCD)
*     QCD second
         call DeriveQCD(dfevQCD)
*     Rotate QCD evolved PDFs from QCD evolution to physical basis
         call PDFevQCD2phys(dfevQCD,dfphQCD)
         do alpha=0,nin(igrid)
            do i=-6,6
               dfph(jgrid,i,alpha) = dfphQCD(i,alpha)
            enddo
            dfgamma(jgrid,alpha) = dfgm(alpha)
         enddo
*
************************************************************************
*     ( QCD x QED + QED x QCD ) / 2 evolution
************************************************************************
      elseif(Th.eq."QavDP".or.Th.eq."QavDS")then
         call PDFphys2evQCD(f0ph,dfevQCD)
         do alpha=0,nin(igrid)
            dfgm(alpha) = f0bos(alpha)
         enddo
         call switchGluonPhoton
         call PDFphys2evQED(f0ph,dfevQED)
         do alpha=0,nin(igrid)
            dfgl(alpha) = f0bos(alpha)
         enddo
*        QCD x QED
         call DeriveQCD(dfevQCD)
         call PDFevQCD2phys(dfevQCD,dfphQCD)
         do alpha=0,nin(igrid)
            dfglb(alpha)     = dfphQCD(0,alpha)
            dfphQCD(0,alpha) = dfgm(alpha)
         enddo
         call PDFphys2evQED(dfphQCD,dfevQEDb)
         call DeriveQED(dfevQEDb)
         call PDFevQED2phys(dfevQEDb,dfphQEDb)
         do alpha=0,nin(igrid)
            dfgm(alpha)       = dfphQEDb(0,alpha)
            dfphQEDb(0,alpha) = dfglb(alpha)
         enddo
*        QCD x QED
         call DeriveQED(dfevQED)
         call PDFevQED2phys(dfevQED,dfphQED)
         do alpha=0,nin(igrid)
            dfgmb(alpha)     = dfphQED(0,alpha)
            dfphQED(0,alpha) = dfgl(alpha)
         enddo
         call PDFphys2evQCD(dfphQED,dfevQCDb)
         call DeriveQCD(dfevQCDb)
         call PDFevQCD2phys(dfevQCDb,dfphQCDb)
         do alpha=0,nin(igrid)
            dfgl(alpha)       = dfphQCDb(0,alpha)
            dfphQCDb(0,alpha) = dfgmb(alpha)
         enddo
*        Take the average
         do alpha=0,nin(igrid)
            do i=1,6
               dfph(jgrid,i,alpha)  = ( dfphQCDb(i,alpha) 
     1                              + dfphQEDb(i,alpha) ) / 2d0
               dfph(jgrid,-i,alpha) = ( dfphQCDb(-i,alpha) 
     1                              + dfphQEDb(-i,alpha) ) / 2d0
            enddo
            dfph(jgrid,0,alpha) =( dfgm(alpha) + dfphQCDb(0,alpha) )/2d0
            dfgamma(jgrid,alpha)=( dfgl(alpha) + dfphQEDb(0,alpha) )/2d0
         enddo
      endif
*
      return
      end
