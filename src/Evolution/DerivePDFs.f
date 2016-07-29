************************************************************************
*
*     DerivePDFs.f:
*
*     This routine combines the derivative operators computed by:
*
*     - DerivativeOperatorsQCD
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
      double precision dfevQCD(0:13,0:nint_max),dfphQCD(-6:6,0:nint_max)
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
      endif
*
      return
      end
