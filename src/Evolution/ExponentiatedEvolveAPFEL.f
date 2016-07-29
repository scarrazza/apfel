************************************************************************
*
*     ExponentiatedEvolveAPFEL.f:
*
*     This ruotine computes the evolved PDFs on the grids.
*
************************************************************************
      subroutine ExponentiatedEvolveAPFEL(Q20,Q2)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/Th.h"
      include "../commons/FastEvol.h"
**
*     Input Variables
*
      double precision Q20,Q2
*
      if(FastEvol)then
*     Evolve directly PDFs on the grid
         do igrid=1,ngrid
            if(Th.eq."QCD")then
               call EvolutionQCD(Q20,Q2)
            elseif(Th.eq."QUniD")then
               call EvolutionUnified(Q20,Q2)
            else
               write(6,*) "The fast evolution is currently available"
               write(6,*) "only for the 'QCD', 'QUniD' evolutions."
               write(6,*) "  "
               call exit(-10)
            endif
         enddo
      else
         do igrid=1,ngrid
*     Evaluate evolution operators on the grid
            if(Th.eq."QCD")then
               call EvolutionOperatorsQCD(Q20,Q2)
            elseif(Th.eq."QUniD")then
               call EvolutionOperatorsUnified(Q20,Q2)
            endif
*     Initialize PDFs at the initial scale on the grid
            call initPDFs(Q20)
*     Convolute intial PDFs with the evolution operators
            call EvolvePDFs(igrid)
         enddo
      endif
*     Join all the subgrids.
*     Join also the operators if the production of the evolution operators
*     has been enabled (For the moment available only for QCD).
      call JoinGrids
*
      return
      end
