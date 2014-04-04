************************************************************************
*
*     EvolveAPFEL.f:
*
*     This ruotine computes the evolved PDFs on the grids.
*
************************************************************************
      subroutine EvolveAPFEL(Q0,Q)
*
      implicit none
*
      include "../commons/scales.h"
      include "../commons/grid.h"
      include "../commons/Th.h"
**
*     Input Variables
*
      double precision Q0,Q
**
*     Internal Variables
*
      double precision Q20,Q2
      double precision t1,t2
*
      Q20   = Q0 * Q0
      Q2    = Q * Q
      Q2ini = Q20
      Q2fin = Q2
*
      if(Q20.lt.Q2min.or.Q20.gt.Q2max)then
         write(6,*) "Initial energy out of range:"
         write(6,*) "- Q0   =",Q0
         write(6,*) "- Qmin =",dsqrt(Q2min)
         write(6,*) "- Qmax =",dsqrt(Q2max)
         call exit(-10)
      elseif(Q2.lt.Q2min.or.Q2.gt.Q2max)then
         write(6,*) "Final energy out of range:"
         write(6,*) "- Q    =",Q
         write(6,*) "- Qmin =",dsqrt(Q2min)
         write(6,*) "- Qmax =",dsqrt(Q2max)
         call exit(-10)
      endif
*
      call cpu_time(t1)
      do igrid=1,ngrid
*     Evaluate evolution operators on the grid
         if(Th.eq."QCD")then
            call EvolutionOperatorsQCD(Q20,Q2)
         elseif(Th.eq."QED")then
            call EvolutionOperatorsQED(Q20,Q2)
         elseif(Th.eq."QCEDP".or.Th.eq."QCEDS".or.
     1          Th.eq."QECDP".or.Th.eq."QECDS".or.
     2          Th.eq."QavDP".or.Th.eq."QavDS")then
            call EvolutionOperatorsQCD(Q20,Q2)
            call EvolutionOperatorsQED(Q20,Q2)
         elseif(Th.eq."QUniD")then
            call EvolutionOperatorsUnified(Q20,Q2)
         endif
*     Initialize PDFs at the initial scale on the grid
         call initPDFs(Q20)
*     Convolute intial PDFs with the evolution operators
         call EvolvePDFs(igrid)
      enddo
*     Join all the subgrids
*     (For the moment available only for QCD)
      call JoinGrids
*
      call cpu_time(t2)
*
c      write(6,*) "Evolution done in",t2-t1," s"
c      write(6,*) " "
*
      return
      end
