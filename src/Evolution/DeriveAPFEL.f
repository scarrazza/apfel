************************************************************************
*
*     DeriveAPFEL.f:
*
*     This ruotine computes the derivative of PDFs on the grids.
*
************************************************************************
      subroutine DeriveAPFEL(Q)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/Th.h"
**
*     Input Variables
*
      double precision Q
**
*     Internal Variables
*
      double precision Q2
*
      Q2 = Q * Q
*
      do igrid=1,ngrid
*     Initialize PDFs at the initial scale on the grid
         call initPDFs(Q2)
*     Evaluate evolution operators on the grid
         if(Th.eq."QCD")then
            call DerivativeOperatorsQCD(Q2)
         elseif(Th.eq."QED")then
            call DerivativeOperatorsQED(Q2)
         elseif(Th.eq."QCEDP".or.Th.eq."QCEDS".or.
     1          Th.eq."QECDP".or.Th.eq."QECDS".or.
     2          Th.eq."QavDP".or.Th.eq."QavDS")then
            call DerivativeOperatorsQCD(Q2)
            call DerivativeOperatorsQED(Q2)
         elseif(Th.eq."QUniD")then
            write(6,*) "Derivative for the unified evolution"
            write(6,*) "not implemented yet."
         endif
*     Convolute intial PDFs with the evolution operators
         call DerivePDFs(igrid)
      enddo
*
      return
      end
