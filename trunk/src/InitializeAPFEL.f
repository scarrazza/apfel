************************************************************************
*
*     InitializeAPFEL.f:
*
*     This routine initializes the integrals needed to evolve PDFs
*     on the grids.
*
************************************************************************
      subroutine InitializeAPFEL
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/m2th.h"
      include "../commons/Th.h"
      include "../commons/Evs.h"
      include "../commons/Nf_FF.h"
      include "../commons/scales.h"
*
*     Variables
*
      integer inf,nfi,nff
      double precision t1,t2
*
      call cpu_time(t1)
*
*     Print welcome message
*
      call WelcomeMessage
*
*     Read input parameters
*
      call initParameters
*
      do igrid=1,ngrid
*
*     Initialize grid
*
         call initGrid
*
*     Initalize spliting functions integral matrices for the given nf
*
*     Fixed Flavour Number Scheme
         if(Evs.eq."FF")then
            nfi = Nf_FF
            nff = Nf_FF
*     Variable Flavour Number Scheme
         elseif(Evs.eq."VF")then
            if(Q2max.gt.m2th(6))then
               nff = 6
            elseif(Q2max.gt.m2th(5))then
               nff = 5
            elseif(Q2max.gt.m2th(4))then
               nff = 4
            else
               nff = 3
            endif
*
            if(Q2min.gt.m2th(6))then
               nfi = 6
            elseif(Q2min.gt.m2th(5))then
               nfi = 5
            elseif(Q2min.gt.m2th(4))then
               nfi = 4
            else
               nfi = 3
            endif
*     Initialize matching conditions in can the VFNS is chosen
            if(Th.ne."QED") call initIntegralsMatching
         endif
*
*     Evaluate evolution operators on the grid
*
         if(Th.eq."QCD")then
            do inf=nfi,nff
               call initIntegralsQCD(inf)
            enddo
         elseif(Th.eq."QED")then
            do inf=nfi,nff
               call initIntegralsQED(inf)
            enddo
         elseif(Th.eq."QCEDP".or.Th.eq."QCEDS".or.
     1          Th.eq."QECDP".or.Th.eq."QECDS".or.
     2          Th.eq."QavDP".or.Th.eq."QavDS")then
            do inf=nfi,nff
               call initIntegralsQCD(inf)
               call initIntegralsQED(inf)
            enddo
         endif
      enddo
      call cpu_time(t2)
*
      write(6,*) "Initialization done in",t2-t1," s"
      write(6,*) " "
*
      return
      end
