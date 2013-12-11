************************************************************************
*
*     EvolveGridAPFEL.f:
*
*     This ruotine computes the PDFs independent interpolation grids.
*
************************************************************************
      subroutine EvolveGridAPFEL(name,Q0,Q)
*
      implicit none
*
      include "../commons/scales.h"
      include "../commons/grid.h"
      include "../commons/Th.h"
      include "../commons/f0ph.h"
      include "../commons/fph.h"
**
*     Input Variables
*
      double precision Q0,Q
      character*100 name
**
*     Internal Variables
*
      integer alpha,beta,gamma,i,n
      double precision x
      double precision Q20,Q2
      double precision t1,t2
      double precision w_int
      double precision xGridPDF(ngrid_max,-6:6,0:nint_max,0:nint_max)
      double precision xGridGamma(ngrid_max,0:nint_max,0:nint_max)
*
      Q20 = Q0 * Q0
      Q2  = Q * Q
      Q2ini = Q20
      Q2fin = Q2
*
      if(Q20.lt.Q2min.or.Q20.gt.Q2max)then
         write(6,*) "Initial Energy out of range:"
         write(6,*) "- Q0   =",Q0
         write(6,*) "- Qmin =",dsqrt(Q2min)
         write(6,*) "- Qmax =",dsqrt(Q2max)
         call exit(-10)
      elseif(Q2.lt.Q2min.or.Q2.gt.Q2max)then
         write(6,*) "Final Energy out of range:"
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
         endif
*
*     Construction of the PDF-independent interpolation tables
*
         n = inter_degree(igrid)
         do beta=0,nin(igrid)
*     Set as a PDF set the "delta" PDFs
            call SetPDFSet("delta")
*     Call the alpha-th replica
            call SetReplica(beta)
*     Initialize PDFs at the initial scale on the grid
            call initPDFs(Q20)
*     Convolute intial PDFs with the evolution operators
            call EvolvePDFs(igrid)
            do alpha=0,nin(igrid)
               x = xg(igrid,alpha)
               do i=-6,6
                  xGridPDF(igrid,i,alpha,beta) = 0d0
                  do gamma=0,nin(igrid)
                     xGridPDF(igrid,i,alpha,beta) = 
     1                    xGridPDF(igrid,i,alpha,beta) 
     2                    + w_int(n,gamma,x) * fph(igrid,i,gamma)
                  enddo
               enddo
               xGridGamma(igrid,alpha,beta) = 0d0
               do gamma=0,nin(igrid)
                  xGridGamma(igrid,alpha,beta) = 
     1                 xGridGamma(igrid,alpha,beta) 
     2                 + w_int(n,gamma,x) * fgamma(igrid,gamma)
               enddo
            enddo
         enddo
      enddo
*
      do igrid=1,ngrid
*     Set the PDF set
         call SetPDFSet(name)
*     Initialize PDFs at the initial scale on the grid
         call initPDFs(Q20)
*     Multply the evolution grids with initial PDFs
         do alpha=0,nin(igrid)
            do i=-6,6
               fph(igrid,i,alpha) = 0d0
            enddo
            fgamma(igrid,alpha) = 0d0
         enddo
         do alpha=0,nin(igrid)
            do beta=0,nin(igrid)
               do i=-6,6
                  fph(igrid,i,alpha) = fph(igrid,i,alpha)
     1                               + xGridPDF(igrid,i,alpha,beta)
     2                               * f0ph(i,beta)
               enddo
               fgamma(igrid,alpha) = fgamma(igrid,alpha)
     1                               + xGridGamma(igrid,alpha,beta)
     2                               * f0bos(beta)
            enddo
         enddo
      enddo
      call cpu_time(t2)
*
c      write(6,*) "Evolution done in",t2-t1," s"
c      write(6,*) " "
*
      return
      end
