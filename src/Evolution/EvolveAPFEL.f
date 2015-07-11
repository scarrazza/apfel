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
      include "../commons/InAPFEL.h"
      include "../commons/scales.h"
      include "../commons/PDFEvolution.h"
**
*     Input Variables
*
      double precision Q0,Q
**
*     Internal Variables
*
      double precision Q20,Q2
      double precision t1,t2
      double precision tol
      parameter(tol=1d-2)
*
*     Check that the APFEL evolution has been initialized
*
      if(InAPFEL.ne."done")then
         write(6,*) "EvolveAPFEL: impossible to perform the evolution,",
     1              " APFEL has not been initialized."
         write(6,*) "Call 'InitializeAPFEL' before calling",
     1              " 'EvolveAPFEL'"
         write(6,*) "   "
         call exit(-10)
      endif
*
      Q20   = Q0 * Q0
      Q2    = Q * Q
      Q2ini = Q20
      Q2fin = Q2
*
      if(Q20.lt.Q2min-tol.or.Q20.gt.Q2max+tol)then
         write(6,*) "Initial energy out of range:"
         write(6,*) "- Q0   =",Q0," GeV"
         write(6,*) "- Qmin =",dsqrt(Q2min)," GeV"
         write(6,*) "- Qmax =",dsqrt(Q2max)," GeV"
         write(6,*) "  "
         call exit(-10)
      elseif(Q2.lt.Q2min-tol.or.Q2.gt.Q2max+tol)then
         write(6,*) "Final energy out of range:"
         write(6,*) "- Q    =",Q," GeV"
         write(6,*) "- Qmin =",dsqrt(Q2min)," GeV"
         write(6,*) "- Qmax =",dsqrt(Q2max)," GeV"
         write(6,*) "  "
         call exit(-10)
      endif
*
      call cpu_time(t1)
*
      if(PDFevol(1:9).eq."truncated")then
         call TruncatedEvolveAPFEL(Q20,Q2)
      else
         call ExponentiatedEvolveAPFEL(Q20,Q2)
      endif
*
      call cpu_time(t2)
*
c      write(6,"(a,f7.3,a)") " Evolution completed in",t2-t1," s"
c      write(6,*) " "
*
      return
      end
