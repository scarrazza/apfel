************************************************************************
*
*     ComputeStructureFunctionsAPFEL.f:
*
*     This routine computes the structure on the grid by convoluting 
*     the PDFs on the grid with the precomputed coefficient functions.
*
************************************************************************
      subroutine ComputeStructureFunctionsAPFEL(Q0,Q)
*
      implicit none
*
      include "../commons/InAPFELDIS.h"
      include "../commons/kfacQ.h"
      include "../commons/MassScheme.h"
      include "../commons/ipt.h"
      include "../commons/DynScVar.h"
**
*     Input Variables
*
      double precision Q0,Q
      double precision t1,t2
**
*     Internal Variables
*
      integer iptbkp
      double precision muF0,muF
*
*     Check that the APFEL DIS module has been initialized
*
      if(InAPFELDIS.ne."done")then
         write(6,*) "ComputeStructureFunctionsAPFEL: impossible to",
     1              " compute structure functions,",
     2              " APFEL DIS has not been initialized."
         write(6,*) "Call 'InitializeAPFEL_DIS' before calling",
     1              " 'ComputeStructureFunctionsAPFEL'"
         write(6,*) "   "
         call exit(-10)
      endif
*
      call cpu_time(t1)
*
*     If the dynamical scale variation is enabled set muR = muF
*
      if(DynScVar) call SetRenQRatio(dsqrt(kfacQ))
*
*     Compute PDF evolution in the factorization scales
*
      muF0 = Q0
      muF  = dsqrt(kfacQ) * Q
*
*     Scale down the perturbative order of the PDF evolution if one
*     of the FFNSs has been chosen
*
      if(MassScheme(1:3).eq."FFN")then
         iptbkp = ipt
         call SetPerturbativeOrder(max(0,ipt-1))
         call EvolveAPFEL(muF0,muF)
         call SetPerturbativeOrder(iptbkp)
      else
         call EvolveAPFEL(muF0,muF)
      endif
*
*     Convolute evolved PDFs with the DIS coefficient functions
*
      call ComputeDISOperators(Q)
      call ConvolutePDFsWithDISOperators
*
      call cpu_time(t2)
*
c      write(6,"(a,a,f9.5,a)") " Computation of the DIS structure ",
c     1                        "functions completed in",t2-t1," s"
c      write(6,*) " "
*
      return
      end
