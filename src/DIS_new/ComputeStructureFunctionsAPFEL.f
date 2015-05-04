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
**
*     Input Variables
*
      double precision Q0,Q
      double precision t1,t2
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
*     Compute PDF evolution
*
      call EvolveAPFEL(Q0,Q)
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
