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
**
*     Input Variables
*
      double precision Q0,Q
      double precision t1,t2
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
      write(6,"(a,a,f9.5,a)") " Computation of the DIS structure ",
     1                        "functions with PDFs completed in",
     2                        t2-t1," s"
      write(6,*) " "
*
      return
      end
