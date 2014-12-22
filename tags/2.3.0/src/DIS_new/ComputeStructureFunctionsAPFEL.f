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
c      call ConvolutePDFsWithCFs(Q)
*
      return
      end
