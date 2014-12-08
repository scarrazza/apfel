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
*     Force fast evolution
*
      call SetFastEvolution(.true.)
*
*     Compute PDF evolution (if needed)
*
      call EvolveAPFEL(Q0,Q)
*
*     Convolute evolved PDFs with the DIS coefficient functions
*
      call ConvolutePDFsWithCFs(Q)
*
      return
      end
