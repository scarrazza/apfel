************************************************************************
*
*     Example program that produces an LHgrid file
*
************************************************************************
      program LHgridProduction
*
      implicit none
*
      call SetPDFSet("NNPDF31_nnlo_as_0118")      
      call LHAPDFgrid(1,1.65d0,"ApfelPDFs")
*
      end
