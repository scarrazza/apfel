************************************************************************
*
*     Example program that produces an LHgrid file with the derivivative
*     of the input set.
*
************************************************************************
      program LHgridDerivativeProduction
*
      implicit none
*
      call SetPerturbativeOrder(0)
      call SetPDFSet("NNPDF23_nlo_as_0118")
*
      call LHAPDFgridDerivative(0,"NNPDF23_nlo_as_0118_derived")
*
      end
