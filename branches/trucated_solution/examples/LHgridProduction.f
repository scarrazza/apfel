************************************************************************
*
*     Example program that produces an LHgrid file
*
************************************************************************
      program LHgridProduction
*
      implicit none
*
      double precision Qin
*
      call SetFastEvolution(.true.)
      call SetPerturbativeOrder(1)
      call SetAlphaQCDRef(0.118d0,91.2d0)
      call SetAlphaEvolution("expanded")
      call SetPDFEvolution("expandalpha")
      call SetPoleMasses(dsqrt(2d0),4.75d0,175d0)
      call SetPDFset("NNPDF23_nlo_as_0118.LHgrid")
*
      Qin = dsqrt(2d0)
      call LHAPDFgrid(100,Qin,"NNPDF23_nlo_as_0118_expanded")
*
      end
