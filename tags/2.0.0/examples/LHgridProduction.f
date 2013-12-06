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
c      call SetFFNS(3)
c      call SetPerturbativeOrder(1)
c      call SetPoleMasses(dsqrt(2d0),1d5,1d5)
*
      Qin = dsqrt(2d0)
      call LHAPDFgrid(0,Qin,"ApfelPDFs")
*
      end
