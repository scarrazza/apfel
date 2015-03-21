************************************************************************
*
*     SetPDFSet.f:
*
*     This subroutine sets the name of the PDF set to be used at the
*     initial scale.
*
************************************************************************
      subroutine SetPDFSet(name)
*
      implicit none
*
      include "../commons/pdfset.h"
*
*     Variables
*
      integer ln
      character*100 name
      logical islhapdf6
*
c      ln = index(name,char(0)) - 1
*
*     Internal PDFs
*
      if(name(1:5).eq."ToyLH")then
         ln = 5
      elseif(name(1:7).eq."private")then
         ln = 7
      elseif(name(1:5).eq."apfel")then
         ln = 5
      elseif(name(1:8).eq."external")then
         ln = 8
      elseif(name(1:12).eq."leptexternal")then
         ln = 12
*
*     Kretzer's parametrization at Q2 = 0.4 GeV^2 of the light partons
*     for pi+ taken at NLO from hep-ph/0003177.
*
      elseif(name(1:7).eq."kretzer")then
         ln = 7
*
*     HKNS parametrization at Q2 = 1 GeV^2 of the light partons
*     for pi+ taken at NLO from hep-ph/0702250 (used for the benchmark against MELA).
*
      elseif(name(1:4).eq."MELA")then
         ln = 4
*
*     External LHAPDF grids
*
      else
         ln = index(name,"LHgrid") + 5
         if(ln.eq.5)then
            write(6,*) "Unknown input set:"
            write(6,*) "set = ",name(1:10)
            write(6,*) "  "
            call exit(-10)
         endif
         if (islhapdf6().eqv..true.) ln = ln - 7         
      endif
      pdfset = name(1:ln)
      InPDFs = "done"
*
      return
      end
