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
c      character name*(*)
      character*100 name
      logical islhapdf6
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
         if(name(9:9).eq."1")then
            ln = 9
         else
            ln = 8
         endif
      elseif(name(1:12).eq."leptexternal")then
         ln = 12
      elseif(name(1:11).eq."repexternal")then
         ln = 11
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
*     Pretabulated PDFs (for internal use).
*
      elseif(name(1:12).eq."pretabulated")then
         if(name(13:13).eq."1")then
            ln = 13
         else
            ln = 12
         endif
*
*     External LHAPDF grids
*
      else
         if(index(name,"LHgrid").eq.0)then
           ln = len_trim(name)
         else
           ln = index(name,"LHgrid") + 5
           if (islhapdf6().eqv..true.) ln = ln - 7
         endif
      endif
      pdfset = name(1:ln)
      pdfsetlen = ln
      InPDFs = "done"
*
      return
      end
