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
*
*     Hardcoded fragmentation functions
*
*     HKNS
*
      elseif(name(1:15).eq."hknsff07_pip_lo")then
         ln = 15
      elseif(name(1:16).eq."hknsff07_pip_nlo")then
         ln = 16
      elseif(name(1:15).eq."hknsff07_pim_lo")then
         ln = 15
      elseif(name(1:16).eq."hknsff07_pim_nlo")then
         ln = 16
      elseif(name(1:14).eq."hknsff07_Kp_lo")then
         ln = 14
      elseif(name(1:15).eq."hknsff07_Kp_nlo")then
         ln = 15
      elseif(name(1:14).eq."hknsff07_Km_lo")then
         ln = 14
      elseif(name(1:15).eq."hknsff07_Km_nlo")then
         ln = 15
      elseif(name(1:13).eq."hknsff07_p_lo")then
         ln = 13
      elseif(name(1:14).eq."hknsff07_p_nlo")then
         ln = 14
      elseif(name(1:14).eq."hknsff07_pb_lo")then
         ln = 14
      elseif(name(1:15).eq."hknsff07_pb_nlo")then
         ln = 15
*
*     DSS
*
      elseif(name(1:10).eq."dss_pip_lo")then
         ln = 10
      elseif(name(1:11).eq."dss_pip_nlo")then
         ln = 11
      elseif(name(1:10).eq."dss_pim_lo")then
         ln = 10
      elseif(name(1:11).eq."dss_pim_nlo")then
         ln = 11
      elseif(name(1:9).eq."dss_Kp_lo")then
         ln = 9
      elseif(name(1:10).eq."dss_Kp_nlo")then
         ln = 10
      elseif(name(1:9).eq."dss_Km_lo")then
         ln = 9
      elseif(name(1:10).eq."dss_Km_nlo")then
         ln = 10
      elseif(name(1:8).eq."dss_p_lo")then
         ln = 8
      elseif(name(1:9).eq."dss_p_nlo")then
         ln = 9
      elseif(name(1:9).eq."dss_pb_lo")then
         ln = 9
      elseif(name(1:10).eq."dss_pb_nlo")then
         ln = 10
      elseif(name(1:8).eq."dss_h_lo")then
         ln = 8
      elseif(name(1:9).eq."dss_h_nlo")then
         ln = 9
      elseif(name(1:9).eq."dss_hb_lo")then
         ln = 9
      elseif(name(1:10).eq."dss_hb_nlo")then
         ln = 10
*
*     Kretzer's parametrization at Q2 = 0.4 GeV^2 of the light partons
*     for pi+ taken at NLO from hep-ph/0003177 (used for the benchmark).
*
      elseif(name(1:7).eq."kretzer")then
         ln = 7
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
