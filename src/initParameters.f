************************************************************************
*
*     initParameters.f:
*
*     It sets all the evoution parameters if they were not set exernally
*     before by the user.
*
************************************************************************
      subroutine initParameters
*
      implicit none
*
      include "../commons/scales.h"
      include "../commons/Evs.h"
      include "../commons/Nf_FF.h"
      include "../commons/ipt.h"
      include "../commons/Th.h"
      include "../commons/alpha_ref_QCD.h"
      include "../commons/alpha_ref_QED.h"
      include "../commons/kren.h"
      include "../commons/mass_scheme.h"
      include "../commons/m2th.h"
      include "../commons/MaxFlavourPDFs.h"
      include "../commons/MaxFlavourAlpha.h"
      include "../commons/grid.h"
      include "../commons/pdfset.h"
      include "../commons/Replica.h"
*
*     Initialize default parameters (those that were not initialized before)
*
      if(InScales.ne."done") call SetQLimits(0.5d0,400d0)
      if(InPt.ne."done")     call SetPerturbativeOrder(2)
      if(InEvs.ne."done")    call SetVFNS
      if(InTheory.ne."done") call SetTheory("QCD")
      if(InAlpQCD.ne."done") call SetAlphaQCDRef(0.35d0,dsqrt(2d0))
      if(InAlpQED.ne."done") call SetAlphaQEDRef(7.496252d-3,1.777d0)
      if(InKren.ne."done")   call SetRenFacRatio(1d0)
      if(InMasses.ne."done") call SetPoleMasses(dsqrt(2d0),4.5d0,175d0)
      if(InMFP.ne."done")    call SetMaxFlavourPDFs(6)
      if(InMFA.ne."done")    call SetMaxFlavourAlpha(6)
      if(InPDFs.ne."done")   call SetPDFset("ToyLH")
      if(InRep.ne."done")    call SetReplica(0)
      if(InGrid.ne."done")then
         call SetNumberOfGrids(3)
         call SetGridParameters(1,80,3,1d-5)
         call SetGridParameters(2,50,5,1d-1)
         call SetGridParameters(3,40,5,8d-1)
      endif
*
*     Report of the evolution parameters and check consistency
*
      write(6,*) "Report of the evolution parameters:"
      write(6,*) "  "
*
      if(Th.ne."QCD".and.Th.ne."QED".and.
     1   Th.ne."QCEDP".and.Th.ne."QCEDS".and.
     2   Th.ne."QECDP".and.Th.ne."QECDS".and.
     3   Th.ne."QavDP".and.Th.ne."QavDS")then
         write(6,*) "Theory not allowed:"
         write(6,*) "Theory = ",Th
         write(6,*) "Check it in input.dat."
         write(6,*) "  "
         call exit(-10)
      endif
      write(6,"(a,a,a)") " ",Th," evolution"
*
      if(Evs.ne."FF".and.Evs.ne."VF")then
         write(6,*) "Evolution scheme not allowed:"
         write(6,*) "Evolution scheme = ",Evs
         write(6,*) "Check it in input.dat."
         write(6,*) "  "
         call exit(-10)
      elseif(Evs.eq."FF")then
         if(Nf_FF.lt.3.or.Nf_FF.gt.6)then
            write(6,*) "Number of active flavours not allowed:"
            write(6,*) "Number of active =",Nf_FF
            write(6,*) "Check it in input.dat."
            write(6,*) "  "
            call exit(-10)
         endif
      endif
*
      if(ipt.lt.0.or.ipt.gt.2)then
         write(6,*) "Perturbative order not allowed:"
         write(6,*) "Perturbative order =",ipt
         write(6,*) "Check it in input.dat."
         write(6,*) "  "
         call exit(-10)
      endif
      write(6,"(a,a,a,i1,a)") " Evolution scheme = ",Evs,
     1                        "NS at N",ipt,"LO"
*
      write(6,"(a,f10.4,a,f10.4,a)") " Evolution range [ ",dsqrt(Q2min),
     1                              " : ",dsqrt(Q2max)," ] GeV"
*
      write(6,*) "Coupling reference values:"
      write(6,"(a,f8.4,a,f10.6)") " - AlphaQCD(",dsqrt(q2_ref_qcd),
     1                            " GeV) = ",alpha_ref_qcd
      write(6,"(a,f8.4,a,f10.6)") " - AlphaQED(",dsqrt(q2_ref_qed),
     1                            " GeV) = ",alpha_ref_qed
*
      if(mass_scheme.ne."Pole".and.mass_scheme.ne."MSbar")then
         write(6,*) "Mass scheme not allowed:"
         write(6,*) "Mass scheme = ",mass_scheme
         write(6,*) "Check it in input.dat."
         write(6,*) "  "
         call exit(-10)
      endif
      write(6,"(a,a,a)") " ",mass_scheme," heavy quark thresholds:"
      write(6,"(a,f10.4,a)") " - mc = ",dsqrt(m2th(4))," GeV"
      write(6,"(a,f10.4,a)") " - mb = ",dsqrt(m2th(5))," GeV"
      write(6,"(a,f10.4,a)") " - mt = ",dsqrt(m2th(6))," GeV"
*
      write(6,*) " "
*
      return
      end
