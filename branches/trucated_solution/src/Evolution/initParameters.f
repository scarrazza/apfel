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
      include "../commons/Welcome.h"
      include "../commons/EvolOp.h"
      include "../commons/scales.h"
      include "../commons/Evs.h"
      include "../commons/Nf_FF.h"
      include "../commons/ipt.h"
      include "../commons/Th.h"
      include "../commons/alpha_ref_QCD.h"
      include "../commons/alpha_ref_QED.h"
      include "../commons/lambda_ref_QCD.h"
      include "../commons/AlphaEvolution.h"
      include "../commons/PDFEvolution.h"
      include "../commons/kren.h"
      include "../commons/mass_scheme.h"
      include "../commons/m2th.h"
      include "../commons/MaxFlavourPDFs.h"
      include "../commons/MaxFlavourAlpha.h"
      include "../commons/grid.h"
      include "../commons/pdfset.h"
      include "../commons/Replica.h"
      include "../commons/lock.h"
      include "../commons/TimeLike.h"
      include "../commons/Smallx.h"
      include "../commons/FastEvol.h"
      include "../commons/EpsTrunc.h"
*
*     Initialize default parameters (those that were not initialized before)
*
      if(InWelcome.ne."done")   call EnableWelcomeMessage(.true.)
      if(InScales.ne."done")    call SetQLimits(0.5d0,1000d0)
      if(InPt.ne."done")        call SetPerturbativeOrder(2)
      if(InEvs.ne."done")       call SetVFNS
      if(InTheory.ne."done")    call SetTheory("QCD")
      if(InFastEvol.ne."done")  call SetFastEvolution(.false.)
      if(InTimeLike.ne."done")  call SetTimeLikeEvolution(.false.)
      if(InSmallx.ne."done")    call SetSmallxResummation(.false.,"NLL")
      if(InAlpQCD.ne."done")    call SetAlphaQCDRef(0.35d0,dsqrt(2d0))
      if(InAlpQED.ne."done")    call SetAlphaQEDRef(7.496252d-3,1.777d0)
      if(InAlphaEvol.ne."done") call SetAlphaEvolution("exact")
      if(InPDFEvol.ne."done")   call SetPDFEvolution("exactmu")
      if(InLambdaQCD.ne."done") call SetLambdaQCDRef(0.220d0,5)
      if(InKren.ne."done")      call SetRenFacRatio(1d0)
      if(InMasses.ne."done")  call SetPoleMasses(dsqrt(2d0),4.5d0,175d0)
      if(InMFP.ne."done")       call SetMaxFlavourPDFs(6)
      if(InMFA.ne."done")       call SetMaxFlavourAlpha(6)
      if(InPDFs.ne."done")      call SetPDFset("ToyLH")
      if(InRep.ne."done")       call SetReplica(0)
      if(InEvolOp.ne."done")    call EnableEvolutionOperator(.false.)
      if(InLock.ne."done")      call LockGrids(.false.)
      if(InEpsTrunc.ne."done")  call SetEpsTrunc(1d0)
      if(InGrid.ne."done")then
         call SetNumberOfGrids(3)
         call SetGridParameters(1,80,3,1d-5)
         call SetGridParameters(2,50,5,1d-1)
         call SetGridParameters(3,40,5,8d-1)
      endif
*
*     Security switchs
*
*     If one of the combined solutions QCD x QED is chose
*     switch of the fast evolution.
*
      if(Th.eq."QCEDP".or.Th.eq."QCEDS".or.
     1   Th.eq."QECDP".or.Th.eq."QECDS".or.
     2   Th.eq."QavDP".or.Th.eq."QavDS")then
         call SetFastEvolution(.false.)
      endif
*
*     If the fast evolution is enabled, disable automatically
*     the computation of the evolution operator
*
      if(FastEvol) call EnableEvolutionOperator(.false.)
*
*     If there is more than one subgrid and one of them is an external grid
*     the external evolution operator cannot be computed
*
      if(ngrid.gt.1.and.ThereAreExtGrids) 
     1     call EnableEvolutionOperator(.false.)
*
*     When the computation of the Evolution Operator is enabled
*     lock the grids by default.
*
      if(EvolOp) call LockGrids(.true.)
*
*     If there are external grids the grids cannot be locked
*
      if(ThereAreExtGrids) call LockGrids(.false.)
*
*     Check the consistency of the input parameters
*
      if(Th.ne."QCD".and.Th.ne."QED".and.
     1   Th.ne."QCEDP".and.Th.ne."QCEDS".and.
     2   Th.ne."QECDP".and.Th.ne."QECDS".and.
     3   Th.ne."QavDP".and.Th.ne."QavDS".and.
     4   Th.ne."QUniD")then
         write(6,*) "Theory unknown:"
         write(6,*) "Theory = ",Th
         write(6,*) "  "
         write(6,*) "The options are:"
         write(6,*) "- 'QCD'"
         write(6,*) "- 'QED'"
         write(6,*) "- 'QCEDP'"
         write(6,*) "- 'QCEDS'"
         write(6,*) "- 'QECDP'"
         write(6,*) "- 'QECDS'"
         write(6,*) "- 'QavDP'"
         write(6,*) "- 'QavDS'"
         write(6,*) "- 'QUniD'"
         write(6,*) "  "
         call exit(-10)
      endif
*
      if(Evs.ne."FF".and.Evs.ne."VF")then
         write(6,*) "Evolution scheme unknown:"
         write(6,*) "Evolution scheme = ",Evs
         write(6,*) "  "
         write(6,*) "The options are:"
         write(6,*) "- 'FF'"
         write(6,*) "- 'VF'"
         write(6,*) "  "
         call exit(-10)
      elseif(Evs.eq."FF")then
         if(Nf_FF.lt.3.or.Nf_FF.gt.6)then
            write(6,*) "Number of active flavours not allowed:"
            write(6,*) "Number of active =",Nf_FF
            write(6,*) "  "
            write(6,*) "The allowed range is [3:6]"
            write(6,*) "  "
            call exit(-10)
         endif
      endif
*
      if(ipt.lt.0.or.ipt.gt.2)then
         write(6,*) "Perturbative order not allowed:"
         write(6,*) "Perturbative order =",ipt
         write(6,*) "  "
         write(6,*) "The allowed range is [0:2]"
         write(6,*) "  "
         call exit(-10)
      endif
*
      if(mass_scheme.ne."Pole".and.mass_scheme.ne."MSbar")then
         write(6,*) "Mass scheme unknown:"
         write(6,*) "Mass scheme = ",mass_scheme
         write(6,*) "  "
         write(6,*) "The options are:"
         write(6,*) "- 'Pole'"
         write(6,*) "- 'MSbar'"
         write(6,*) "  "
         call exit(-10)
      endif
*
      if(AlphaEvol(1:5).ne."exact".and.
     1   AlphaEvol(1:8).ne."expanded".and.
     2   AlphaEvol(1:6).ne."lambda")then
         write(6,*) "Alpha evolution unknown:"
         write(6,*) "Alpha evolution = ",AlphaEvol
         write(6,*) "  "
         write(6,*) "The options are:"
         write(6,*) "- 'exact'"
         write(6,*) "- 'expanded'"
         write(6,*) "- 'lambda'"
         write(6,*) "  "
         call exit(-10)
      endif
*
*     If the alpha solution is "lambda", compute values of LambdaQCD
*     for all the number of flavours.     
*
      if(AlphaEvol(1:6).eq."lambda") call LambdaQCDnf
*
      if(PDFEvol(1:7).ne."exactmu".and.
     1   PDFEvol(1:10).ne."exactalpha".and.
     2   PDFEvol(1:11).ne."expandalpha")then
         write(6,*) "PDF evolution unknown:"
         write(6,*) "PDF evolution = ",PDFEvol
         write(6,*) "  "
         write(6,*) "The options are:"
         write(6,*) "- 'exactmu'"
         write(6,*) "- 'exactalpha'"
         write(6,*) "- 'expandalpha'"
         write(6,*) "  "
         call exit(-10)
      elseif((PDFEvol(1:10).eq."exactalpha".or.
     2        PDFEvol(1:11).eq."expandalpha").and.Th.eq."QUniD")then
         write(6,*) "The unified solution cannot be used with the"
         write(6,*) "'alpha' solution of the DGLAP equation."
         write(6,*) "  "
         call exit(-10)
      endif
*
      if(Smallx)then
         if(TimeLike)then
            write(6,*) "Timelike evolution and Small-x resummation"
            write(6,*) "cannot be combined, switch off one of them."
            write(6,*) "  "
            call exit(-10)
         endif
         if(kren.ne.1d0)then
            write(6,*) "Renormalization scale variation not allowed"
            write(6,*) "if the small-x resummation is enabled."
            write(6,*) "  "
            call exit(-10)
         endif
         if(PDFEvol.eq."expandalpha")then
            write(6,*) "The 'expandalpha' solution of the DGLAP"
            write(6,*) "equation cannot be used if the small-x "
            write(6,*) "resummation is enabled."
            write(6,*) "  "
            call exit(-10)
         endif
         if(LogAcc.ne.0.and.LogAcc.ne.1)then
            write(6,*) "Logarithmic accuracy not allowed:"
            write(6,*) "LogAcc =",LogAcc
            write(6,*) " "
         write(6,*) "The options are:"
         write(6,*) "- 'LL'"
         write(6,*) "- 'NLL'"
         write(6,*) "  "
         call exit(-10)
         endif
      endif
*
*     Print welcome message and report of the parameters (if enabled)
*
      if(Welcome)then
*
         call WelcomeMessage
*
         write(6,*) "Report of the evolution parameters:"
         write(6,*) "  "
*
         write(6,"(a,a,a)") " ",Th," evolution"
*
         if(FastEvol)then
            write(6,*) "Fast evolution enabled"
         endif
*
         if(TimeLike)then
            write(6,*) "Time-like evolution (fragmentation functions)"
         else
            write(6,*) "Space-like evolution (PDFs)"
         endif
         write(6,*) "Solution of the DGLAP equation: ",PDFEvol
         write(6,"(a,a,a,i1,a)") " Evolution scheme = ",Evs,
     1                           "NS at N",ipt,"LO"
*
         write(6,"(a,f7.2,a,f8.2,a)") " Evolution range [",
     1               dsqrt(Q2min)," :",dsqrt(Q2max)," ] GeV"
*
         if(Smallx)then
            if(LogAcc.eq.0) 
     1           write(6,*) "Small-x resummation at LL enabled"
            if(LogAcc.eq.1) 
     1           write(6,*) "Small-x resummation at NLL enabled"
         endif
*
         write(6,*) "Solution of the coupling equations: ",AlphaEvol
         if(AlphaEvol(1:6).eq."lambda")then
            write(6,*) "Lambda reference value:"
            write(6,"(a,i1,a,f10.6,a)")" - LambdaQCD(",n_ref_qcd,
     1                                 ") = ",lambda_ref_qcd," GeV"
         else
            if(Th.eq."QCD")then
               write(6,*) "Coupling reference value:"
               write(6,"(a,f8.4,a,f9.6)")" - AlphaQCD(",
     1              dsqrt(q2_ref_qcd)," GeV) = ",alpha_ref_qcd
            elseif(Th.eq."QED")then
               write(6,*) "Coupling reference value:"
               write(6,"(a,f8.4,a,f9.6)")" - AlphaQED(",
     1              dsqrt(q2_ref_qed)," GeV) = ",alpha_ref_qed
            else
               write(6,*) "Coupling reference values:"
               write(6,"(a,f8.4,a,f9.6)")" - AlphaQCD(",
     1              dsqrt(q2_ref_qcd)," GeV) = ",alpha_ref_qcd
               write(6,"(a,f8.4,a,f9.6)")" - AlphaQED(",
     1              dsqrt(q2_ref_qed)," GeV) = ",alpha_ref_qed
            endif
         endif
*
         write(6,"(a,a,a)") " ",mass_scheme," heavy quark thresholds:"
         write(6,"(a,f6.2,a)") " - mc = ",dsqrt(m2th(4))," GeV"
         write(6,"(a,f6.2,a)") " - mb = ",dsqrt(m2th(5))," GeV"
         write(6,"(a,f6.2,a)") " - mt = ",dsqrt(m2th(6))," GeV"
*
         write(6,*) " "
      endif
*
      return
      end
