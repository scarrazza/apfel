************************************************************************
*
*     ReportParameters.f:
*
*     This subrotine reports the evoution parameters.
*
************************************************************************
      subroutine ReportParameters
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
      include "../commons/lock.h"
      include "../commons/TimeLike.h"
      include "../commons/Smallx.h"
      include "../commons/FastEvol.h"
      include "../commons/MassRunning.h"
      include "../commons/TauMass.h"
      include "../commons/LeptEvol.h"
      include "../commons/EpsTrunc.h"
*
*     Report of the parameters (if enabled)
*
      if(Welcome)then
*
         call WelcomeMessage
*
         write(6,*) achar(27)//"[34m"//
     1              "Report of the evolution parameters:"
         write(6,*) "  "
*
         write(6,"(a,a,a)") " ",trim(Th)," evolution"
*
         if(LeptEvol.and.Th.eq."QUniD")then
            write(6,*) "Lepton evolution enabled"
         endif
*
         if(TimeLike)then
            write(6,*) "Time-like evolution (fragmentation functions)"
         else
            write(6,*) "Space-like evolution (PDFs)"
         endif
*
         if(Evs.eq."VF")then
            write(6,"(a,a,i1,a)") " Evolution scheme: ",
     1                            "VFNS at N",ipt,"LO"
         elseif(Evs.eq."FF")then
            write(6,"(a,a,i1,a,i1,a)") " Evolution scheme: ",
     1                                 "FFNS with ",Nf_FF,
     2                                 " active flavours at N",ipt,"LO"
         endif
*
         if(Evs.eq."VF")then
            write(6,"(a,a,a,i1,a)")
     1           " Solution of the DGLAP equation: '",trim(PDFEvol),
     2           "' with maximum ",nfMaxPDFs," active flavours"
         elseif(Evs.eq."FF")then
            write(6,"(a,a,a)")
     1           " Solution of the DGLAP equation: '",trim(PDFEvol),"'"
         endif
         if(PDFevol(1:9).eq."truncated")then
            write(6,"(a,a,es10.3)") " - value of the",
     1           " truncation parameter epsilon =",EpsTrunc
         endif
*
         if(Evs.eq."VF")then
            write(6,"(a,a,a,a,i1,a)") " Solution of the coupling",
     1           " equations: '",trim(AlphaEvol),"' with maximum ",
     2           nfMaxAlpha," active flavours"
         elseif(Evs.eq."FF")then
            write(6,"(a,a,a,a)") " Solution of the coupling",
     1           " equations: '",trim(AlphaEvol),"'"
         endif
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
         if(mass_scheme.eq."MSbar")then
            write(6,*) "MSbar heavy quark reference values:"
            write(6,"(a,f7.3,a,f7.3,a)") " - mc(",dsqrt(Q2th(4)),
     1           " GeV) = ",dsqrt(m2q(4))," GeV"
            write(6,"(a,f7.3,a,f7.3,a)") " - mb(",dsqrt(Q2th(5)),
     1           " GeV) = ",dsqrt(m2q(5))," GeV"
            write(6,"(a,f7.3,a,f7.3,a)") " - mt(",dsqrt(Q2th(6)),
     1           " GeV) = ",dsqrt(m2q(6))," GeV"
            write(6,*) "MSbar heavy quark thresholds:"
            write(6,"(a,f7.3,a)") " - mc(mc) = ",dsqrt(m2th(4))," GeV"
            write(6,"(a,f7.3,a)") " - mb(mb) = ",dsqrt(m2th(5))," GeV"
            write(6,"(a,f7.3,a)") " - mt(mt) = ",dsqrt(m2th(6))," GeV"
            if(MassRunning)then
               write(6,*) "Running of the masses enabled"
            else
               write(6,*) "Running of the masses disabled"
            endif
         elseif(mass_scheme(1:4).eq."Pole")then
            write(6,*) "Pole heavy quark thresholds:"
            write(6,"(a,f7.3,a)") " - Mc = ",dsqrt(m2th(4))," GeV"
            write(6,"(a,f7.3,a)") " - Mb = ",dsqrt(m2th(5))," GeV"
            write(6,"(a,f7.3,a)") " - Mt = ",dsqrt(m2th(6))," GeV"
         endif
*     
         write(6,"(a,f7.4)") " muR / muF = ",dsqrt(kren)
*
         if(LeptEvol.and.Th.eq."QUniD")then
            write(6,"(a,f6.3,a)") " Mass of the tau lepton =",MTau,
     1                            " GeV"
         endif
*
         if(Smallx)then
            if(LogAcc.eq.0) 
     1           write(6,*) "Small-x resummation at LL enabled"
            if(LogAcc.eq.1) 
     1           write(6,*) "Small-x resummation at NLL enabled"
         endif
*
         write(6,*) " "
*
         write(6,"(a,f7.2,a,f10.2,a)") " Allowed evolution range [",
     1               dsqrt(Q2min)," :",dsqrt(Q2max)," ] GeV"
*
         if(Lock)then
            write(6,*) "The internal subgrids will be locked"
         endif
*
         if(FastEvol)then
            write(6,*) "Fast evolution enabled"
         endif
*
         if(EvolOp)then
            write(6,*) "Computation of the evolution operator enabled"
         endif
*
         write(6,*) achar(27)//"[0m"
      endif
*
      return
      end
