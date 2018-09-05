************************************************************************
*
*     initParameters.f:
*
*     It sets all the evoution parameters if they were not set externally
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
      include "../commons/gridQ.h"
      include "../commons/pdfset.h"
      include "../commons/Replica.h"
      include "../commons/lock.h"
      include "../commons/TimeLike.h"
      include "../commons/Polarized.h"
      include "../commons/Smallx.h"
      include "../commons/FastEvol.h"
      include "../commons/MassRunning.h"
      include "../commons/TauMass.h"
      include "../commons/LeptEvol.h"
      include "../commons/LHAgrid.h"
      include "../commons/EpsTrunc.h"
      include "../commons/NLOQEDCorrections.h"
*
*     Initialize default parameters (those that were not initialized before)
*
      if(InWelcome.ne."done")     call EnableWelcomeMessage(.true.)
      if(InScales.ne."done")      call SetQLimits(1d0,10000d0)
      if(InPt.ne."done")          call SetPerturbativeOrder(2)
      if(InEvs.ne."done")         call SetVFNS
      if(InTheory.ne."done")      call SetTheory("QCD")
      if(InNLOQED.ne."done")      call EnableNLOQEDCorrections(.false.)
      if(InFastEvol.ne."done")    call SetFastEvolution(.true.)
      if(InTimeLike.ne."done")    call SetTimeLikeEvolution(.false.)
      if(InPolarized.ne."done")   call SetPolarizedEvolution(.false.)
      if(InSmallx.ne."done")      call SetSmallxResummation(.false.,
     1                                                      "NLL")
      if(InAlpQCD.ne."done")      call SetAlphaQCDRef(0.35d0,dsqrt(2d0))
      if(InAlpQED.ne."done")      call SetAlphaQEDRef(7.496252d-3,
     1                                                1.777d0)
      if(InLambdaQCD.ne."done")   call SetLambdaQCDRef(0.220d0,5)
      if(InEpsTrunc.ne."done")    call SetEpsilonTruncation(1d-2)
      if(InAlphaEvol.ne."done")   call SetAlphaEvolution("exact")
      if(InPDFEvol.ne."done")     call SetPDFEvolution("exactalpha")
      if(InKren.ne."done")        call SetRenFacRatio(1d0)
      if(InMasses.ne."done")      call SetPoleMasses(dsqrt(2d0),4.5d0,
     1                                               175d0)
      if(InThrRatios.ne."done")   call SetMassMatchingScales(1d0,1d0,
     1                                                       1d0)
      if(InMassRef.ne."done")     call SetMassScaleReference(
     1                                                   dsqrt(m2ph(4)),
     2                                                   dsqrt(m2ph(5)),
     3                                                   dsqrt(m2ph(6)))
      if(InMTau.ne."done")        call SetTauMass(1.777d0)
      if(InMassRunning.ne."done") call EnableMassRunning(.true.)
      if(InMFP.ne."done")         call SetMaxFlavourPDFs(6)
      if(InMFA.ne."done")         call SetMaxFlavourAlpha(6)
      if(InPDFs.ne."done")        call SetPDFset("ToyLH")
      if(InRep.ne."done")         call SetReplica(0)
      if(InEvolOp.ne."done")      call EnableEvolutionOperator(.false.)
      if(InLeptEvol.ne."done")    call EnableLeptonEvolution(.false.)
      if(InLock.ne."done")        call LockGrids(.true.)
      if(InLHgrid.ne."done")      call SetLHgridParameters(100,50,1d-9,
     1                                            1d-1,1d0,50,1d0,1d10)
      if(InQGrid.ne."done")       call SetQGridParameters(70,3)
      if(InGrid.ne."done")then
         call SetNumberOfGrids(3)
         call SetGridParameters(1,80,3,1d-5)
         call SetGridParameters(2,50,5,1d-1)
         call SetGridParameters(3,40,5,8d-1)
      endif
*
*     Check the consistency of the input parameters
*
      if(Th.ne."QCD".and.Th.ne."QUniD")then
         write(6,*) achar(27)//"[31mERROR:"
         write(6,*) "Theory unknown:"
         write(6,*) "Theory = ",Th
         write(6,*) "  "
         write(6,*) "The options are:"
         write(6,*) "- 'QCD'"
         write(6,*) "- 'QUniD'"
         write(6,*) achar(27)//"[0m"
         call exit(-10)
      endif
*
      if(Evs.ne."FF".and.Evs.ne."VF")then
         write(6,*) achar(27)//"[31mERROR:"
         write(6,*) "Evolution scheme unknown:"
         write(6,*) "Evolution scheme = ",Evs
         write(6,*) "  "
         write(6,*) "The options are:"
         write(6,*) "- 'FF'"
         write(6,*) "- 'VF'"
         write(6,*) achar(27)//"[0m"
         call exit(-10)
      endif
*     
      if(Evs.eq."FF".and.(Nf_FF.lt.3.or.Nf_FF.gt.6))then
         write(6,*) achar(27)//"[31mERROR:"
         write(6,*) "Number of active flavours not allowed:"
         write(6,*) "Number of active =",Nf_FF
         write(6,*) "  "
         write(6,*) "The allowed range is [3:6]"
         write(6,*) achar(27)//"[0m"
         call exit(-10)
      endif
*
      if(ipt.lt.0.or.ipt.gt.2)then
         write(6,*) achar(27)//"[31mERROR:"
         write(6,*) "Perturbative order not allowed:"
         write(6,*) "Perturbative order =",ipt
         write(6,*) "  "
         write(6,*) "The allowed range is [0:2]"
         write(6,*) achar(27)//"[0m"
         call exit(-10)
      endif
*
      if(mass_scheme.ne."Pole".and.mass_scheme.ne."MSbar")then
         write(6,*) achar(27)//"[31mERROR:"
         write(6,*) "Mass scheme unknown:"
         write(6,*) "Mass scheme = ",mass_scheme
         write(6,*) "  "
         write(6,*) "The options are:"
         write(6,*) "- 'Pole'"
         write(6,*) "- 'MSbar'"
         write(6,*) achar(27)//"[0m"
         call exit(-10)
      endif
*
      if(AlphaEvol(1:5).ne."exact".and.
     1   AlphaEvol(1:8).ne."expanded".and.
     2   AlphaEvol(1:6).ne."lambda")then
         write(6,*) achar(27)//"[31mERROR:"
         write(6,*) "Alpha evolution unknown:"
         write(6,*) "Alpha evolution = ",AlphaEvol
         write(6,*) "  "
         write(6,*) "The options are:"
         write(6,*) "- 'exact'"
         write(6,*) "- 'expanded'"
         write(6,*) "- 'lambda'"
         write(6,*) achar(27)//"[0m"
         call exit(-10)
      elseif(AlphaEvol(1:5).eq."lambda".and.Th.eq."QUniD")then
         write(6,*) achar(27)//"[31mERROR:"
         write(6,*) "The unified solution cannot be used with the"
         write(6,*) "'lambda' solution of the coupling equations."
         write(6,*) achar(27)//"[0m"
         call exit(-10)
      endif
*
      if(PDFEvol(1:7).ne."exactmu".and.
     1   PDFEvol(1:10).ne."exactalpha".and.
     2   PDFEvol(1:11).ne."expandalpha".and.
     3   PDFEvol(1:9).ne."truncated")then
         write(6,*) achar(27)//"[31mERROR:"
         write(6,*) "PDF evolution unknown:"
         write(6,*) "PDF evolution = ",PDFEvol
         write(6,*) "  "
         write(6,*) "The options are:"
         write(6,*) "- 'exactmu'"
         write(6,*) "- 'exactalpha'"
         write(6,*) "- 'expandalpha'"
         write(6,*) "- 'truncated'"
         write(6,*) achar(27)//"[0m"
         call exit(-10)
      endif
*
      if(Smallx)then
         if(TimeLike)then
            write(6,*) achar(27)//"[31mERROR:"
            write(6,*) "Time-like evolution and Small-x resummation"
            write(6,*) "cannot be combined, switch off one of them."
            write(6,*) achar(27)//"[0m"
            call exit(-10)
         endif
         if(Polarized)then
            write(6,*) achar(27)//"[31mERROR:"
            write(6,*) "Polarized evolution and Small-x resummation"
            write(6,*) "cannot be combined, switch off one of them."
            write(6,*) achar(27)//"[0m"
            call exit(-10)
         endif
         if(kren.ne.1d0)then
            write(6,*) achar(27)//"[31mERROR:"
            write(6,*) "Renormalization scale variation not allowed"
            write(6,*) "if small-x resummation is enabled."
            write(6,*) achar(27)//"[0m"
            call exit(-10)
         endif
         if(LogAcc.ne.0.and.LogAcc.ne.1)then
            write(6,*) achar(27)//"[31mERROR:"
            write(6,*) "Logarithmic accuracy not allowed:"
            write(6,*) "LogAcc =",LogAcc
            write(6,*) " "
            write(6,*) "The options are:"
            write(6,*) "- 'LL'"
            write(6,*) "- 'NLL'"
            write(6,*) achar(27)//"[0m"
            call exit(-10)
         endif
      endif
*
*     Check that there are no identical masses and that
*     masses are correctly ordered
*
      if(m2ph(4).eq.m2ph(5).or.
     1   m2ph(4).eq.m2ph(6).or.
     2   m2ph(5).eq.m2ph(6))then
         write(6,*) achar(27)//"[31mERROR:"
         write(6,*) "There cannot be equal heavy quark masses."
         write(6,*) "This condition is not fulfilled with:"
         write(6,*) "  "
         write(6,"(a,f8.3,a)") " Mc = ",dsqrt(m2ph(4))," GeV"
         write(6,"(a,f8.3,a)") " Mb = ",dsqrt(m2ph(5))," GeV"
         write(6,"(a,f8.3,a)") " Mt = ",dsqrt(m2ph(6))," GeV"
         write(6,*) achar(27)//"[0m"
         call exit(-10)
      endif
*
      if(m2ph(4).gt.m2ph(5).or.
     1   m2ph(4).gt.m2ph(6).or.
     2   m2ph(5).gt.m2ph(6))then
         write(6,*) achar(27)//"[31mERROR:"
         write(6,*) "The heavy quark masses must be ordered, i.e.:"
         write(6,*) "- Mc < Mb < Mt"
         write(6,*) "This condition is not fulfilled with:"
         write(6,*) "  "
         write(6,"(a,f8.3,a)") " Mc = ",dsqrt(m2ph(4))," GeV"
         write(6,"(a,f8.3,a)") " Mb = ",dsqrt(m2ph(5))," GeV"
         write(6,"(a,f8.3,a)") " Mt = ",dsqrt(m2ph(6))," GeV"
         write(6,*) achar(27)//"[0m"
         call exit(-10)
      endif
*
*     Check that also the thresholds are correctly ordered
*
      if(k2th(4)*m2ph(4).gt.k2th(5)*m2ph(5).or.
     1   k2th(4)*m2ph(4).gt.k2th(6)*m2ph(6).or.
     2   k2th(5)*m2ph(5).gt.k2th(6)*m2ph(6))then
         write(6,*) achar(27)//"[31mERROR:"
         write(6,*) "The heavy quark thresholds must be ordered, i.e.:"
         write(6,*) "- Mthc < Mthb < Mtht"
         write(6,*) "This condition is not fulfilled with:"
         write(6,*) "  "
         write(6,"(a,f8.3,a)") " Mthc = ",dsqrt(k2th(4)*m2ph(4))," GeV"
         write(6,"(a,f8.3,a)") " Mthb = ",dsqrt(k2th(5)*m2ph(5))," GeV"
         write(6,"(a,f8.3,a)") " Mtht = ",dsqrt(k2th(6)*m2ph(6))," GeV"
         write(6,*) achar(27)//"[0m"
         call exit(-10)
      endif
*
*     Check that the scales at which heavy quark masses have been defined
*     is above the value of the masses themselves.
*
      if(mass_scheme.eq."MSbar")then
         if(m2q(4).gt.q2th(4).or.
     1      m2q(5).gt.q2th(5).or.
     1      m2q(6).gt.q2th(6))then
            write(6,*) achar(27)//"[31mERROR:"
            write(6,*) "Each heavy quark mass reference scale must be"
            write(6,*) "larger than the value of the mass itself, i.e.:"
            write(6,*) "- Qc > mc(Qc)"
            write(6,*) "- Qb > mc(Qb)"
            write(6,*) "- Qt > mc(Qt)"
            write(6,*) "This condition is not fulfilled with:"
            write(6,*) "  "
            write(6,"(a,f8.3,a)") " mc(Qc) = ",dsqrt(m2q(4))," GeV"
            write(6,"(a,f8.3,a)") " mb(Qb) = ",dsqrt(m2q(5))," GeV"
            write(6,"(a,f8.3,a)") " mt(Qt) = ",dsqrt(m2q(6))," GeV"
            write(6,*) "   "
            write(6,"(a,f8.3,a)") " Qc = ",dsqrt(q2th(4))," GeV"
            write(6,"(a,f8.3,a)") " Qb = ",dsqrt(q2th(5))," GeV"
            write(6,"(a,f8.3,a)") " Qt = ",dsqrt(q2th(6))," GeV"
            write(6,*) achar(27)//"[0m"
            call exit(-10)
         endif
      endif
*
*     Check that the parameters of the LHA grid are within the maximum allowed.
*
      if(nq2LHA.gt.nq2max)then
         write(6,*) achar(27)//"[31mERROR:"
         write(6,*) "The number of points of the LHA Q2 grid exceeds",
     1              " the maximum:"
         write(6,*) "- input number =",nq2LHA
         write(6,*) "- maximum number =",nq2max
         write(6,*) achar(27)//"[0m"
         call exit(-10)
      endif
*
      if(nxLHA.gt.nxmax)then
         write(6,*) achar(27)//"[31mERROR:"
         write(6,*) "The number of points of the LHA x grid exceeds",
     1              " the maximum:"
         write(6,*) "- input number =",nxLHA
         write(6,*) "- maximum number =",nxmax
         write(6,*) achar(27)//"[0m"
         call exit(-10)
      endif
*
*     Make sure that time-like and polarized evolutions are not invoked
*     at the same time.
*
      if(Timelike.and.Polarized)then
         write(6,*) achar(27)//"[31mERROR:"
         write(6,*) "Time-like polarized evolution not available yet."
         write(6,*) "Switch off either the time-like or the polarized",
     1              " evolution."
         write(6,*) achar(27)//"[0m"
         call exit(-10)
      endif
c*
c*     Make sure that the time-like evolution at NNLO is done in the FFNS
c*
c      if(TimeLike.and.Evs.eq."VF".and.ipt.ge.2)then
c         write(6,*) achar(27)//"[31mERROR:"
c         write(6,*) "The time-like evolution at NNLO is available only",
c     1              " in the FFNS."
c         write(6,*) "(Unknown matching conditions)"
c         write(6,*) "Use 'SetFFNS(nf)' to use the FFNS."
c         write(6,*) achar(27)//"[0m"
c         call exit(-10)
c      endif
*
*     Make sure that the polarized evolution at NNLO is done in the FFNS
*
      if(Polarized.and.Evs.eq."VF".and.ipt.ge.2)then
         write(6,*) achar(27)//"[31mERROR:"
         write(6,*) "The polarized evolution at NNLO is available only",
     1              " in the FFNS."
         write(6,*) "(Unknown matching conditions)"
         write(6,*) "Use 'SetFFNS(nf)' to use the FFNS."
         write(6,*) achar(27)//"[0m"
         call exit(-10)
      endif
*
*     Security switches
*
*     If the fast evolution is enabled, disable automatically
*     the computation of the evolution operator
*
      if(EvolOp.and.FastEvol)then
         write(6,*) achar(27)//"[33m"//
     1              "WARNING: Computation of the evolution operator",
     2              " not possible if the fast evolution is enabled"
         write(6,*) "         ... disabling fast evolution"
     1              //achar(27)//"[0m"
         call SetFastEvolution(.false.)
      endif
*
*     If there is more than one subgrid and one of them is an external grid
*     the external evolution operator cannot be computed
*
      if(ngrid.gt.1.and.ThereAreExtGrids)then
         write(6,*) achar(27)//"[33m"//
     1              "WARNING: Computation of the evolution operator",
     2              " not possible if there is more than one subgrid"
         write(6,*) "             and one of the is external"
         write(6,*) "         ... disabling evolution operator",
     1              " computation"
     2              //achar(27)//"[0m"
         call EnableEvolutionOperator(.false.)
      endif
*
*     When the computation of the Evolution Operator is enabled
*     lock the grids by default (if there are no external grids).
*
      if(EvolOp.and..not.lock.and..not.ThereAreExtGrids)then
         write(6,*) achar(27)//"[33m"//
     1              "WARNING: computation of the evolution operator",
     2              " possible only on locked subgrids"
         write(6,*) "         ... locking subgrids"
     1              //achar(27)//"[0m"
         call LockGrids(.true.)
      endif
*
*     When enabling the evolution operator with the "QUniD" solution,
*     make use that the tau threshold is never crossed. Thus set the
*     tau mass to zero (infinity would do as well).
*
      if(EvolOp.and.Th.eq."QUniD")then
         write(6,*) achar(27)//"[33m"//
     1              "WARNING: no tau mass threshold crossing allowed ",
     2              "if the evolution operator is enabled"
         write(6,*) "         ... setting tau mass to zero"
     2              //achar(27)//"[0m"
         call SetTauMass(0d0)
      endif
*
*     If there are external grids the subgrids cannot be locked
*
      if(ThereAreExtGrids.and.lock)then
         write(6,*) achar(27)//"[33m"//
     1              "WARNING: if there are external grids they cannot ",
     2              "be locked"
         write(6,*) "         ... unlocking subgrids"
     1              //achar(27)//"[0m"
         call LockGrids(.false.)
      endif
*
*     If the polarized evolution is invoked at NNLO, give a warning
*     to inform the user that P^(2,v) is not known yet and that
*     the code uses P^(2,v) = P^(2,-).
*
      if(Polarized.and.ipt.ge.2)then
         write(6,*) achar(27)//"[33m"//
     1              "WARNING: the polarized evolution at NNLO is ",
     2              "incomplete."
         write(6,*) "         APFEL assumes P^(2,v) = P^(2,minus). "
     1              //achar(27)//"[0m"
      endif
*
*     Make sure that the NLO QED corrections are included only if
*     the 'QUniD' solution is used.
*
      if(NLOQED.and.Th(1:5).ne."QUniD")then
         write(6,*) achar(27)//"[33m"//
     1              "WARNING: NLO QED corrections can be included only"
         write(6,*) "         when using the 'QUniD' solution of the"
         write(6,*) "         DGLAP equation"
         write(6,*) "         ... turning off NLO QED corrections"
     1              //achar(27)//"[0m"
         call EnableNLOQEDCorrections(.false.)
      endif
*
*     Compute the heavy quark thresholds
*
      call ComputeHeavyQuarkThresholds
*
*     If the MSbar are used and the heavy quark masses are not given
*     at the mass scales themselves, compute the RG invariant masses
*     i.e. mc(mc), mb(mb) and mt(mt).
*
      if(mass_scheme.eq."MSbar") call ComputeRGInvariantMasses
*
*     If the alpha solution is "lambda", compute values of LambdaQCD
*     for all the number of flavours.     
*
      if(AlphaEvol(1:6).eq."lambda") call LambdaQCDnf
*
*     Compute alphas at the thresholds.
*     (Needed if the unified solution has been chosen with solution of
*     the DGLAP equation different from 'exactmu').
*
      call ThresholdAlphaQCD
*
*     Initialize HELL with the righ accuracies if the small-x
*     resummation is enabled
*
      if(Smallx) call initHELL
*
      return
      end
