************************************************************************
*
*     ReportParametersDIS.f:
*
*     This subroutine reports the DIS parameters.
*
************************************************************************
      subroutine ReportParametersDIS
*
      implicit none
*
      include "../commons/Welcome.h"
      include "../commons/MassScheme.h"
      include "../commons/ProcessDIS.h"
      include "../commons/PolarizationDIS.h"
      include "../commons/ProjectileDIS.h"
      include "../commons/TargetDIS.h"
      include "../commons/ZedMass.h"
      include "../commons/WMass.h"
      include "../commons/ProtonMass.h"
      include "../commons/Sin2ThetaW.h"
      include "../commons/GFermi.h"
      include "../commons/CKM.h"
      include "../commons/TMC.h"
      include "../commons/DampingFONLL.h"
      include "../commons/SelectedCharge.h"
      include "../commons/krenQ.h"
      include "../commons/kfacQ.h"
      include "../commons/PropagatorCorrection.h"
      include "../commons/EWCouplings.h"
      include "../commons/DynScVar.h"
      include "../commons/IntrinsicCharm.h"
      include "../commons/NLOQEDCorrections.h"
      include "../commons/ScaleVariationProcedure.h"
*
*     Print welcome message and report of the parameters (if enabled)
*
      if(Welcome)then
         write(6,*) achar(27)//"[34m"//
     1              "Report of the electroweak parameters:"
         write(6,*) "  "
         write(6,"(a,f7.3,a)")  " Mass of the Z =",MZ," GeV"
         write(6,"(a,f7.3,a)")  " Mass of the W =",MW," GeV"
         write(6,"(a,f7.4,a)")  " Mass of the proton =",MProton," GeV"
         write(6,"(a,f7.4)")    " sin^2(thetaW) =",Sin2ThetaW
         write(6,"(a,es12.5)")   " GFermi =",GFermi
         write(6,"(a,3f7.4,a)") "       |",V_ud,V_us,V_ub," |"
         write(6,"(a,3f7.4,a)") " CKM = |",V_cd,V_cs,V_cb," |"
         write(6,"(a,3f7.4,a)") "       |",V_td,V_ts,V_tb," |"
         write(6,"(a,f7.5)") " Z propagator correction = ",DeltaR
         if(ExtCoup) write(6,*) "External EW couplings will be used"
         write(6,*) "  "
*
         write(6,*) "Report of the DIS parameters:"
         write(6,*) "  "
*
         write(6,"(a,a,a)") " Computation in the ",trim(MassScheme),
     1                      " mass scheme"
*
         if(ProcessDIS.eq."EM")then
            write(6,"(a,a,a)") " Electromagnetic (EM) process"
         elseif(ProcessDIS.eq."NC")then
            write(6,"(a,a,a)") " Neutral Current (NC) process"
         elseif(ProcessDIS.eq."CC")then
            write(6,"(a,a,a)") " Charged Current (CC) process"
         endif
*
         write(6,"(a,a,a,a)") " Scattering ",trim(ProjectileDIS),
     1                        " - ",TargetDIS
*
         if(PolarizationDIS.ne.0d0) 
     1   write(6,"(a,f7.3)") " Polarization fraction =",PolarizationDIS
*
         if(SelectedCharge(1:3).ne."all") 
     1   write(6,*) "Selected Charge: ",trim(SelectedCharge)
*
         if(DynScVAr)then
            write(6,*) "Dynamical scale variations enabled"
         else
            write(6,"(a,f7.4)") " muR / Q = ",dsqrt(krenQ)
            write(6,"(a,f7.4)") " muF / Q = ",dsqrt(kfacQ)
         endif
*
         if(ScVarProc.eq.1)then
            write(6,*) "No scale variations in the evolution"
         endif
*
         if(TMC)then
            write(6,*) "Target Mass corrections enabled"
         else
            write(6,*) "Target Mass corrections disabled"
         endif
*
         if(MassScheme(1:5).eq."FONLL")then
            if(DampingFONLL)then
               if(DampPowerFONLL(4).gt.0)then
                  write(6,"(a,a,i1)") " FONLL damping factor ",
     1                 "for charm enabled with suppression power = ",
     2                 DampPowerFONLL(4)
               elseif(DampPowerFONLL(4).eq.0)then
                  write(6,"(a,a,i1)") " FONLL damping factor ",
     1                 "for charm disbled"
               else
                  write(6,*) "Using BGMPU prescription for charm"
               endif
               if(DampPowerFONLL(5).gt.0)then
                  write(6,"(a,a,i1)") " FONLL damping factor ",
     1                 "for bottom enabled with suppression power = ",
     2                 DampPowerFONLL(5)
               elseif(DampPowerFONLL(5).eq.0)then
                  write(6,"(a,a,i1)") " FONLL damping factor ",
     1                 "for bottom disbled"
               else
                  write(6,*) "Using BGMPU prescription for bottom"
               endif
               if(DampPowerFONLL(6).gt.0)then
                  write(6,"(a,a,i1)") " FONLL damping factor ",
     1                 "for top enabled with suppression power = ",
     2                 DampPowerFONLL(6)
               elseif(DampPowerFONLL(6).eq.0)then
                  write(6,"(a,a,i1)") " FONLL damping factor ",
     1                 "for top disbled"
               else
                  write(6,*) "Using BGMPU prescription for top"
               endif
            else
               write(6,*) "FONLL damping factor disabled",
     1              " for all heavy quarks"
            endif
         endif
*
         if(IntrinsicCharm)then
            write(6,*) "Intrinsic charm enabled"
         else
            write(6,*) "Intrinsic charm disabled"
         endif
*
         if(.not.SFNLOQED)then
            write(6,*) "NLO QED Corrections disabled",
     1           " in the stucture functions"
         endif
         write(6,*) achar(27)//"[0m"
      endif
*
      return
      end
