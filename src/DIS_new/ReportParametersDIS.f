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
      include "../commons/SinThetaW.h"
      include "../commons/GFermi.h"
      include "../commons/CKM.h"
      include "../commons/TMC.h"
      include "../commons/DampingFONLL.h"
      include "../commons/SelectedCharge.h"
      include "../commons/krenQ.h"
      include "../commons/kfacQ.h"
      include "../commons/PropagatorCorrection.h"
      include "../commons/EWCouplings.h"
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
         write(6,"(a,f7.4)")    " sin(thetaW) =",SinThetaW
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
         write(6,"(a,f7.4)") " muR / Q = ",dsqrt(krenQ)
         write(6,"(a,f7.4)") " muF / Q = ",dsqrt(kfacQ)
*
         if(TMC)then
            write(6,*) "Target Mass corrections enabled"
         else
            write(6,*) "Target Mass corrections disabled"
         endif
*
         if(MassScheme(1:5).eq."FONLL")then
            if(DampingFONLL)then
               write(6,*) "FONLL damping factor enabled"
            else
               write(6,*) "FONLL damping factor disabled"
            endif
         endif
         write(6,*) achar(27)//"[0m"
      endif
*
      return
      end
