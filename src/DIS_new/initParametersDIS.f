************************************************************************
*
*     initParametersDIS.f:
*
*     It sets all the parameters for the computation of the DIS structure
*     functions if they were not set exernally before by the user.
*
************************************************************************
      subroutine initParametersDIS
*
      implicit none
*
      include "../commons/Welcome.h"
      include "../commons/MassScheme.h"
      include "../commons/ProcessDIS.h"
      include "../commons/PolarizationDIS.h"
      include "../commons/ProjectileDIS.h"
      include "../commons/TargetDIS.h"
      include "../commons/ipt.h"
      include "../commons/Nf_FF.h"
      include "../commons/ZedMass.h"
      include "../commons/WMass.h"
      include "../commons/ProtonMass.h"
      include "../commons/SinThetaW.h"
      include "../commons/GFermi.h"
      include "../commons/CKM.h"
      include "../commons/TMC.h"
*
*     Initialize default parameters (those that were not initialized before)
*
      if(InWelcome.ne."done")         call EnableWelcomeMessage(.true.)
      if(InMassScheme.ne."done")      call SetMassScheme("ZM-VFNS")
      if(InProcessDIS.ne."done")      call SetProcessDIS("EM")
      if(InPolarizationDIS.ne."done") call SetPolarizationDIS(0d0)
      if(InProjectileDIS.ne."done")   call SetProjectileDIS("electron")
      if(InTargetDIS.ne."done")       call SetTargetDIS("proton")
      if(InTMC.ne."done")      call EnableTargetMassCorrections(.false.)
*
      if(InMZ.ne."done")              call SetZMass(91.1876d0)
      if(InMW.ne."done")              call SetWMass(80.385d0)
      if(InMProton.ne."done")         call SetProtonMass(0.938272046d0)
      if(InSinThetaW.ne."done")       call SetSinThetaW(0.23126d0)
      if(InGFermi.ne."done")          call SetGFermi(1.1663787d-5)
      if(InCKM.ne."done")             call SetCKM(
     1                                0.97427d0, 0.22536d0, 0.00355d0,
     2                                0.22522d0, 0.97343d0, 0.04140d0,
     3                                0.00886d0, 0.04050d0, 0.99914d0)
*
*     Check the consistency of the input parameters
*
      write(6,*) "  "
      if(MassScheme.ne."ZM-VFNS".and.
     1   MassScheme(1:4).ne."FFNS".and.
     2   MassScheme(1:4).ne."FFN0".and.
     3   MassScheme.ne."FONLL-A".and.
     4   MassScheme.ne."FONLL-B".and.
     5   MassScheme.ne."FONLL-C")then
         write(6,*) "Mass scheme unknown:"
         write(6,*) "MassScheme = ",MassScheme
         write(6,*) "  "
         write(6,*) "The options are:"
         write(6,*) "- 'ZM-VFNS'"
         write(6,*) "- 'FFNS' (default NF=3)"
         write(6,*) "- 'FFNS3'"
         write(6,*) "- 'FFNS4'"
         write(6,*) "- 'FFNS5'"
         write(6,*) "- 'FFNS6'"
         write(6,*) "- 'FFN0' (default NF=3)"
         write(6,*) "- 'FFN03'"
         write(6,*) "- 'FFN04'"
         write(6,*) "- 'FFN05'"
         write(6,*) "- 'FFN06'"
         write(6,*) "- 'FONLL-A'"
         write(6,*) "- 'FONLL-B'"
         write(6,*) "- 'FONLL-C'"
         write(6,*) "  "
         call exit(-10)
      endif
*
      if(ProcessDIS.ne."EM".and.
     1   ProcessDIS.ne."NC".and.
     2   ProcessDIS.ne."CC")then
         write(6,*) "DIS process unknown:"
         write(6,*) "ProcessDIS = ",ProcessDIS
         write(6,*) "  "
         write(6,*) "The options are:"
         write(6,*) "- 'EM'"
         write(6,*) "- 'NC'"
         write(6,*) "- 'CC'"
         write(6,*) "  "
         call exit(-10)
      endif
*
      if(dabs(PolarizationDIS).gt.1d0)then
         write(6,*) "Polarization fraction not allowed:"
         write(6,*) "PolarizationDIS = ",PolarizationDIS
         write(6,*) "  "
         write(6,*) "PolarizationDIS must be between 1 and -1"
         write(6,*) "  "
         call exit(-10)
      endif
*
      if(ProjectileDIS(1:8).ne."electron".and.
     1   ProjectileDIS(1:8).ne."positron".and.
     2   ProjectileDIS(1:8).ne."neutrino".and.
     3   ProjectileDIS.ne."antineutrino")then
         write(6,*) "Projectile unknown:"
         write(6,*) "ProjectileDIS = ",ProjectileDIS
         write(6,*) "  "
         write(6,*) "The options are:"
         write(6,*) "- 'electron'"
         write(6,*) "- 'positron'"
         write(6,*) "- 'neutrino'"
         write(6,*) "- 'antineutrino'"
         write(6,*) "  "
         call exit(-10)
      endif
*
      if(TargetDIS(1:6).ne."proton".and.
     1   TargetDIS(1:7).ne."neutron".and.
     2   TargetDIS.ne."isoscalar".and.
     3   TargetDIS(1:4).ne."iron")then
         write(6,*) "Target unknown:"
         write(6,*) "TargetDIS = ",TargetDIS
         write(6,*) "  "
         write(6,*) "The options are:"
         write(6,*) "- 'proton'"
         write(6,*) "- 'neutron'"
         write(6,*) "- 'isoscalar'"
         write(6,*) "- 'iron'"
         write(6,*) "  "
         call exit(-10)
      endif
*
*     Additional settings
*
      if(MassScheme(1:4).eq."FFNS".or.
     1   MassScheme(1:4).eq."FFN0")then
         if(MassScheme(5:5).eq."3")then
            write(6,*) "INFO: Setting NF = 3 FFNS PDF evolution"
            call SetFFNS(3)
         elseif(MassScheme(5:5).eq."4")then
            write(6,*) "INFO: Setting NF = 4 FFNS PDF evolution"
            call SetFFNS(4)
         elseif(MassScheme(5:5).eq."5")then
            write(6,*) "INFO: Setting NF = 5 FFNS PDF evolution"
            call SetFFNS(5)
         elseif(MassScheme(5:5).eq."6")then
            write(6,*) "INFO: Setting NF = 6 FFNS PDF evolution"
            call SetFFNS(6)
         else
            write(6,*) "INFO: Setting NF = 3 FFNS PDF evolution"
            call SetFFNS(3)
         endif
      else
*     If the number of active flavour has not been specified (by means of SetFFNS)
*     set it automatically to 3.
         if(Nf_FF.lt.3.or.Nf_FF.gt.6) Nf_FF = 3
         call SetVFNS
         write(6,*) "INFO: Setting VFNS PDF evolution"
         if(MassScheme(1:5).eq."FONLL".and.ipt.eq.0)then
            call SetMassScheme("ZM-VFNS")
            write(6,*) "INFO: Any of the FONLL schemes at LO",
     1                 " concides with the the ZM-VFNS"
         endif


         if(MassScheme.eq."FONLL-A".and.ipt.eq.2)then
            call SetPerturbativeOrder(1)
            write(6,*) "INFO: For the FONLL-A scheme the perturbative",
     1                 " order has been automatically set to NLO"
         endif


         if(MassScheme.eq."FONLL-B".and.ipt.eq.2)then
            call SetMassScheme("FONLL-C")
            write(6,*) "INFO: The FONLL-B scheme at NNLO concides",
     1                 " with the the FONLL-C scheme"
         endif
         if(MassScheme.eq."FONLL-C".and.ipt.eq.1)then
            call SetMassScheme("FONLL-A")
            write(6,*) "INFO: The FONLL-C scheme at NLO concides",
     1                 " with the the FONLL-A scheme"
         endif
      endif
      write(6,*) " "
*
*     Print welcome message and report of the parameters (if enabled)
*
      if(Welcome)then
         write(6,*) "Report of the electroweak parameters:"
         write(6,*) "  "
         write(6,"(a,f7.3,a)")  " Mass of the Z =",MZ," GeV"
         write(6,"(a,f7.3,a)")  " Mass of the W =",MW," GeV"
         write(6,"(a,f7.4,a)")  " Mass of the proton =",MProton," GeV"
         write(6,"(a,f7.4)")    " sin(thetaW) =",SinThetaW
         write(6,"(a,es12.5)")   " GFermi =",GFermi
         write(6,"(a,3f7.4,a)") "       |",V_ud,V_us,V_ub," |"
         write(6,"(a,3f7.4,a)") " CKM = |",V_cd,V_cs,V_cb," |"
         write(6,"(a,3f7.4,a)") "       |",V_td,V_ts,V_tb," |"
         write(6,*) "  "
*
c$$$         write(6,*) "Report of the DIS parameters:"
c$$$         write(6,*) "  "
*
         write(6,"(a,a,a)") " Computation in the ",MassScheme,
     1                      " mass scheme"
c$$$*
c$$$         if(ProcessDIS.eq."EM")then
c$$$            write(6,"(a,a,a)") " Electromagnetic (EM) process"
c$$$         elseif(ProcessDIS.eq."NC")then
c$$$            write(6,"(a,a,a)") " Neutral Current (NC) process"
c$$$         elseif(ProcessDIS.eq."CC")then
c$$$            write(6,"(a,a,a)") " Charged Current (CC) process"
c$$$         endif
c$$$*
c$$$         if(ProjectileDIS(1:8).eq."electron".or.
c$$$     1      ProjectileDIS(1:8).eq."positron".or.
c$$$     2      ProjectileDIS(1:8).eq."neutrino")then
c$$$            write(6,"(a,a,a,a)") " Scattering ",
c$$$     1                           ProjectileDIS(1:8),
c$$$     2                           " - ",TargetDIS
c$$$         else
c$$$            write(6,"(a,a,a,a)") " Scattering ",
c$$$     1                           ProjectileDIS,
c$$$     3                           " - ",TargetDIS
c$$$         endif
c$$$*
c$$$         if(PolarizationDIS.ne.0d0) 
c$$$     1   write(6,"(a,f7.3)") " Polarization fraction =",PolarizationDIS
c$$$*
         if(TMC)then
            write(6,*) "Target Mass corrections enabled"
         else
            write(6,*) "Target Mass corrections disabled"
         endif
      endif
*
      return
      end
