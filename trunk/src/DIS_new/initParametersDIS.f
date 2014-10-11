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
*
*     Initialize default parameters (those that were not initialized before)
*
      if(InMassScheme.ne."done")      call SetMassScheme("ZM-VFNS")   ! "ZM-VFNS", "FFNS", "FONLL-A", "FONLL-B" or "FONLL-C"
      if(InProcessDIS.ne."done")      call SetProcessDIS("EM")        ! "EM", "NC" or "CC"
      if(InPolarizationDIS.ne."done") call SetPolarizationDIS(0d0)
*
*     Check the consistency of the input parameters
*
      if(MassScheme.ne."ZM-VFNS".and.
     1   MassScheme(1:4).ne."FFNS".and.
     2   MassScheme.ne."FONLL-A".and.
     3   MassScheme.ne."FONLL-B".and.
     4   MassScheme.ne."FONLL-C")then
         write(6,*) "Mass scheme unknown:"
         write(6,*) "MassScheme = ",MassScheme
         write(6,*) "  "
         write(6,*) "The options are:"
         write(6,*) "- 'ZM-VFNS'"
         write(6,*) "- 'FFNS'"
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
*     Print welcome message and report of the parameters (if enabled)
*
      if(Welcome)then
         write(6,*) "Report of the DIS parameters:"
         write(6,*) "  "
*
         write(6,"(a,a,a)") " Computation in the ",MassScheme,
     1                      " mass scheme"
*
         if(ProcessDIS.eq."EM")then
            write(6,"(a,a,a)") " Electromagnetic process"
         elseif(ProcessDIS.eq."NC")then
            write(6,"(a,a,a)") " Neutral current process"
         elseif(ProcessDIS.eq."CC")then
            write(6,"(a,a,a)") " Charged current process"
         endif
*
         write(6,"(a,f7.3)") " Polarization fraction =",PolarizationDIS
*
         write(6,*) " "
      endif
*
      return
      end
