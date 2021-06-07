************************************************************************
*
*     InitializeAPFEL_DIS.f:
*
*     This routine initializes the integrals needed to compute the 
*     structure functions.
*
************************************************************************
      subroutine InitializeAPFEL_DIS
*
      implicit none
*
      include "../commons/Welcome.h"
      include "../commons/grid.h"
      include "../commons/TimeLike.h"
      include "../commons/Smallx.h"
      include "../commons/Polarized.h"
      Include "../commons/InAPFELDIS.h"
*
*     Variables
*
      double precision t1,t2
*
*     Read input parameters
*
      call initParametersDIS
*
*     Initialize the APFEL evolution
*
      call InitializeAPFEL
*
*     Report DIS parameters
*
      call ReportParametersDIS
*
*     Evaluate DIS or SIA integrals on the grid
*
      call cpu_time(t1)
      if(TimeLike)then
         do igrid=1,ngrid
            call initIntegralsSIA
         enddo
      else
         do igrid=1,ngrid
            if(Polarized)then
               call initIntegralspDIS
            else
               call initIntegralsDIS
            endif
            if(Smallx) call initIntegralsDISRes
         enddo
      endif
      call cpu_time(t2)
*
      if(Welcome)then
         write(6,"(a,a,f8.3,a)") " Initialization of the DIS module",
     1                           " completed in",t2-t1," s"
         write(6,*) " "
      endif
*
*     Initialization of the DIS module complete
*
      InAPFELDIS = "done"
*
      return
      end
