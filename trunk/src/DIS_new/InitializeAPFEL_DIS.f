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
*
*     Variables
*
      integer inf,nfi,nff
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
*     Evaluate DIS integrals on the grid
*
      call cpu_time(t1)
      do igrid=1,ngrid
         call initIntegralsDIS
      enddo
      call cpu_time(t2)
*
      if(Welcome)then
         write(6,"(a,a,f8.3,a)") " Initialization of the DIS module",
     1                           " completed in",t2-t1," s"
         write(6,*) " "
      endif
*
      return
      end
