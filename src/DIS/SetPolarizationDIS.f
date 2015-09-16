************************************************************************
*
*     SetPolarizationDIS.f:
*
*     This subroutine sets the mass scheme to be used for the computation
*     of the DIS observables.
*
************************************************************************
      subroutine SetPolarizationDIS(pol)
*
      implicit none
*
      include "../commons/PolarizationDIS.h"
*
*     Variables
*
      double precision pol
*
      PolarizationDIS = pol
      InPolarizationDIS = "done"
*
      return
      end
