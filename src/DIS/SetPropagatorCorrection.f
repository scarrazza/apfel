************************************************************************
*
*     SetPropagatorCorrection.f:
*
*     This subroutine sets the correction tothe Z propagator involved in
*     the computation of the NC DIS structure functions.
*
************************************************************************
      subroutine SetPropagatorCorrection(dr)
*
      implicit none
*
      include "../commons/PropagatorCorrection.h"
*
*     Variables
*
      double precision dr
*
      DeltaR   = dr
      InDeltaR = "done"
*
      return
      end
