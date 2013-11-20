************************************************************************
*
*     SetQLimits.f:
*
*     This subroutine sets the minimum and the maximum energy allowed
*     for the evolution.
*
************************************************************************
      subroutine SetQLimits(Qmin,Qmax)
*
      implicit none
*
      include "../commons/scales.h"
*
*     Variables
*
      double precision Qmin,Qmax
*
      Q2min = Qmin * Qmin
      Q2max = Qmax * Qmax
      InScales = "done"
*
      return
      end
