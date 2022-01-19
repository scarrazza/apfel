************************************************************************
*
*     AlphaQCD.f:
*
*     This function returns the value of alpha_QCD at the given scale
*     using the parameters of the evolution.
*     Be careful because, if kren.ne.1, the matching is not done at the
*     heavy quark thresholds.
*
************************************************************************
      function AlphaQCD(Q)
*
      implicit none
*
      include "../commons/consts.h"
**
*     Input Variables
*
      double precision Q
**
*     Internal Variables
*
      double precision a_QCD
**
*     Output Variables
*
      double precision AlphaQCD
*
      AlphaQCD = 4d0 * pi * a_QCD(Q * Q)
*
      return
      end
