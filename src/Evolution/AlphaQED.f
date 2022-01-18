************************************************************************
*
*     AlphaQED.f:
*
*     This function returns the value of alpha_QED at the given scale
*     using the parameters of the evolution.
*     Be careful because, if kren.ne.1, the matching is not done at the
*     heavy quark thresholds.
*
************************************************************************
      function AlphaQED(Q)
*
      implicit none
*
      include "../commons/consts.h"
      include "../commons/krenalpha.h"
**
*     Input Variables
*
      double precision Q
**
*     Internal Variables
*
      double precision Q2
      double precision a_QED
**
*     Output Variables
*
      double precision AlphaQED
*
      Q2 = Q * Q / krena
*
      AlphaQED = 4d0 * pi * a_QED(Q2)
*
      return
      end
