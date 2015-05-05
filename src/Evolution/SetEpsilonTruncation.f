************************************************************************
*
*     SetEpsilonTruncation.f:
*
*     This subroutine sets the value of the small parameter epsilon
*     used for the truncation of the DGLAP solution.
*
************************************************************************
      subroutine SetEpsilonTruncation(eps)
*
      implicit none
*
      include "../commons/EpsTrunc.h"
*
*     Variables
*
      double precision eps
*
      EpsTrunc   = eps
      EpsEff     = 1d0
      InEpsTrunc = "done"
*
      return
      end
