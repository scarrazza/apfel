************************************************************************
*
*     SetGridParameters.f:
*
*     This subroutine sets the the parameters of the i-th x-space grid.
*
************************************************************************
      subroutine SetGridParameters(i,np,deg,x)
*
      implicit none
*
      include "../commons/grid.h"
*
*     Variables
*
      integer i,np,deg
      double precision x
      double precision eps
      parameter(eps=1d-12)
*
      nin(i)          = np
      inter_degree(i) = deg
      IsExt(i)        = .false.
      xmin(i)         = x
      if(i.eq.ngrid) xmin(i+1) = xmax + eps
*
      return
      end
