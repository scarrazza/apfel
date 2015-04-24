************************************************************************
*
*     SetLHgridParameters.f:
*
*     This subroutine sets the the parameters of grid over which
*     the LHAPDF set will be tabulated.
*
************************************************************************
      subroutine SetLHgridParameters(nx,nxm,xmin,xm,xmax,
     1                               nq2,q2min,q2max)
*
      implicit none
*
      include "../commons/LHAgrid.h"
*
*     Variables
*
      integer nx,nxm,nq2
      double precision q2min,q2max
      double precision xmin,xm,xmax
*
      nxLHA    = nx
      nxmLHA   = nxm
      xminLHA  = xmin
      xmLHA    = xm
      xmaxLHA  = xmax
*
      nq2LHA   = nq2
      q2minLHA = q2min
      q2maxLHA = q2max
*
      InLHgrid = "done"
*
      return
      end
