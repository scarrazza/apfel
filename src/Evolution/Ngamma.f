************************************************************************
*
*     Ngamma.f:
*
*     This function returns the N-th Mellin moment of the photon PDF 
*     in the physical basis at the final scale.
*
************************************************************************
      function Ngamma(N)
*
      implicit none
*
      include "../commons/grid.h"
**
*     Input Variables
*
      integer N
**
*     Input Variables
*
      integer M
      double precision xgammawrap
      double precision dgauss,a,b,eps
      external xgammawrap

      common / gammaindex / M
**
*     Output Variables
*
      double precision Ngamma
*
      M = N
*
      a   = xmin(1)
      b   = xmax
      eps = 1d-7
      Ngamma = dgauss(xgammawrap,a,b,eps)
*
      return
      end
************************************************************************
*
*     Wrapping of the function xPDF
*
************************************************************************
      function xgammawrap(x)
*
      implicit none
**
*     Input variables
*
      double precision x
**
*     Internal variables
*
      integer M
      double precision xgamma

      common / gammaindex / M
**
*     Output variables
*
      double precision xgammawrap
*
      xgammawrap = x**(M-2) * xgamma(x)
*
      return
      end
