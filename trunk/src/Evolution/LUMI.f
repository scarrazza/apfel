************************************************************************
*
*     LUMI.f:
*
*     This function returns the luminosity for i,j = [-6,7] where
*     i = [-6,6] \ {0} are quarks, 0 gluon and 7 photon. 
*     S the square of the CM energy.
*
************************************************************************
      function LUMI(i,j,S)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/scales.h"
**
*     Input Variables
*
      integer i,j
      double precision S, MX
**
*     Internal Variables
*      
      double precision dgauss,a,b,eps,tau
      double precision LUMIwrap
      external LUMIwrap

      integer ix,jx
      double precision taux
      common / LUMIindex / ix,jx,taux      
**
*     Output Variables
*
      double precision LUMI
*
      if (i.gt.7 .or. i.lt.-6) then
         write(6,*) "LUMI.f Index out of range"
         call exit(-10)
      endif

      if (j.gt.7 .or. j.lt.-6) then
         write(6,*) "LUMI.f Index out of range"
         call exit(-10)
      endif

      ix = i
      jx = j
      MX = dsqrt(Q2fin)
      tau = MX * MX / S
      taux = tau
*      
      a   = tau
      b   = 1d0
      eps = 1d-5
      LUMI = dgauss(LUMIwrap,a,b,eps) / S
*
      return
      end
************************************************************************
*
*     Wrapping of the function xPDF
*
************************************************************************
      function LUMIwrap(x)
*
      implicit none
**
*     Input variables
*
      double precision x
**
*     Internal variables
*
      integer ix,jx
      double precision taux,X1,X2,xPDF,xgamma
      common / LUMIindex / ix,jx,taux
**
*     Output variables
*
      double precision LUMIwrap,xpdf1,xpdf2
*
      X1 = x
      X2 = taux / X1     

      if (ix.eq.7) then
         xpdf1 = xgamma(X1)
      else
         xpdf1 = xPDF(ix,X1)
      endif

      if (jx.eq.7) then
         xpdf2 = xgamma(X2)
      else
         xpdf2 = xPDF(jx,X2)
      endif         

      LUMIwrap = (1d0/X1)*(xpdf1/X1)*(xpdf2/X2)
*
      return
      end
