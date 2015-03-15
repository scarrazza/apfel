************************************************************************
*
*     NEUTRAL CURRENT
*
************************************************************************
*
*     MassiveCoefficientFunctions.f:
*
*     This file contains all the ZM coefficient functions.
*
*
*     O(as):
*
*        - icoef = 1  F_2,g (NLO, gluon)
*        - icoef = 2  F_L,g (NLO, gluon)
*
*     O(as^2):
*
*        - icoef = 3  F_2,g (NNLO, gluon)
*        - icoef = 4  F_2,q (NNLO, singlet)
*        - icoef = 5  F_L,g (NNLO, gluon)
*        - icoef = 6  F_L,q (NNLO, singlet)
*
*        - icoef = 7  F_2,q (NNLO, non-singlet for the light structure function)
*        - icoef = 8  F_L,q (NNLO, non-singlet for the light structure function)
*
*     Coefficients of the terms proportional to ln(q2/m2):
*
*        - icoef = 9  F_2,g (NNLO, gluon, bar term)
*        - icoef = 10 F_2,q (NNLO, singlet, bar term)
*        - icoef = 11 F_L,g (NNLO, gluon, bar term)
*        - icoef = 12 F_L,q (NNLO, singlet, bar term)
*
*        - icoef = 13 F_2,q (NNLO, non-singlet for the light structure function, bar term)
*        - icoef = xx F_L,q (NNLO, non-singlet for the light structure function, bar term) = 0
*
*     "ixi" runs on the xi grid where xi = q2 / m2. 
*
*     Phys.Lett.B594:299-307,2004
*     e-Print: hep-ph/0404034
*
************************************************************************
      function MassiveCF(icoef,ixi,z)
*
      implicit none
*
      include "../commons/consts.h"
      include "../commons/coeffhqmellin.h"
      include "../commons/ColorFactors.h"
      include "../commons/mass_scheme.h"
**
*     input variables
*
      integer icoef,ixi
      double precision z
**
*     internal variables
*
      integer i
      double precision rho
      double precision cm21g,cml1g
      double precision cm22q,cml2q
      double precision dcm21g,dcml1g
      double precision h1
**
*     output variables
*
      double precision MassiveCF
*
      rho = 1d0 / ( 4d0 / xigrid(ixi) + 1d0 )
*
      MassiveCF = 0d0
      if(z.ge.rho) return
*
*     Exact expressions
*
      if(icoef.eq.1)then
         MassiveCF = 2d0 * cm21g(xigrid(ixi),z)
         return
      elseif(icoef.eq.2)then
         MassiveCF = 2d0 * cml1g(xigrid(ixi),z)
         return
      elseif(icoef.eq.7)then
         MassiveCF = cm22q(xigrid(ixi),z)
         return
      elseif(icoef.eq.8)then
         MassiveCF = cml2q(xigrid(ixi),z)
         return
      endif
*
      do i=0,m_coef(icoef)-1
         MassiveCF = MassiveCF + coef(icoef,ixi,i+1) * z**i
      enddo
      MassiveCF = 16d0 * pi * xigrid(ixi) 
     1          * z**(-coef_p1(icoef)-1)
     2          * ( rho - z )**(-coef_p2(icoef))
     3          * MassiveCF
*
      if(icoef.eq.1.or.icoef.eq.2) MassiveCF = MassiveCF / 16d0 / pi**2
*
*     If the MSbar masses are to be used, add the appropriate term to
*     the NNLO gluon coefficient functions.
*
      if(mass_scheme.eq."MSbar")then
         if(icoef.eq.3)then
            h1 = CF * ( 4d0 + dlog(xigrid(ixi)) )
            MassiveCF = MassiveCF + h1 * dcm21g(xigrid(ixi),z)
         elseif(icoef.eq.5)then
            h1 = CF * ( 4d0 + dlog(xigrid(ixi)) )
            MassiveCF = MassiveCF + h1 * dcml1g(xigrid(ixi),z)
         endif
      endif
*
      return
      end
*
************************************************************************
*
*     Order alphas gluon coefficient functions for F2 and FL in the 
*     massive scheme.
*
************************************************************************
      function cm21g(xi,z)
*
      implicit none
**
*     Input Variables
*
      double precision xi,z
**
*     Internal Variables
*
      double precision zmax
      double precision eps,v
**
*     Output Variables
*
      double precision cm21g
*
      zmax = 1d0 / ( 1d0 + 4d0 / xi )
*
      cm21g = 0d0
      if(z.gt.zmax) return
*
      eps = 1d0 / xi
      v   = dsqrt( 1d0 - 4d0 * z / ( 1d0 - z ) / xi )
*
      cm21g = ( z**2d0 + ( 1d0 - z )**2d0 
     1      + 4d0 * eps * z * ( 1d0 - 3d0 * z ) 
     2      - 8d0 * eps**2d0 * z**2d0 ) 
     3      * dlog( ( 1d0 + v ) / ( 1d0 - v ) )
     4      + ( 8d0 * z * ( 1d0 - z ) - 1d0 
     5      - 4d0 * eps * z * ( 1d0 - z ) ) * v
*
      return
      end
*
************************************************************************
*
*     Logarithmic derivative with rescpect to the mass M of the function
*     "cm21g": dcm21g = M * d(cm21g)/dM.
*     (Needed for the MSbar mass implementation).
*
************************************************************************
      function dcm21g(xi,z)
*
      implicit none
**
*     Input Variables
*
      double precision xi,z
**
*     Internal Variables
*
      double precision zmax
      double precision eps,v
**
*     Output Variables
*
      double precision dcm21g
*
      zmax = 1d0 / ( 1d0 + 4d0 / xi )
*
      dcm21g = 0d0
      if(z.gt.zmax) return
*
      eps = 1d0 / xi
      v   = dsqrt( 1d0 - 4d0 * z / ( 1d0 - z ) / xi )
*
      dcm21g = 4d0 * eps * ( ( - 6d0 - 8d0 * eps ) * z**2d0 + 2d0 * z )
     1       * dlog( ( 1d0 + v ) / ( 1d0 - v ) ) 
     2       + 8d0 * eps * z * ( z - 1d0 ) * v 
     3       - 2d0 * ( ( - 8d0 * eps**2d0 - 12d0 * eps + 2d0 ) * z**2d0
     4       - ( 2d0 - 4d0 * eps ) * z + 1d0 ) / v 
     5       + 4d0 * eps * z * ( - 4d0 * ( 2d0 - eps ) * z**2d0
     6       + 4d0 * ( 2d0 - eps ) * z - 1d0 ) / v / ( z - 1d0 )
*
      return
      end
*
************************************************************************
      function cml1g(xi,z)
*
      implicit none
**
*     Input Variables
*
      double precision xi,z
**
*     Internal Variables
*
      double precision zmax
      double precision eps,v
**
*     Output Variables
*
      double precision cml1g
*
      zmax = 1d0 / ( 1d0 + 4d0 / xi )
*
      cml1g = 0d0
      if(z.gt.zmax) return
*
      eps = 1d0 / xi
      v   = dsqrt( 1d0 - 4d0 * z / ( 1d0 - z ) / xi )
*
      cml1g = - 8d0 * eps * z**2d0 
     1      * dlog( ( 1d0 + v ) / ( 1d0 - v ) )
     2      + 4d0 * v * z * ( 1d0 - z )
*
      return
      end
*
************************************************************************
*
*     Logarithmic derivative with rescpect to the mass M of the function
*     "cml1g": dcml1g = M * d(cml1g)/dM.
*     (Needed for the MSbar mass implementation).
*
************************************************************************
      function dcml1g(xi,z)
*
      implicit none
**
*     Input Variables
*
      double precision xi,z
**
*     Internal Variables
*
      double precision zmax
      double precision eps,v
**
*     Output Variables
*
      double precision dcml1g
*
      zmax = 1d0 / ( 1d0 + 4d0 / xi )
*
      dcml1g = 0d0
      if(z.gt.zmax) return
*
      eps = 1d0 / xi
      v   = dsqrt( 1d0 - 4d0 * z / ( 1d0 - z ) / xi )
*
      dcml1g = - 16d0 * eps * z**2d0 * dlog( ( 1d0 + v ) / ( 1d0 - v ) ) 
     1       + 16d0 * eps * z**2d0 / v 
     2       + 4d0 * eps * z * ( - 4d0 * z**2d0 + 4d0 * z ) 
     3       / v / ( z - 1d0 )
*
      return
      end
*
************************************************************************
*
*     Order alphas^2 pure-singlet coefficient functions for F2 and FL
*     in the massive scheme.
*     Reference: Appendix A of hep-ph/9601302
*
************************************************************************
      function cm22q(xi,z)
*
      implicit none
**
*     Input Variables
*
      double precision xi,z
**
*     Internal Variables
*
      double precision zmax
      double precision sq1,sq2,l1,l2,l3,dil1,dil2,dil3,dil4
      double precision ddilog
**
*     Output Variables
*
      double precision cm22q
*
      zmax = 1d0 / ( 1d0 + 4d0 / xi )
*
      cm22q = 0d0
      if(z.ge.zmax) return
*
      sq1 = dsqrt( 1d0 - 4d0 * z / xi / ( 1d0 - z ) )
      sq2 = dsqrt( 1d0 - 4d0 * z / xi )
*
      l1 = dlog( ( 1d0 + sq1 ) / ( 1d0 - sq1 ) )
      l2 = dlog( ( 1d0 + sq2 ) / ( 1d0 - sq2 ) )
      l3 = dlog( ( sq2 + sq1 ) / ( sq2 - sq1 ) )
*
      dil1 = ddilog( ( 1d0 - z ) * ( 1d0 + sq1 ) / ( 1d0 + sq2 ) )
      dil2 = ddilog( ( 1d0 - sq2 ) / ( 1d0 + sq1 ) )
      dil3 = ddilog( ( 1d0 - sq1 ) / ( 1d0 + sq2 ) )
      dil4 = ddilog( ( 1d0 + sq1 ) / ( 1d0 + sq2 ) )
*
      cm22q =2d0/3d0*((4d0/3d0*(1d0+z**2d0)/(1d0-z)
     1     -16d0/(1d0-z)*(z/xi)**2d0*(1d0-9d0*z+9d0*z**2d0))
     2     *(dlog((1d0-z)/z**2d0)*l1+l1*l2+2d0*(-dil1+dil2+dil3-dil4))
     3     +(-8d0/3d0+4d0/(1d0-z)+(z/(1d0-z)/xi)**2d0*(128d0-432d0*z
     4     +288d0*z**2d0-8d0/(1d0-z)))*l1+(88d0/9d0+136d0/9d0*z
     5     -152d0/9d0/(1d0-z)+(z/(1d0-z)/xi)*(464d0/9d0-512d0/3d0*z
     6     +2048d0/9d0*z**2d0)+(z/(1d0-z)/xi)**2d0*(-832d0/9d0
     7     +6208d0/9d0*z-11392d0/9d0*z**2d0+6016d0/9d0*z**3d0))*l3/sq2
     8     +(-272d0/27d0-1244d0/27d0*z+718d0/27d0/(1d0-z)+(z/(1d0-z)/xi)
     9     *(-3424d0/27d0+15608d0/27d0*z-4304d0/9d0*z**2d0
     1     +20d0/27d0/(1d0-z)))*sq1)
*
      return
      end
*
************************************************************************
      function cml2q(xi,z)
*
      implicit none
**
*     Input Variables
*
      double precision xi,z
**
*     Internal Variables
*
      double precision zmax
      double precision sq1,sq2,l1,l2,l3,dil1,dil2,dil3,dil4
      double precision ddilog
**
*     Output Variables
*
      double precision cml2q
*
      zmax = 1d0 / ( 1d0 + 4d0 / xi )
*
      cml2q = 0d0
      if(z.ge.zmax) return
*
      sq1 = dsqrt( 1d0 - 4d0 * z / xi / ( 1d0 - z ) )
      sq2 = dsqrt( 1d0 - 4d0 * z / xi )
*
      l1 = dlog( ( 1d0 + sq1 ) / ( 1d0 - sq1 ) )
      l2 = dlog( ( 1d0 + sq2 ) / ( 1d0 - sq2 ) )
      l3 = dlog( ( sq2 + sq1 ) / ( sq2 - sq1 ) )
*
      dil1 = ddilog( ( 1d0 - z ) * ( 1d0 + sq1 ) / ( 1d0 + sq2 ) )
      dil2 = ddilog( ( 1d0 - sq2 ) / ( 1d0 + sq1 ) )
      dil3 = ddilog( ( 1d0 - sq1 ) / ( 1d0 + sq2 ) )
      dil4 = ddilog( ( 1d0 + sq1 ) / ( 1d0 + sq2 ) )
*
      cml2q =2d0/3d0*(96d0*z*(z/xi)**2d0*(dlog((1d0-z)/z**2d0)*l1
     1     +l1*l2+2d0*(-dil1+dil2+dil3-dil4))+(z/xi/(1d0-z))**2d0
     2     *(64d0-288d0*z+192d0*z**2d0)*l1+z*(16d0/3d0-416d0*z/3d0/xi
     3     +1408d0*z**2d0/3d0/xi**2d0)*l3/sq2+(16d0/3d0-400d0*z/18d0
     4     +z*(-160d0/3d0+3824d0*z/9d0-992d0*z**2d0/3d0)/(1d0-z)/xi)
     5     *sq1)
*
      return
      end
*
************************************************************************
*
*     Evaluation of the integral in eq. (97) of ArXiv.1001.2312 needed 
*     for the Adler sum rule to be fulfilled.
*
************************************************************************
      function cm22q_adler(z)
*
      implicit none
*
      include "../commons/wrapDIS.h"
      include "../commons/coeffhqmellin.h"
**
*     Input Variables
*
      double precision z
**
*     Input Variables
*
      double precision xi
      double precision cm22q
**
*     Output Variables
*
      double precision cm22q_adler
*
      xi = xigrid(wixi)
      cm22q_adler = cm22q(xi,z)
*
      return
      end
*
************************************************************************
*
*     CHARGED CURRENT
*
************************************************************************
*
*     Order alphas^1 coeficient functions (NLO)
*     Expansion parameter alphas/4*pi
*
************************************************************************
c$$$      function c2ns1cca(xi,z)
c$$$*
c$$$      implicit none
c$$$*
c$$$      include "../commons/ColorFactors.h"
c$$$**
c$$$*     Input Variables
c$$$*
c$$$      double precision xi,z
c$$$**
c$$$*     Internal Variables
c$$$*
c$$$      double precision lambda
c$$$      double precision kQF2
c$$$**
c$$$*     Output Variables
c$$$*
c$$$      double precision c2ns1cca
c$$$*
c$$$      kQF2 = 1d0                ! Q2 / muF2
c$$$      lambda = xi / ( 1d0 + xi )
c$$$*
c$$$      c2ns1cca = 2d0 * CF * ( - ( 1d0 + z**2d0 ) * dlog(z) / ( 1d0 - z )
c$$$     1         + ( 2d0 - dlog( kQF2 / lambda ) 
c$$$     2         -  2d0 * dlog( 1d0 - z ) + dlog( 1d0 - lambda * z ) ) 
c$$$     3         * ( 1d0 + z ) + 1d0 / lambda )
c$$$*
c$$$      return
c$$$      end
c$$$*
c$$$************************************************************************
c$$$      function c2ns1ccb(xi,z)
c$$$*
c$$$      implicit none
c$$$*
c$$$      include "../commons/ColorFactors.h"
c$$$**
c$$$*     Input Variables
c$$$*
c$$$      double precision xi,z
c$$$**
c$$$*     Internal Variables
c$$$*
c$$$      double precision lambda
c$$$      double precision kQF2
c$$$**
c$$$*     Output Variables
c$$$*
c$$$      double precision c2ns1ccb
c$$$*
c$$$      kQF2 = 1d0                ! Q2 / muF2
c$$$      lambda = xi / ( 1d0 + xi )
c$$$*
c$$$      c2ns1ccb = 2d0 * CF * ( 2d0 * ( 2d0 * dlog( 1d0 - z ) 
c$$$     1         - dlog( 1d0 - lambda * z ) ) / ( 1d0 - z ) 
c$$$     2         + 2d0 * ( - 1d0 + dlog( kQF2 / lambda ) ) / ( 1d0 - z )
c$$$     3         + ( 2d0 * lambda**2d0 - lambda - 1d0 ) / lambda 
c$$$     4         / ( 1d0 - lambda * z ) + ( 1d0 - z ) 
c$$$     5         / ( 1d0 - lambda * z )**2d0 / 2d0 )
c$$$*
c$$$      return
c$$$      end
c$$$*
c$$$************************************************************************
c$$$      function c2ns1ccc(Rf,xi,z)
c$$$*
c$$$      implicit none
c$$$*
c$$$      include "../commons/ColorFactors.h"
c$$$      include "../commons/consts.h"
c$$$**
c$$$*     Input Variables
c$$$*
c$$$      double precision Rf,xi,z
c$$$**
c$$$*     Internal Variables
c$$$*
c$$$      double precision lambda
c$$$      double precision KA,ln1mz,ln1mlz
c$$$      double precision kQF2
c$$$**
c$$$*     Output Variables
c$$$*
c$$$      double precision c2ns1ccc
c$$$*
c$$$      kQF2 = 1d0                ! Q2 / muF2
c$$$      lambda = xi / ( 1d0 + xi )
c$$$      KA = ( 1d0 - lambda ) * dlog( 1d0 - lambda ) / lambda
c$$$      ln1mz  = dlog( 1d0 - z )
c$$$      ln1mlz = dlog( 1d0 - lambda * z )
c$$$*
c$$$      c2ns1ccc = 2d0 * CF * ( - 4d0 - 1d0 / 2d0 / lambda - 2d0 * zeta2 
c$$$     1         - ( 1d0 + lambda ) * KA / 2d0 / lambda 
c$$$     2         + 3d0 * dlog( kQF2 / lambda ) / 2d0
c$$$     3         + 2d0 * ln1mz**2d0 - 2d0 * Rf 
c$$$     4         + 2d0 * ( - 1d0 + dlog( kQF2 / lambda ) ) * ln1mz
c$$$     5         + ( 2d0 * lambda**2d0 - lambda - 1d0 / 2d0 ) 
c$$$     6         * ln1mlz / lambda**2d0 + ( 1d0 - lambda ) * z 
c$$$     7         / 2d0 / lambda / ( 1d0 - lambda * z ) )
c$$$*
c$$$      return
c$$$      end
c$$$*
c$$$************************************************************************
c$$$      function c2g1cca(xi,z)
c$$$*
c$$$      implicit none
c$$$*
c$$$      include "../commons/ColorFactors.h"
c$$$**
c$$$*     Input Variables
c$$$*
c$$$      double precision xi,z
c$$$**
c$$$*     Internal Variables
c$$$*
c$$$      double precision lambda
c$$$      double precision kQF2
c$$$**
c$$$*     Output Variables
c$$$*
c$$$      double precision c2g1cca
c$$$*
c$$$      kQF2 = 1d0                ! Q2 / muF2
c$$$      lambda = xi / ( 1d0 + xi )
c$$$*
c$$$      c2g1cca = 4d0 * TR * ( ( z**2d0 + ( 1d0 - z )**2d0 )
c$$$     1        * ( dlog( ( 1d0 -  z ) / z ) - dlog( 1d0 - lambda ) / 2d0 
c$$$     2        + dlog( kQF2 / lambda ) / 2d0 ) + 8d0 * z * ( 1d0 - z )
c$$$     3        - 1d0 + ( 1d0 - lambda ) * ( - 6d0 * ( 1d0 + 2d0*lambda )
c$$$     4        * z * ( 1d0 - z ) + 1d0 / ( 1d0 - lambda * z ) 
c$$$     5        + 6d0 * lambda * z * ( 1d0 - 2d0 * lambda * z ) 
c$$$     6        * dlog( ( 1d0 - lambda * z ) / ( 1d0 - lambda ) / z ) ) )
c$$$*
c$$$      return
c$$$      end
c$$$*
c$$$************************************************************************
c$$$      function clns1cca(xi,z)
c$$$*
c$$$      implicit none
c$$$*
c$$$      include "../commons/ColorFactors.h"
c$$$**
c$$$*     Input Variables
c$$$*
c$$$      double precision xi,z
c$$$**
c$$$*     Internal Variables
c$$$*
c$$$      double precision lambda
c$$$      double precision kQF2
c$$$**
c$$$*     Output Variables
c$$$*
c$$$      double precision clns1cca
c$$$*
c$$$      kQF2 = 1d0                ! Q2 / muF2
c$$$      lambda = xi / ( 1d0 + xi )
c$$$*
c$$$      clns1cca = 2d0 * CF * ( ( - ( 1d0 + z**2d0 ) * dlog(z) 
c$$$     1         / ( 1d0 - z ) + ( - dlog( kQF2 / lambda ) 
c$$$     2         -  2d0 * dlog( 1d0 - z ) + dlog( 1d0 - lambda * z ) ) 
c$$$     3         * ( 1d0 + z ) + 3d0 ) * ( 1d0 - lambda ) 
c$$$     4         + ( 1d0 + lambda ) * z )
c$$$*
c$$$      return
c$$$      end
c$$$*
c$$$************************************************************************
c$$$      function clns1ccb(xi,z)
c$$$*
c$$$      implicit none
c$$$*
c$$$      include "../commons/ColorFactors.h"
c$$$**
c$$$*     Input Variables
c$$$*
c$$$      double precision xi,z
c$$$**
c$$$*     Internal Variables
c$$$*
c$$$      double precision lambda
c$$$      double precision kQF2
c$$$**
c$$$*     Output Variables
c$$$*
c$$$      double precision clns1ccb
c$$$*
c$$$      kQF2 = 1d0                ! Q2 / muF2
c$$$      lambda = xi / ( 1d0 + xi )
c$$$*
c$$$      clns1ccb = 2d0 * CF * ( 2d0 * ( 2d0 * dlog( 1d0 - z ) 
c$$$     1         - dlog( 1d0 - lambda * z ) ) / ( 1d0 - z )
c$$$     2         + 2d0 * ( - 1d0 + dlog( kQF2 / lambda ) ) / ( 1d0 - z ) 
c$$$     3         - 2d0 / ( 1d0 - lambda * z ) + ( 1d0 - z ) 
c$$$     4         / ( 1d0 - lambda * z )**2d0 / 2d0 ) * ( 1d0 - lambda )
c$$$*
c$$$      return
c$$$      end
c$$$*
c$$$************************************************************************
c$$$      function clns1ccc(Rf,xi,z)
c$$$*
c$$$      implicit none
c$$$*
c$$$      include "../commons/ColorFactors.h"
c$$$      include "../commons/consts.h"
c$$$**
c$$$*     Input Variables
c$$$*
c$$$      double precision Rf,xi,z
c$$$**
c$$$*     Internal Variables
c$$$*
c$$$      double precision lambda
c$$$      double precision KA,ln1mz,ln1mlz
c$$$      double precision kQF2
c$$$**
c$$$*     Output Variables
c$$$*
c$$$      double precision clns1ccc
c$$$*
c$$$      kQF2 = 1d0                ! Q2 / muF2
c$$$      lambda = xi / ( 1d0 + xi )
c$$$      KA = ( 1d0 - lambda ) * dlog( 1d0 - lambda ) / lambda
c$$$      ln1mz  = dlog( 1d0 - z )
c$$$      ln1mlz = dlog( 1d0 - lambda * z )
c$$$*
c$$$      clns1ccc = 2d0 * CF * ( ( - 4d0 - 1d0 / 2d0 / lambda - 2d0 * zeta2
c$$$     1         - ( 1d0 + lambda ) * KA / 2d0 / lambda 
c$$$     2         + 3d0 * dlog( kQF2 / lambda ) / 2d0
c$$$     3         + 2d0 * ln1mz**2d0 - 2d0 * Rf 
c$$$     4         + 2d0 * ( - 1d0 + dlog( kQF2 / lambda ) ) * ln1mz
c$$$     5         + ( - 2d0 * lambda + 1d0 / 2d0 ) 
c$$$     6         * ln1mlz / lambda**2d0 + ( 1d0 - lambda ) * z 
c$$$     7         / 2d0 / lambda / ( 1d0 - lambda * z ) )
c$$$     8         * ( 1d0 - lambda ) + lambda * KA )
c$$$*
c$$$      return
c$$$      end
c$$$*
c$$$************************************************************************
c$$$      function clg1cca(xi,z)
c$$$*
c$$$      implicit none
c$$$*
c$$$      include "../commons/ColorFactors.h"
c$$$**
c$$$*     Input Variables
c$$$*
c$$$      double precision xi,z
c$$$**
c$$$*     Internal Variables
c$$$*
c$$$      double precision lambda
c$$$      double precision kQF2
c$$$**
c$$$*     Output Variables
c$$$*
c$$$      double precision clg1cca
c$$$*
c$$$      kQF2 = 1d0                ! Q2 / muF2
c$$$      lambda = xi / ( 1d0 + xi )
c$$$*
c$$$      clg1cca = 4d0 * TR * ( ( 1d0 - lambda ) 
c$$$     1        * ( z**2d0 + ( 1d0 - z )**2d0 ) 
c$$$     1        * ( dlog( ( 1d0 -  z ) / z ) - dlog( 1d0 - lambda ) / 2d0 
c$$$     2        + dlog(kQF2 / lambda ) / 2d0 ) 
c$$$     1        + 4d0 * ( 2d0 - lambda ) * z * ( 1d0 - z )
c$$$     3        + ( 1d0 - lambda ) * ( - 2d0 * ( 3d0 + 4d0 * lambda ) 
c$$$     4        * z * ( 1d0 - z ) 
c$$$     5        + 4d0 * lambda * z * ( 1d0 - 2d0 * lambda * z ) 
c$$$     6        * dlog( ( 1d0 - lambda * z ) / ( 1d0 - lambda ) / z ) ) )
c$$$*
c$$$      return
c$$$      end
c$$$*
c$$$************************************************************************
c$$$      function c3ns1cca(xi,z)
c$$$*
c$$$      implicit none
c$$$*
c$$$      include "../commons/ColorFactors.h"
c$$$**
c$$$*     Input Variables
c$$$*
c$$$      double precision xi,z
c$$$**
c$$$*     Internal Variables
c$$$*
c$$$      double precision lambda
c$$$      double precision kQF2
c$$$**
c$$$*     Output Variables
c$$$*
c$$$      double precision c3ns1cca
c$$$*
c$$$      kQF2 = 1d0                ! Q2 / muF2
c$$$      lambda = xi / ( 1d0 + xi )
c$$$*
c$$$      c3ns1cca = 2d0 * CF * ( - ( 1d0 + z**2d0 ) * dlog(z) / ( 1d0 - z ) 
c$$$     1         + ( 1d0 - dlog( kQF2 / lambda ) 
c$$$     2         -  2d0 * dlog( 1d0 - z ) + dlog( 1d0 - lambda * z ) ) 
c$$$     3         * ( 1d0 + z ) + 1d0 / lambda )
c$$$*
c$$$      return
c$$$      end
c$$$*
c$$$************************************************************************
c$$$      function c3ns1ccb(xi,z)
c$$$*
c$$$      implicit none
c$$$*
c$$$      include "../commons/ColorFactors.h"
c$$$**
c$$$*     Input Variables
c$$$*
c$$$      double precision xi,z
c$$$**
c$$$*     Internal Variables
c$$$*
c$$$      double precision lambda
c$$$      double precision kQF2
c$$$**
c$$$*     Output Variables
c$$$*
c$$$      double precision c3ns1ccb
c$$$*
c$$$      kQF2 = 1d0                ! Q2 / muF2
c$$$      lambda = xi / ( 1d0 + xi )
c$$$*
c$$$      c3ns1ccb = 2d0 * CF * ( 2d0 * ( 2d0 * dlog( 1d0 - z ) 
c$$$     1         - dlog( 1d0 - lambda * z ) ) / ( 1d0 - z )
c$$$     2         + 2d0 * ( - 1d0 + dlog( kQF2 / lambda ) ) / ( 1d0 - z ) 
c$$$     3         + ( lambda - 1d0 ) / lambda 
c$$$     4         / ( 1d0 - lambda * z ) + ( 1d0 - z ) 
c$$$     5         / ( 1d0 - lambda * z )**2d0 / 2d0 )
c$$$*
c$$$      return
c$$$      end
c$$$*
c$$$************************************************************************
c$$$      function c3ns1ccc(Rf,xi,z)
c$$$*
c$$$      implicit none
c$$$*
c$$$      include "../commons/ColorFactors.h"
c$$$      include "../commons/consts.h"
c$$$**
c$$$*     Input Variables
c$$$*
c$$$      double precision Rf,xi,z
c$$$**
c$$$*     Internal Variables
c$$$*
c$$$      double precision lambda
c$$$      double precision KA,ln1mz,ln1mlz
c$$$      double precision kQF2
c$$$**
c$$$*     Output Variables
c$$$*
c$$$      double precision c3ns1ccc
c$$$*
c$$$      kQF2 = 1d0                ! Q2 / muF2
c$$$      lambda = xi / ( 1d0 + xi )
c$$$      KA = ( 1d0 - lambda ) * dlog( 1d0 - lambda ) / lambda
c$$$      ln1mz  = dlog( 1d0 - z )
c$$$      ln1mlz = dlog( 1d0 - lambda * z )
c$$$*
c$$$      c3ns1ccc = 2d0 * CF * ( - 4d0 - 1d0 / 2d0 / lambda - 2d0 * zeta2
c$$$     1         - ( 1d0 + 3d0 * lambda ) * KA / 2d0 / lambda 
c$$$     2         + 3d0 * dlog( kQF2 / lambda ) / 2d0 
c$$$     3         + 2d0 * ln1mz**2d0 - 2d0 * Rf 
c$$$     4         + 2d0 * ( - 1d0 + dlog( kQF2 / lambda ) ) * ln1mz
c$$$     5         + ( lambda - 1d0 / 2d0 ) 
c$$$     6         * ln1mlz / lambda**2d0 + ( 1d0 - lambda ) * z 
c$$$     7         / 2d0 / lambda / ( 1d0 - lambda * z ) )
c$$$*
c$$$      return
c$$$      end
c$$$*
c$$$************************************************************************
c$$$      function c3g1cca(xi,z)
c$$$*
c$$$      implicit none
c$$$*
c$$$      include "../commons/ColorFactors.h"
c$$$**
c$$$*     Input Variables
c$$$*
c$$$      double precision xi,z
c$$$**
c$$$*     Internal Variables
c$$$*
c$$$      double precision lambda
c$$$      double precision kQF2
c$$$**
c$$$*     Output Variables
c$$$*
c$$$      double precision c3g1cca
c$$$*
c$$$      kQF2 = 1d0                ! Q2 / muF2
c$$$      lambda = xi / ( 1d0 + xi )
c$$$*
c$$$      c3g1cca = 4d0 * TR * ( ( z**2d0 + ( 1d0 - z )**2d0 ) 
c$$$     1        * ( dlog( ( 1d0 -  z ) / ( 1d0 - lambda * z ) )
c$$$     2        + dlog( 1d0 - lambda ) / 2d0 
c$$$     2        + dlog(kQF2 / lambda ) / 2d0 )
c$$$     3        + ( 1d0 - lambda ) * ( 2d0 * z * ( 1d0 - z )
c$$$     5        - 2d0 * z * ( 1d0 - ( 1d0 + lambda ) * z ) 
c$$$     6        * dlog( ( 1d0 - lambda * z ) / ( 1d0 - lambda ) / z ) ) )
c$$$*
c$$$      return
c$$$      end
*
************************************************************************
      function Rfun(xi,x)
*
      implicit none
**
*     Input Variables
*
      double precision xi,x
**
*     Internal Variables
*
      double precision dgauss
      double precision FunLam
      double precision eps
      parameter(eps=1d-5)
      external FunLam

      double precision clam
      common / IntLam / clam
**
*     Output Variables
*
      double precision Rfun
*
      clam = xi / ( 1d0 + xi )
      Rfun = - dgauss(FunLam,0d0,x,eps)
*
      return
      end
*
************************************************************************
      function FunLam(x)
*
      implicit none
**
*     Input Variables
*
      double precision x
**
*     Internal Variables
*
      double precision clam
      common / IntLam / clam
**
*     Output Variables
*
      double precision FunLam
*
      FunLam = dlog( 1d0 - clam * x ) / ( 1d0 - x )
*
      return
      end
*
************************************************************************
*
*     Set of Coefficient functions to reproduce the FKgenerator
*     results.
*
************************************************************************
      function c2ns1cca(xi,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
**
*     Input Variables
*
      double precision xi,z
**
*     Internal Variables
*
      double precision lambda
      double precision kQF2
**
*     Output Variables
*
      double precision c2ns1cca
*
      kQF2 = 1d0                ! Q2 / muF2
      lambda = xi / ( 1d0 + xi )
*
      c2ns1cca = 2d0 * CF * ( - ( 1d0 + z**2d0 ) * dlog(z) / ( 1d0 - z )
     1         + ( 2d0 - dlog( kQF2 / lambda ) 
     2         -  2d0 * dlog( 1d0 - z ) + dlog( 1d0 - lambda * z ) ) 
     3         * ( 1d0 + z ) + 1d0 / lambda 
     &         + ( 2d0 * lambda**2d0 - lambda - 1d0 ) / lambda 
     &         / ( 1d0 - lambda * z ) )
*
      return
      end
*
************************************************************************
      function c2ns1ccb(xi,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
**
*     Input Variables
*
      double precision xi,z
**
*     Internal Variables
*
      double precision lambda
      double precision kQF2
**
*     Output Variables
*
      double precision c2ns1ccb
*
      kQF2 = 1d0                ! Q2 / muF2
      lambda = xi / ( 1d0 + xi )
*
      c2ns1ccb = 2d0 * CF * ( 2d0 * ( 2d0 * dlog( 1d0 - z ) 
     1         - dlog( 1d0 - lambda * z ) ) / ( 1d0 - z ) 
     2         + 2d0 * ( - 1d0 + dlog( kQF2 / lambda ) ) / ( 1d0 - z )
     3         + ( 1d0 - z ) / ( 1d0 - lambda * z )**2d0 / 2d0 )
*
      return
      end
*
************************************************************************
      function c2ns1ccc(Rf,xi,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
      include "../commons/consts.h"
**
*     Input Variables
*
      double precision Rf,xi,z
**
*     Internal Variables
*
      double precision lambda
      double precision KA,ln1mz,ln1mlz
      double precision kQF2
**
*     Output Variables
*
      double precision c2ns1ccc
*
      kQF2 = 1d0                ! Q2 / muF2
      lambda = xi / ( 1d0 + xi )
      KA = ( 1d0 - lambda ) * dlog( 1d0 - lambda ) / lambda
      ln1mz  = dlog( 1d0 - z )
      ln1mlz = dlog( 1d0 - lambda * z )
*
      c2ns1ccc = 2d0 * CF * ( - 4d0 - 1d0 / 2d0 / lambda - 2d0 * zeta2 
     1         - ( 1d0 + lambda ) * KA / 2d0 / lambda 
     2         + 3d0 * dlog( kQF2 / lambda ) / 2d0
     3         + 2d0 * ln1mz**2d0 - 2d0 * Rf 
     4         + 2d0 * ( - 1d0 + dlog( kQF2 / lambda ) ) * ln1mz
     5         + ( 1d0 / 2d0 ) 
     6         * ln1mlz / lambda**2d0 + ( 1d0 - lambda ) * z 
     7         / 2d0 / lambda / ( 1d0 - lambda * z ) )
*
      return
      end
*
************************************************************************
      function c2g1cca(xi,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
**
*     Input Variables
*
      double precision xi,z
**
*     Internal Variables
*
      double precision lambda
      double precision kQF2
**
*     Output Variables
*
      double precision c2g1cca
*
      kQF2 = 1d0                ! Q2 / muF2
      lambda = xi / ( 1d0 + xi )
*
      c2g1cca = 4d0 * TR * ( ( z**2d0 + ( 1d0 - z )**2d0 )
     1        * ( dlog( ( 1d0 -  z ) / z ) - dlog( 1d0 - lambda ) / 2d0 
     2        + dlog( kQF2 / lambda ) / 2d0 ) + 8d0 * z * ( 1d0 - z )
     3        - 1d0 + ( 1d0 - lambda ) * ( - 6d0 * ( 1d0 + 2d0*lambda )
     4        * z * ( 1d0 - z ) + 1d0 / ( 1d0 - lambda * z ) 
     5        + 6d0 * lambda * z * ( 1d0 - 2d0 * lambda * z ) 
     6        * dlog( ( 1d0 - lambda * z ) / ( 1d0 - lambda ) / z ) ) )
*
      return
      end
*
************************************************************************
      function clns1cca(xi,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
**
*     Input Variables
*
      double precision xi,z
**
*     Internal Variables
*
      double precision lambda
      double precision kQF2
**
*     Output Variables
*
      double precision clns1cca
*
      kQF2 = 1d0                ! Q2 / muF2
      lambda = xi / ( 1d0 + xi )
*
      clns1cca = 2d0 * CF * ( ( - ( 1d0 + z**2d0 ) * dlog(z) 
     1         / ( 1d0 - z ) + ( - dlog( kQF2 / lambda ) 
     2         -  2d0 * dlog( 1d0 - z ) + dlog( 1d0 - lambda * z ) ) 
     3         * ( 1d0 + z ) + 3d0 - 2d0 / ( 1d0 - lambda * z ) )
     4         * ( 1d0 - lambda ) + ( 1d0 + lambda ) * z )
*
      return
      end
*
************************************************************************
      function clns1ccb(xi,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
**
*     Input Variables
*
      double precision xi,z
**
*     Internal Variables
*
      double precision lambda
      double precision kQF2
**
*     Output Variables
*
      double precision clns1ccb
*
      kQF2 = 1d0                ! Q2 / muF2
      lambda = xi / ( 1d0 + xi )
*
      clns1ccb = 2d0 * CF * ( 2d0 * ( 2d0 * dlog( 1d0 - z ) 
     1         - dlog( 1d0 - lambda * z ) ) / ( 1d0 - z )
     2         + 2d0 * ( - 1d0 + dlog( kQF2 / lambda ) ) / ( 1d0 - z ) 
     3         + ( 1d0 - z ) 
     4         / ( 1d0 - lambda * z )**2d0 / 2d0 ) * ( 1d0 - lambda )
*
      return
      end
*
************************************************************************
      function clns1ccc(Rf,xi,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
      include "../commons/consts.h"
**
*     Input Variables
*
      double precision Rf,xi,z
**
*     Internal Variables
*
      double precision lambda
      double precision KA,ln1mz,ln1mlz
      double precision kQF2
**
*     Output Variables
*
      double precision clns1ccc
*
      kQF2 = 1d0                ! Q2 / muF2
      lambda = xi / ( 1d0 + xi )
      KA = ( 1d0 - lambda ) * dlog( 1d0 - lambda ) / lambda
      ln1mz  = dlog( 1d0 - z )
      ln1mlz = dlog( 1d0 - lambda * z )
*
      clns1ccc = 2d0 * CF * ( ( - 4d0 - 1d0 / 2d0 / lambda - 2d0 * zeta2
     1         - ( 1d0 + lambda ) * KA / 2d0 / lambda 
     2         + 3d0 * dlog( kQF2 / lambda ) / 2d0
     3         + 2d0 * ln1mz**2d0 - 2d0 * Rf 
     4         + 2d0 * ( - 1d0 + dlog( kQF2 / lambda ) ) * ln1mz
     5         + ( + 1d0 / 2d0 ) 
     6         * ln1mlz / lambda**2d0 + ( 1d0 - lambda ) * z 
     7         / 2d0 / lambda / ( 1d0 - lambda * z ) )
     8         * ( 1d0 - lambda ) + lambda * KA )
*
      return
      end
*
************************************************************************
      function clg1cca(xi,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
**
*     Input Variables
*
      double precision xi,z
**
*     Internal Variables
*
      double precision lambda
      double precision kQF2
**
*     Output Variables
*
      double precision clg1cca
*
      kQF2 = 1d0                ! Q2 / muF2
      lambda = xi / ( 1d0 + xi )
*
      clg1cca = 4d0 * TR * ( ( 1d0 - lambda ) 
     1        * ( z**2d0 + ( 1d0 - z )**2d0 ) 
     1        * ( dlog( ( 1d0 -  z ) / z ) - dlog( 1d0 - lambda ) / 2d0 
     2        + dlog(kQF2 / lambda ) / 2d0 ) 
     1        + 4d0 * ( 2d0 - lambda ) * z * ( 1d0 - z )
     3        + ( 1d0 - lambda ) * ( - 2d0 * ( 3d0 + 4d0 * lambda ) 
     4        * z * ( 1d0 - z ) 
     5        + 4d0 * lambda * z * ( 1d0 - 2d0 * lambda * z ) 
     6        * dlog( ( 1d0 - lambda * z ) / ( 1d0 - lambda ) / z ) ) )
*
      return
      end
*
************************************************************************
      function c3ns1cca(xi,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
**
*     Input Variables
*
      double precision xi,z
**
*     Internal Variables
*
      double precision lambda
      double precision kQF2
**
*     Output Variables
*
      double precision c3ns1cca
*
      kQF2 = 1d0                ! Q2 / muF2
      lambda = xi / ( 1d0 + xi )
*
      c3ns1cca = 2d0 * CF * ( - ( 1d0 + z**2d0 ) * dlog(z) / ( 1d0 - z ) 
     1         + ( 1d0 - dlog( kQF2 / lambda ) 
     2         -  2d0 * dlog( 1d0 - z ) + dlog( 1d0 - lambda * z ) ) 
     3         * ( 1d0 + z ) + 1d0 / lambda + ( lambda - 1d0 ) / lambda 
     4         / ( 1d0 - lambda * z ) )
*
      return
      end
*
************************************************************************
      function c3ns1ccb(xi,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
**
*     Input Variables
*
      double precision xi,z
**
*     Internal Variables
*
      double precision lambda
      double precision kQF2
**
*     Output Variables
*
      double precision c3ns1ccb
*
      kQF2 = 1d0                ! Q2 / muF2
      lambda = xi / ( 1d0 + xi )
*
      c3ns1ccb = 2d0 * CF * ( 2d0 * ( 2d0 * dlog( 1d0 - z ) 
     1         - dlog( 1d0 - lambda * z ) ) / ( 1d0 - z )
     2         + 2d0 * ( - 1d0 + dlog( kQF2 / lambda ) ) / ( 1d0 - z ) 
     3         + ( 1d0 - z ) 
     5         / ( 1d0 - lambda * z )**2d0 / 2d0 )
*
      return
      end
*
************************************************************************
      function c3ns1ccc(Rf,xi,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
      include "../commons/consts.h"
**
*     Input Variables
*
      double precision Rf,xi,z
**
*     Internal Variables
*
      double precision lambda
      double precision KA,ln1mz,ln1mlz
      double precision kQF2
**
*     Output Variables
*
      double precision c3ns1ccc
*
      kQF2 = 1d0                ! Q2 / muF2
      lambda = xi / ( 1d0 + xi )
      KA = ( 1d0 - lambda ) * dlog( 1d0 - lambda ) / lambda
      ln1mz  = dlog( 1d0 - z )
      ln1mlz = dlog( 1d0 - lambda * z )
*
      c3ns1ccc = 2d0 * CF * ( - 4d0 - 1d0 / 2d0 / lambda - 2d0 * zeta2
     1         - ( 1d0 + 3d0 * lambda ) * KA / 2d0 / lambda 
     2         + 3d0 * dlog( kQF2 / lambda ) / 2d0 
     3         + 2d0 * ln1mz**2d0 - 2d0 * Rf 
     4         + 2d0 * ( - 1d0 + dlog( kQF2 / lambda ) ) * ln1mz
     5         + ( + 1d0 / 2d0 ) 
     6         * ln1mlz / lambda**2d0 + ( 1d0 - lambda ) * z 
     7         / 2d0 / lambda / ( 1d0 - lambda * z ) )
*
      return
      end
*
************************************************************************
      function c3g1cca(xi,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
**
*     Input Variables
*
      double precision xi,z
**
*     Internal Variables
*
      double precision lambda
      double precision kQF2
**
*     Output Variables
*
      double precision c3g1cca
*
      kQF2 = 1d0                ! Q2 / muF2
      lambda = xi / ( 1d0 + xi )
*
      c3g1cca = 4d0 * TR * ( ( z**2d0 + ( 1d0 - z )**2d0 ) 
     1        * ( dlog( ( 1d0 -  z ) / ( 1d0 - lambda * z ) )
     2        + dlog( 1d0 - lambda ) / 2d0 
     2        + dlog(kQF2 / lambda ) / 2d0 )
     3        + ( 1d0 - lambda ) * ( 2d0 * z * ( 1d0 - z )
     5        - 2d0 * z * ( 1d0 - ( 1d0 + lambda ) * z ) 
     6        * dlog( ( 1d0 - lambda * z ) / ( 1d0 - lambda ) / z ) ) )
*
      return
      end
