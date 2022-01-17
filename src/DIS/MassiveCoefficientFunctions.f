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
      include "../commons/MassRunning.h"
      include "../commons/krenQ.h"
c      include "../commons/Nf_FF.h"
**
*     input variables
*
      integer icoef,ixi
      double precision z
**
*     internal variables
*
      double precision xi,rho,eta
      double precision cm21g,cml1g
      double precision cm22q,cml2q
      double precision dcm21g,dcml1g
      double precision h1
c      double precision APFc2log,APFcllog
      double precision APFc2nlog,APFclnlog,APFc2nlobarg,APFclnlobarg
      double precision APFc2nloq,APFclnloq,APFc2nlobarq,APFclnlobarq
c      double precision APFd2nloq,APFdlnloq
c      double precision beta0apf
**
*     output variables
*
      double precision MassiveCF
*
      xi  = xigrid(ixi)
      rho = 1d0 / ( 4d0 / xi + 1d0 )
      eta = xi * ( 1d0 / z - 1d0 ) / 4d0 - 1d0
*
      MassiveCF = 0d0
      if(z.ge.rho) return
*
*     Exact expressions
*
      if(icoef.eq.1)then
         MassiveCF = 2d0 * cm21g(xi,z)
c         MassiveCF = xi * APFcllog(eta,xi) / pi / z ! From hqcoef.f
      elseif(icoef.eq.2)then
         MassiveCF = 2d0 * cml1g(xi,z)
c         MassiveCF = xi * APFc2log(eta,xi) / pi / z ! From hqcoef.f
      elseif(icoef.eq.7)then
         MassiveCF = cm22q(xi,z)
c         MassiveCF = 16d0 * pi * xi * APFd2nloq(eta,xi) / z ! From hqcoef.f
      elseif(icoef.eq.8)then
         MassiveCF = cml2q(xi,z)
c         MassiveCF = 16d0 * pi * xi * APFdlnloq(eta,xi) / z ! From hqcoef.f
      elseif(icoef.eq.3)then
         MassiveCF = 16d0 * pi * xi * APFc2nlog(eta,xi) / z
      elseif(icoef.eq.4)then
         MassiveCF = 16d0 * pi * xi * APFc2nloq(eta,xi) / z
      elseif(icoef.eq.5)then
         MassiveCF = 16d0 * pi * xi * APFclnlog(eta,xi) / z
      elseif(icoef.eq.6)then
         MassiveCF = 16d0 * pi * xi * APFclnloq(eta,xi) / z
      elseif(icoef.eq.9)then
         MassiveCF = 16d0 * pi * xi * APFc2nlobarg(eta,xi) / z
      elseif(icoef.eq.10)then
         MassiveCF = 16d0 * pi * xi * APFc2nlobarq(eta,xi) / z
      elseif(icoef.eq.11)then
         MassiveCF = 16d0 * pi * xi * APFclnlobarg(eta,xi) / z
      elseif(icoef.eq.12)then
         MassiveCF = 16d0 * pi * xi * APFclnlobarq(eta,xi) / z
      endif
*
*     If the MSbar masses are to be used, add the appropriate term to
*     the NNLO gluon coefficient functions.
*
      if(mass_scheme.eq."MSbar")then
         if(icoef.eq.3)then
            h1 = CF * 4d0
            if(MassRunning) h1 = h1 + CF * 3d0 * dlog(xi*krenQ)
            MassiveCF = MassiveCF + 2d0 * h1 * dcm21g(xi,z)
c            if(.not.MassRunning)then
c               MassiveCF = MassiveCF + 2d0 * beta0apf(Nf_FF)
c     1                   * cm21g(xi,z) * dlog(xi)
c            endif
         elseif(icoef.eq.5)then
            h1 = CF * 4d0
            if(MassRunning) h1 = h1 + CF * 3d0 * dlog(xi*krenQ)
            MassiveCF = MassiveCF + 2d0 * h1 * dcml1g(xi,z)
c            if(.not.MassRunning)then
c               MassiveCF = MassiveCF + 2d0 * beta0apf(Nf_FF)
c     1                   * cml1g(xi,z) * dlog(xi)
c            endif
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
      cm21g = ( z**2 + ( 1d0 - z )**2 
     1      + 4d0 * eps * z * ( 1d0 - 3d0 * z ) 
     2      - 8d0 * eps**2 * z**2 ) 
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
*     "cm21g": dcm21g = M * d(cm21g)/dM = - 2 * xi * d(cm21g)/d(xi).
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
      dcm21g = 4d0 * eps * ( ( - 6d0 - 8d0 * eps ) * z**2 + 2d0 * z )
     1       * dlog( ( 1d0 + v ) / ( 1d0 - v ) ) 
     2       + 8d0 * eps * z * ( z - 1d0 ) * v 
     3       - 2d0 * ( ( - 8d0 * eps**2 - 12d0 * eps + 2d0 ) * z**2
     4       - ( 2d0 - 4d0 * eps ) * z + 1d0 ) / v 
     5       + 4d0 * eps * z * ( - 4d0 * ( 2d0 - eps ) * z**2
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
      cml1g = - 8d0 * eps * z**2 
     1      * dlog( ( 1d0 + v ) / ( 1d0 - v ) )
     2      + 4d0 * v * z * ( 1d0 - z )
*
      return
      end
*
************************************************************************
*
*     Logarithmic derivative with rescpect to the mass M of the function
*     "cml1g": dcml1g = M * d(cml1g)/dM = - 2 * xi * d(cml1g)/d(xi).
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
      dcml1g = - 16d0 * eps * z**2 * dlog( ( 1d0 + v ) / ( 1d0 - v ) ) 
     1       + 16d0 * eps * z**2 / v 
     2       + 4d0 * eps * z * ( - 4d0 * z**2 + 4d0 * z ) 
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
      if(sq1.eq.sq2) return
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
      cm22q=2d0/3d0*((4d0/3d0*(1d0+z**2)/(1d0-z)
     1     -16d0/(1d0-z)*(z/xi)**2*(1d0-9d0*z+9d0*z**2))
     2     *(dlog((1d0-z)/z**2)*l1+l1*l2+2d0*(-dil1+dil2+dil3-dil4))
     3     +(-8d0/3d0+4d0/(1d0-z)+(z/(1d0-z)/xi)**2*(128d0-432d0*z
     4     +288d0*z**2-8d0/(1d0-z)))*l1+(88d0/9d0+136d0/9d0*z
     5     -152d0/9d0/(1d0-z)+(z/(1d0-z)/xi)*(464d0/9d0-512d0/3d0*z
     6     +2048d0/9d0*z**2)+(z/(1d0-z)/xi)**2*(-832d0/9d0
     7     +6208d0/9d0*z-11392d0/9d0*z**2+6016d0/9d0*z**3))*l3/sq2
     8     +(-272d0/27d0-1244d0/27d0*z+718d0/27d0/(1d0-z)+(z/(1d0-z)/xi)
     9     *(-3424d0/27d0+15608d0/27d0*z-4304d0/9d0*z**2
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
      if(sq1.eq.sq2) return
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
      cml2q=2d0/3d0*(96d0*z*(z/xi)**2*(dlog((1d0-z)/z**2)*l1
     1     +l1*l2+2d0*(-dil1+dil2+dil3-dil4))+(z/xi/(1d0-z))**2
     2     *(64d0-288d0*z+192d0*z**2)*l1+z*(16d0/3d0-416d0*z/3d0/xi
     3     +1408d0*z**2/3d0/xi**2)*l3/sq2+(16d0/3d0-400d0*z/18d0
     4     +z*(-160d0/3d0+3824d0*z/9d0-992d0*z**2/3d0)/(1d0-z)/xi)
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
      function c2ns1cca(xi,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
      include "../commons/kfacQ.h"
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
      kQF2 = 1d0 / kfacQ                ! Q2 / muF2
      lambda = xi / ( 1d0 + xi )
*
      c2ns1cca = 2d0 * CF * ( - ( 1d0 + z**2 ) * dlog(z) / ( 1d0 - z )
     1         + ( 2d0 - dlog( kQF2 / lambda ) 
     2         -  2d0 * dlog( 1d0 - z ) + dlog( 1d0 - lambda * z ) ) 
     3         * ( 1d0 + z ) + 1d0 / lambda
     4         + ( 2d0 * lambda**2 - lambda - 1d0 ) / lambda 
     5         / ( 1d0 - lambda * z ) ) 
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
      include "../commons/kfacQ.h"
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
      kQF2 = 1d0 / kfacQ                ! Q2 / muF2
      lambda = xi / ( 1d0 + xi )
*
      c2ns1ccb = 2d0 * CF * ( 2d0 * ( 2d0 * dlog( 1d0 - z ) 
     1         - dlog( 1d0 - lambda * z ) ) / ( 1d0 - z ) 
     2         + 2d0 * ( - 1d0 + dlog( kQF2 / lambda ) ) / ( 1d0 - z )
     3         + ( 1d0 - z ) / ( 1d0 - lambda * z )**2 / 2d0 )
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
      include "../commons/kfacQ.h"
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
      kQF2 = 1d0 / kfacQ                ! Q2 / muF2
      lambda = xi / ( 1d0 + xi )
      KA = ( 1d0 - lambda ) * dlog( 1d0 - lambda ) / lambda
      ln1mz  = dlog( 1d0 - z )
      ln1mlz = dlog( 1d0 - lambda * z )
*
      c2ns1ccc = 2d0 * CF * ( - 4d0 - 1d0 / 2d0 / lambda - 2d0 * zeta2 
     1         - ( 1d0 + lambda ) * KA / 2d0 / lambda 
     2         + 3d0 * dlog( kQF2 / lambda ) / 2d0
     3         + 2d0 * ln1mz**2 - 2d0 * Rf 
     4         + 2d0 * ( - 1d0 + dlog( kQF2 / lambda ) ) * ln1mz
     5         + ln1mlz / 2d0 / lambda**2 + ( 1d0 - lambda ) * z 
     6         / 2d0 / lambda / ( 1d0 - lambda * z ) )
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
      include "../commons/kfacQ.h"
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
      kQF2 = 1d0 / kfacQ                ! Q2 / muF2
      lambda = xi / ( 1d0 + xi )
*
      c2g1cca = 4d0 * TR * ( ( z**2 + ( 1d0 - z )**2 )
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
      include "../commons/kfacQ.h"
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
      kQF2 = 1d0 / kfacQ                ! Q2 / muF2
      lambda = xi / ( 1d0 + xi )
*
      clns1cca = 2d0 * CF * ( ( - ( 1d0 + z**2 ) * dlog(z) 
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
      include "../commons/kfacQ.h"
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
      kQF2 = 1d0 / kfacQ                ! Q2 / muF2
      lambda = xi / ( 1d0 + xi )
*
      clns1ccb = 2d0 * CF * ( 2d0 * ( 2d0 * dlog( 1d0 - z ) 
     1         - dlog( 1d0 - lambda * z ) ) / ( 1d0 - z )
     2         + 2d0 * ( - 1d0 + dlog( kQF2 / lambda ) ) / ( 1d0 - z ) 
     3         + ( 1d0 - z ) / ( 1d0 - lambda * z )**2 / 2d0 )
     4         * ( 1d0 - lambda )
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
      include "../commons/kfacQ.h"
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
      kQF2 = 1d0 / kfacQ                ! Q2 / muF2
      lambda = xi / ( 1d0 + xi )
      KA = ( 1d0 - lambda ) * dlog( 1d0 - lambda ) / lambda
      ln1mz  = dlog( 1d0 - z )
      ln1mlz = dlog( 1d0 - lambda * z )
*
      clns1ccc = 2d0 * CF * ( ( - 4d0 - 1d0 / 2d0 / lambda - 2d0 * zeta2
     1         - ( 1d0 + lambda ) * KA / 2d0 / lambda 
     2         + 3d0 * dlog( kQF2 / lambda ) / 2d0
     3         + 2d0 * ln1mz**2 - 2d0 * Rf 
     4         + 2d0 * ( - 1d0 + dlog( kQF2 / lambda ) ) * ln1mz
     5         + ln1mlz / 2d0 / lambda**2 + ( 1d0 - lambda ) * z 
     6         / 2d0 / lambda / ( 1d0 - lambda * z ) )
     7         * ( 1d0 - lambda ) + lambda * KA )
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
      include "../commons/kfacQ.h"
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
      kQF2 = 1d0 / kfacQ                ! Q2 / muF2
      lambda = xi / ( 1d0 + xi )
*
      clg1cca = 4d0 * TR * ( ( 1d0 - lambda ) 
     1        * ( z**2 + ( 1d0 - z )**2 ) 
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
      include "../commons/kfacQ.h"
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
      kQF2 = 1d0 / kfacQ                ! Q2 / muF2
      lambda = xi / ( 1d0 + xi )
*
      c3ns1cca = 2d0 * CF * ( - ( 1d0 + z**2 ) * dlog(z) / ( 1d0 - z ) 
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
      include "../commons/kfacQ.h"
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
      kQF2 = 1d0 / kfacQ                ! Q2 / muF2
      lambda = xi / ( 1d0 + xi )
*
      c3ns1ccb = 2d0 * CF * ( 2d0 * ( 2d0 * dlog( 1d0 - z ) 
     1         - dlog( 1d0 - lambda * z ) ) / ( 1d0 - z )
     2         + 2d0 * ( - 1d0 + dlog( kQF2 / lambda ) ) / ( 1d0 - z ) 
     3         + ( 1d0 - z ) / ( 1d0 - lambda * z )**2 / 2d0 )
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
      include "../commons/kfacQ.h"
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
      kQF2 = 1d0 / kfacQ                ! Q2 / muF2
      lambda = xi / ( 1d0 + xi )
      KA = ( 1d0 - lambda ) * dlog( 1d0 - lambda ) / lambda
      ln1mz  = dlog( 1d0 - z )
      ln1mlz = dlog( 1d0 - lambda * z )
*
      c3ns1ccc = 2d0 * CF * ( - 4d0 - 1d0 / 2d0 / lambda - 2d0 * zeta2
     1         - ( 1d0 + 3d0 * lambda ) * KA / 2d0 / lambda 
     2         + 3d0 * dlog( kQF2 / lambda ) / 2d0 
     3         + 2d0 * ln1mz**2 - 2d0 * Rf 
     4         + 2d0 * ( - 1d0 + dlog( kQF2 / lambda ) ) * ln1mz
     5         + ln1mlz / 2d0 / lambda**2 + ( 1d0 - lambda ) * z 
     6         / 2d0 / lambda / ( 1d0 - lambda * z ) )
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
      include "../commons/kfacQ.h"
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
      kQF2 = 1d0 / kfacQ                ! Q2 / muF2
      lambda = xi / ( 1d0 + xi )
*
      c3g1cca = 4d0 * TR * ( ( z**2 + ( 1d0 - z )**2 ) 
     1        * ( dlog( ( 1d0 -  z ) / ( 1d0 - lambda * z ) )
     2        + dlog( 1d0 - lambda ) / 2d0 
     2        + dlog(kQF2 / lambda ) / 2d0 )
     3        + ( 1d0 - lambda ) * ( 2d0 * z * ( 1d0 - z )
     5        - 2d0 * z * ( 1d0 - ( 1d0 + lambda ) * z ) 
     6        * dlog( ( 1d0 - lambda * z ) / ( 1d0 - lambda ) / z ) ) )
*
      return
      end
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
c$$$************************************************************************
c$$$*
c$$$*     Set of Coefficient functions to reproduce the FKgenerator
c$$$*     results.
c$$$*
c$$$************************************************************************
c$$$      function c2ns1cca(xi,z)
c$$$*
c$$$      implicit none
c$$$*
c$$$      include "../commons/ColorFactors.h"
c$$$      include "../commons/kfacQ.h"
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
c$$$      kQF2 = 1d0 / kfacQ                ! Q2 / muF2
c$$$      lambda = xi / ( 1d0 + xi )
c$$$*
c$$$      c2ns1cca = 2d0 * CF * ( - ( 1d0 + z**2 ) * dlog(z) / ( 1d0 - z )
c$$$     1         + ( 2d0 - dlog( kQF2 / lambda ) 
c$$$     2         -  2d0 * dlog( 1d0 - z ) + dlog( 1d0 - lambda * z ) ) 
c$$$     3         * ( 1d0 + z ) + 1d0 / lambda 
c$$$     &         + ( 2d0 * lambda**2 - lambda - 1d0 ) / lambda 
c$$$     &         / ( 1d0 - lambda * z ) )
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
c$$$      include "../commons/kfacQ.h"
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
c$$$      kQF2 = 1d0 / kfacQ                ! Q2 / muF2
c$$$      lambda = xi / ( 1d0 + xi )
c$$$*
c$$$      c2ns1ccb = 2d0 * CF * ( 2d0 * ( 2d0 * dlog( 1d0 - z ) 
c$$$     1         - dlog( 1d0 - lambda * z ) ) / ( 1d0 - z ) 
c$$$     2         + 2d0 * ( - 1d0 + dlog( kQF2 / lambda ) ) / ( 1d0 - z )
c$$$     3         + ( 1d0 - z ) / ( 1d0 - lambda * z )**2 / 2d0 )
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
c$$$      include "../commons/kfacQ.h"
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
c$$$      kQF2 = 1d0 / kfacQ                ! Q2 / muF2
c$$$      lambda = xi / ( 1d0 + xi )
c$$$      KA = ( 1d0 - lambda ) * dlog( 1d0 - lambda ) / lambda
c$$$      ln1mz  = dlog( 1d0 - z )
c$$$      ln1mlz = dlog( 1d0 - lambda * z )
c$$$*
c$$$      c2ns1ccc = 2d0 * CF * ( - 4d0 - 1d0 / 2d0 / lambda - 2d0 * zeta2 
c$$$     1         - ( 1d0 + lambda ) * KA / 2d0 / lambda 
c$$$     2         + 3d0 * dlog( kQF2 / lambda ) / 2d0
c$$$     3         + 2d0 * ln1mz**2 - 2d0 * Rf 
c$$$     4         + 2d0 * ( - 1d0 + dlog( kQF2 / lambda ) ) * ln1mz
c$$$     5         + ( 1d0 / 2d0 ) 
c$$$     6         * ln1mlz / lambda**2 + ( 1d0 - lambda ) * z 
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
c$$$      include "../commons/kfacQ.h"
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
c$$$      kQF2 = 1d0 / kfacQ                ! Q2 / muF2
c$$$      lambda = xi / ( 1d0 + xi )
c$$$*
c$$$      c2g1cca = 4d0 * TR * ( ( z**2 + ( 1d0 - z )**2 )
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
c$$$      include "../commons/kfacQ.h"
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
c$$$      kQF2 = 1d0 / kfacQ                ! Q2 / muF2
c$$$      lambda = xi / ( 1d0 + xi )
c$$$*
c$$$      clns1cca = 2d0 * CF * ( ( - ( 1d0 + z**2 ) * dlog(z) 
c$$$     1         / ( 1d0 - z ) + ( - dlog( kQF2 / lambda ) 
c$$$     2         -  2d0 * dlog( 1d0 - z ) + dlog( 1d0 - lambda * z ) ) 
c$$$     3         * ( 1d0 + z ) + 3d0 - 2d0 / ( 1d0 - lambda * z ) )
c$$$     4         * ( 1d0 - lambda ) + ( 1d0 + lambda ) * z )
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
c$$$      include "../commons/kfacQ.h"
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
c$$$      kQF2 = 1d0 / kfacQ                ! Q2 / muF2
c$$$      lambda = xi / ( 1d0 + xi )
c$$$*
c$$$      clns1ccb = 2d0 * CF * ( 2d0 * ( 2d0 * dlog( 1d0 - z ) 
c$$$     1         - dlog( 1d0 - lambda * z ) ) / ( 1d0 - z )
c$$$     2         + 2d0 * ( - 1d0 + dlog( kQF2 / lambda ) ) / ( 1d0 - z ) 
c$$$     3         + ( 1d0 - z ) 
c$$$     4         / ( 1d0 - lambda * z )**2 / 2d0 ) * ( 1d0 - lambda )
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
c$$$      include "../commons/kfacQ.h"
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
c$$$      kQF2 = 1d0 / kfacQ                ! Q2 / muF2
c$$$      lambda = xi / ( 1d0 + xi )
c$$$      KA = ( 1d0 - lambda ) * dlog( 1d0 - lambda ) / lambda
c$$$      ln1mz  = dlog( 1d0 - z )
c$$$      ln1mlz = dlog( 1d0 - lambda * z )
c$$$*
c$$$      clns1ccc = 2d0 * CF * ( ( - 4d0 - 1d0 / 2d0 / lambda - 2d0 * zeta2
c$$$     1         - ( 1d0 + lambda ) * KA / 2d0 / lambda 
c$$$     2         + 3d0 * dlog( kQF2 / lambda ) / 2d0
c$$$     3         + 2d0 * ln1mz**2 - 2d0 * Rf 
c$$$     4         + 2d0 * ( - 1d0 + dlog( kQF2 / lambda ) ) * ln1mz
c$$$     5         + ( + 1d0 / 2d0 ) 
c$$$     6         * ln1mlz / lambda**2 + ( 1d0 - lambda ) * z 
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
c$$$      include "../commons/kfacQ.h"
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
c$$$      kQF2 = 1d0 / kfacQ                ! Q2 / muF2
c$$$      lambda = xi / ( 1d0 + xi )
c$$$*
c$$$      clg1cca = 4d0 * TR * ( ( 1d0 - lambda ) 
c$$$     1        * ( z**2 + ( 1d0 - z )**2 ) 
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
c$$$      include "../commons/kfacQ.h"
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
c$$$      kQF2 = 1d0 / kfacQ                ! Q2 / muF2
c$$$      lambda = xi / ( 1d0 + xi )
c$$$*
c$$$      c3ns1cca = 2d0 * CF * ( - ( 1d0 + z**2 ) * dlog(z) / ( 1d0 - z ) 
c$$$     1         + ( 1d0 - dlog( kQF2 / lambda ) 
c$$$     2         -  2d0 * dlog( 1d0 - z ) + dlog( 1d0 - lambda * z ) ) 
c$$$     3         * ( 1d0 + z ) + 1d0 / lambda + ( lambda - 1d0 ) / lambda 
c$$$     4         / ( 1d0 - lambda * z ) )
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
c$$$      include "../commons/kfacQ.h"
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
c$$$      kQF2 = 1d0 / kfacQ                ! Q2 / muF2
c$$$      lambda = xi / ( 1d0 + xi )
c$$$*
c$$$      c3ns1ccb = 2d0 * CF * ( 2d0 * ( 2d0 * dlog( 1d0 - z ) 
c$$$     1         - dlog( 1d0 - lambda * z ) ) / ( 1d0 - z )
c$$$     2         + 2d0 * ( - 1d0 + dlog( kQF2 / lambda ) ) / ( 1d0 - z ) 
c$$$     3         + ( 1d0 - z ) 
c$$$     5         / ( 1d0 - lambda * z )**2 / 2d0 )
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
c$$$      include "../commons/kfacQ.h"
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
c$$$      kQF2 = 1d0 / kfacQ                ! Q2 / muF2
c$$$      lambda = xi / ( 1d0 + xi )
c$$$      KA = ( 1d0 - lambda ) * dlog( 1d0 - lambda ) / lambda
c$$$      ln1mz  = dlog( 1d0 - z )
c$$$      ln1mlz = dlog( 1d0 - lambda * z )
c$$$*
c$$$      c3ns1ccc = 2d0 * CF * ( - 4d0 - 1d0 / 2d0 / lambda - 2d0 * zeta2
c$$$     1         - ( 1d0 + 3d0 * lambda ) * KA / 2d0 / lambda 
c$$$     2         + 3d0 * dlog( kQF2 / lambda ) / 2d0 
c$$$     3         + 2d0 * ln1mz**2 - 2d0 * Rf 
c$$$     4         + 2d0 * ( - 1d0 + dlog( kQF2 / lambda ) ) * ln1mz
c$$$     5         + ( + 1d0 / 2d0 ) 
c$$$     6         * ln1mlz / lambda**2 + ( 1d0 - lambda ) * z 
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
c$$$      include "../commons/kfacQ.h"
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
c$$$      kQF2 = 1d0 / kfacQ                ! Q2 / muF2
c$$$      lambda = xi / ( 1d0 + xi )
c$$$*
c$$$      c3g1cca = 4d0 * TR * ( ( z**2 + ( 1d0 - z )**2 ) 
c$$$     1        * ( dlog( ( 1d0 -  z ) / ( 1d0 - lambda * z ) )
c$$$     2        + dlog( 1d0 - lambda ) / 2d0 
c$$$     2        + dlog(kQF2 / lambda ) / 2d0 )
c$$$     3        + ( 1d0 - lambda ) * ( 2d0 * z * ( 1d0 - z )
c$$$     5        - 2d0 * z * ( 1d0 - ( 1d0 + lambda ) * z ) 
c$$$     6        * dlog( ( 1d0 - lambda * z ) / ( 1d0 - lambda ) / z ) ) )
c$$$*
c$$$      return
c$$$      end
c$$$*
************************************************************************
*
*     O(as) coefficient functions for massive initial state charm.
*     The following expressions hold both for NC and CC processes.
*     The z-independent kinematics is bassed by common.
*     Reference: hep-ph/9805233
*
************************************************************************
*     Regular functions that multiplies 1 / ( 1 - z )
************************************************************************
      function c11ICR(z)
*
      implicit none
*
      include "../commons/wrapIC.h"
**
*     Input Variables
*
      double precision z
**
*     Internal Variables
*
      double precision s1h,s1h2
      double precision Lxi,Ixi
      double precision DeltaFun,Delp,Delp2
      double precision f1hat,N1
**
*     Output Variables
*
      double precision c11ICR
*
      s1h   = ( 1d0 - z ) * ( ( Del - Spm ) * z + Del + Spm ) / 2d0 / z
      s1h2  = s1h * s1h
      Delp  = DeltaFun(m12,s1h+m22,-Q2IC)
      Delp2 = Delp * Delp
      Lxi   = dlog( ( Spp + s1h - Delp ) / ( Spp + s1h + Delp ) )
      Ixi   = ( ( s1h + 2d0 * m22 ) / s1h2 + ( s1h + m22 )
     1      / Delp / s1h2 * Spp * Lxi )
*
      f1hat = 8d0 / Delp2 * ( - Del2 * ( Splus * Spp
     1     - 2d0 * m1 * m2 * Sminus ) * Ixi
     2     + 2d0 * m1 * m2 * Sminus * ( 1d0 / s1h * ( Delp2
     3     + 4d0 * m22 * Spm )
     4     + 2d0 * Spm - Smp + (Spp + s1h) / 2d0
     5     + ( s1h + m22 ) / Delp / s1h * ( Delp2
     6     + 2d0 * Spm * Spp + ( m22 + Q2IC ) * s1h ) * Lxi )
     7     + Splus * ( ( - m22 * Spp ) / ( ( s1h + m22 ) * s1h )
     8     * ( Del2 + 4d0 * m22 * Spm)
     9     - 1d0 / 4d0 / ( s1h + m22 )
     1     * ( 3d0 * Spp**2 * Smp
     2     + 4d0 * m22 * (10d0 * Spp * Spm - Spm * Smp
     3     - m12 * Spp)
     4     + s1h * ( - 7d0 * Spp * Smp + 18d0 * Del2
     5     - 4d0 * m12 * ( 7d0 * Q2IC - 4 * m22
     6     + 7d0 * m12 ) )
     7     + 3d0 * s1h2 * ( Spm - 2d0 * m12 ) - s1h**3 )
     8     + ( s1h + m22 ) / 2d0 / Delp
     9     * ( - 2d0 / s1h * Spp * ( Del2 + 2d0 * Spm * Spp )
     1     + ( 4d0 * m12 * m22 - 7d0 * Spm * Spp )
     2     - 4d0 * Spm * s1h - s1h**2 ) * Lxi ) )
*
      N1 = ( Splus * Spp - 2d0 * m1 * m2 * Sminus ) / 2d0 /  Del
*
      c11ICR = ( 1d0 - z ) * s1h * f1hat / N1 / 8d0 / ( s1h + m22 )
*
      return
      end
*
************************************************************************
      function c21ICR(z)
*
      implicit none
*
      include "../commons/wrapIC.h"
**
*     Input Variables
*
      double precision z
**
*     Internal Variables
*
      double precision s1h,s1h2
      double precision Lxi,Ixi
      double precision DeltaFun,Delp,Delp2
      double precision f2hat,N2
**
*     Output Variables
*
      double precision c21ICR
*
      s1h   = ( 1d0 - z ) * ( ( Del - Spm ) * z + Del + Spm ) / 2d0 / z
      s1h2  = s1h * s1h
      Delp  = DeltaFun(m12,s1h+m22,-Q2IC)
      Delp2 = Delp * Delp
      Lxi   = dlog( ( Spp + s1h - Delp ) / ( Spp + s1h + Delp ) )
      Ixi   = ( ( s1h + 2d0 * m22 ) / s1h2 + ( s1h + m22 )
     1      / Delp / s1h2 * Spp * Lxi )
*
      f2hat = 16d0 / Delp**4 * ( - 2d0 * Del**4 * Splus * Ixi
     1     + 2d0 * m1 * m2 * Sminus * ( ( ( s1h + m22 ) / Delp )
     2     * ( Delp2 - 6d0 * m12 * Q2IC ) * Lxi
     3     - Delp2 * ( s1h + Spp ) / 2d0 / ( s1h + m22 )
     4     + ( 2d0 * Delp2 - 3d0 * Q2IC * ( s1h + Spp) ) )
     5     + Splus * ( - 2d0 * ( Del2 - 6d0 * m12 * Q2IC )
     6     * ( s1h + m22 ) - 2d0 * ( m12 + m22 ) * s1h2
     7     - 9d0 * m22 * Spm**2 + Del2 * ( 2d0 * Spp - m22 )
     8     + 2d0 * s1h * ( 2d0 * Del2 + ( m12 - 5d0 * m22 ) * Spm )
     9     + ( Delp2 - 6d0 * Q2IC * ( m22 + s1h ) )
     1     * Spp * ( s1h + Spp ) / 2d0 / ( s1h + m22 )
     2     - 2d0 * Del2 / s1h * ( Del2
     3     + 2d0 * ( 2d0 * m22 + s1h ) * Spm )
     4     + ( s1h + m22 ) / Delp * ( - 2d0 / s1h * Del2
     5     * ( Del2 + 2d0 * Spm * Spp )
     6     - 2d0 * s1h * ( Del2 - 6d0 * m12 * Q2IC )
     7     - ( Delp2 - 18d0 * m12 * Q2IC ) * Spp
     8     - 2d0 * Del2 * ( Spp + 2d0 * Spm) ) * Lxi ) )
*
      N2 = 2d0 * Splus * Del / Delp2
*
      c21ICR = ( 1d0 - z ) * s1h * f2hat / N2 / 8d0 / ( s1h + m22 )
*
      return
      end
*
************************************************************************
      function c31ICR(z)
*
      implicit none
*
      include "../commons/wrapIC.h"
**
*     Input Variables
*
      double precision z
**
*     Internal Variables
*
      double precision s1h,s1h2
      double precision Lxi,Ixi
      double precision DeltaFun,Delp,Delp2
      double precision f3hat,N3
**
*     Output Variables
*
      double precision c31ICR
*
      s1h   = ( 1d0 - z ) * ( ( Del - Spm ) * z + Del + Spm ) / 2d0 / z
      s1h2  = s1h * s1h
      Delp  = DeltaFun(m12,s1h+m22,-Q2IC)
      Delp2 = Delp * Delp
      Lxi   = dlog( ( Spp + s1h - Delp ) / ( Spp + s1h + Delp ) )
      Ixi   = ( ( s1h + 2d0 * m22 ) / s1h2 + ( s1h + m22 )
     1      / Delp / s1h2 * Spp * Lxi )
*
      f3hat = 16d0 / Delp2 * ( - 2d0 * Del2 * Rplus * Ixi
     1     + 2d0 * m1 * m2 * Rminus * ( 1d0 - Smp / s1h + ( s1h + m22 )
     2     * ( s1h + Spm ) / Delp / s1h * Lxi )
     3     + Rplus * ( Smp - 3d0 * Spm - 2d0 / s1h * ( Del2
     4     + 2d0 * m22 * Spm ) - ( s1h - Smp ) * ( s1h + Spp ) / 2d0
     5     / ( s1h + m22 ) + ( s1h + m22 ) / Delp / s1h * ( - s1h2
     6     + 4d0 * ( m12 * Smp - Del2 ) - 3d0 * s1h * Spm ) * Lxi ) )
*
      N3 = 2d0 * Rplus / Delp
*
      c31ICR = ( 1d0 - z ) * s1h * f3hat / N3 / 8d0 / ( s1h + m22 )
      return
      end
*
************************************************************************
      function cL1ICR(z)
*
      implicit none
*
      include "../commons/wrapIC.h"
**
*     Input Variables
*
      double precision z
**
*     Internal Variables
*
      double precision c11ICR,c21ICR
**
*     Output Variables
*
      double precision cL1ICR
*
      cL1ICR = fact2 * c21ICR(z) - fact1 * c11ICR(z)
*
      return
      end
*
************************************************************************
*     Local terms proportional to delta(1 - z)
*     + the residual term coming frm the plus prescripted term.
************************************************************************
      function c11ICL(z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
      include "../commons/wrapIC.h"
**
*     Input Variables
*
      double precision z
**
*     Internal Variables
*
      double precision c11ICR
**
*     Output Variables
*
      double precision c11ICL
*
      c11ICL = S1 + V1 + c11ICR(one) * dlog( 1d0 - z )
*
      return
      end
*
************************************************************************
      function c21ICL(z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
      include "../commons/wrapIC.h"
**
*     Input Variables
*
      double precision z
**
*     Internal Variables
*
      double precision c21ICR
**
*     Output Variables
*
      double precision c21ICL
*
      c21ICL = S2 + V2 + c21ICR(one) * dlog( 1d0 - z )
*
      return
      end
*
************************************************************************
      function c31ICL(z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
      include "../commons/wrapIC.h"
**
*     Input Variables
*
      double precision z
**
*     Internal Variables
*
      double precision c31ICR
**
*     Output Variables
*
      double precision c31ICL
*
      c31ICL = S3 + V3 + c31ICR(one) * dlog( 1d0 - z )
*
      return
      end
*
************************************************************************
      function cL1ICL(z)
*
      implicit none
*
      include "../commons/wrapIC.h"
**
*     Input Variables
*
      double precision z
**
*     Internal Variables
*
      double precision c11ICL,c21ICL
**
*     Output Variables
*
      double precision cL1ICL
*
      cL1ICL = fact2 * c21ICL(z) - fact1 * c11ICL(z)
*
      return
      end
*
************************************************************************
      function DeltaFun(a,b,c)
*
      implicit none
**
*     Input Variables
*
      double precision a,b,c
**
*     Output Variables
*
      double precision DeltaFun
*
      DeltaFun = dsqrt( a * a + b * b + c * c
     1         - 2d0 * ( a * b + b * c + c * a ) )
*
      return
      end
