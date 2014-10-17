************************************************************************
*
*     MassiveCoefficientFunctions.f:
*
*     This file contains all the ZM coefficient functions.
*
*
*     O(as):
*
*        - icof = 1  F_2,g (NLO, gluon)
*        - icof = 2  F_L,g (NLO, gluon)
*
*     O(as^2):
*
*        - icof = 3  F_2,g (NNLO, gluon)
*        - icof = 4  F_2,q (NNLO, singlet)
*        - icof = 5  F_L,g (NNLO, gluon)
*        - icof = 6  F_L,q (NNLO, singlet)
*
*        - icof = 7  F_2,q (NNLO, non-singlet for the light structure function)
*        - icof = 8  F_L,q (NNLO, non-singlet for the light structure function)
*
*     Coefficients of the terms proportional to ln(q2/m2):
*
*        - icof = 9  F_2,g (NNLO, gluon, bar term)
*        - icof = 10 F_2,q (NNLO, singlet, bar term)
*        - icof = 11 F_L,g (NNLO, gluon, bar term)
*        - icof = 12 F_L,q (NNLO, singlet, bar term)
*
*        - icof = 13 F_2,q (NNLO, non-singlet for the light structure function, bar term)
*        - icof = xx F_L,q (NNLO, non-singlet for the light structure function, bar term) = 0
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
     3      * dlog( ( 1d0 + v ) / ( 1d0 -v ) )
     4      + ( 8d0 * z * (1d0 - z ) - 1d0 
     5      - 4d0 * eps * z * ( 1d0 - z ) ) * v
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













c$$$************************************************************************
c$$$*
c$$$*     Evaluation of the integral in eq. (97) of ArXiv.1001.2312 needed 
c$$$*     for the Adler sum rule to be true
c$$$*
c$$$************************************************************************
c$$$      FUNCTION ADLERSR(Q2,M2)
c$$$*
c$$$      IMPLICIT NONE
c$$$*
c$$$      include "../commons/vfns.h"
c$$$**
c$$$*     Input Variables
c$$$*
c$$$      DOUBLE PRECISION Q2,M2
c$$$**
c$$$*     Internal Variables
c$$$*
c$$$      DOUBLE PRECISION DGAUSS,A,B,EPS
c$$$      DOUBLE PRECISION LN,ZETA2,CF,TR
c$$$      PARAMETER(TR=1D0/2D0)
c$$$      PARAMETER(CF=4D0/3D0)
c$$$
c$$$      DOUBLE PRECISION QQ2,MM2
c$$$      COMMON/SCALES/QQ2,MM2
c$$$
c$$$      DOUBLE PRECISION L22Q
c$$$      EXTERNAL L22Q
c$$$**
c$$$*     Output Variables
c$$$*
c$$$      DOUBLE PRECISION ADLERSR
c$$$*
c$$$      IF(VFNS.EQ."FFNS")THEN
c$$$         QQ2 = Q2
c$$$         MM2 = M2
c$$$*
c$$$         A   = 1D-6
c$$$         B   = 1D0
c$$$         EPS = 1D-5
c$$$         ADLERSR = DGAUSS(L22Q,A,B,EPS)
c$$$      ELSEIF(VFNS.EQ."FFN0")THEN
c$$$         LN = DLOG(Q2/M2)
c$$$         ZETA2 = 1.6449340668D0
c$$$*     Eq (4.10) of hep-ph/9601302 (Appendix D)
c$$$         ADLERSR = CF * TR * ( 2D0 * LN**2D0 - ( 32D0 * ZETA2 / 3D0 
c$$$     1           + 38D0 / 3D0 ) * LN + 268D0 * ZETA2 / 9D0 
c$$$     2           + 265D0 / 9D0 )
c$$$      ENDIF
c$$$*
c$$$      RETURN
c$$$      END
c$$$*
c$$$
c$$$
c$$$
c$$$
c$$$*
c$$$************************************************************************
c$$$*
c$$$*     Order alphas coeficient functions (NLO)
c$$$*     Expansion parameter alphas/4*pi
c$$$*
c$$$************************************************************************
c$$$      subroutine initFFNStables
c$$$*
c$$$      implicit none
c$$$*
c$$$      include "../commons/coeffhqmellin.h"
c$$$      include "../commons/minimax.h"
c$$$*
c$$$      return
c$$$      end
c$$$*
c$$$************************************************************************
c$$$*
c$$$*     Linear Iterpolation
c$$$*
c$$$************************************************************************
c$$$      FUNCTION CH_FFNS_NC_AB(ICO,X,Q2,M2)
c$$$*
c$$$      IMPLICIT NONE
c$$$*
c$$$      include "../commons/coeffhqmellin.h"
c$$$**
c$$$*     Input Variables
c$$$*
c$$$      INTEGER ICO,IXI
c$$$      DOUBLE PRECISION X,Q2,M2
c$$$**
c$$$*     Internal Variables
c$$$*
c$$$      INTEGER I
c$$$      DOUBLE PRECISION XI
c$$$      DOUBLE PRECISION DIFF(NXI),SGN
c$$$      DOUBLE PRECISION MINIMAX
c$$$**
c$$$*     Output Variables
c$$$*
c$$$      DOUBLE PRECISION CH_FFNS_NC_AB
c$$$*
c$$$      XI = Q2 / M2
c$$$*
c$$$*     Find IXI such that XIGRID(IXI) < XI < XIGRID(IXI+1)
c$$$*
c$$$      DIFF(1) = XI - XIGRID(1)
c$$$      DO I=2,NXI
c$$$         DIFF(I) = XI - XIGRID(I)
c$$$         SGN = DIFF(I-1) * DIFF(I)
c$$$         IF (SGN.LT.0.D0) THEN
c$$$             IXI = I - 1
c$$$         ENDIF
c$$$      ENDDO
c$$$*
c$$$*     Check that the value of XI is within the range of the AB CFs or on the borders
c$$$*
c$$$      IF(XI.LT.XIMIN.OR.XI.GT.XIMAX)THEN
c$$$         CH_FFNS_NC_AB = 0D0
c$$$         RETURN
c$$$      ELSEIF(XI.EQ.XIMAX)THEN
c$$$         CH_FFNS_NC_AB = MINIMAX(ICO,NXI,X)
c$$$         RETURN
c$$$      ELSEIF(XI.EQ.XIMIN)THEN
c$$$         CH_FFNS_NC_AB = MINIMAX(ICO,1,X)
c$$$         RETURN
c$$$      ENDIF
c$$$*
c$$$*     Linear Interpolation in log(xi)
c$$$*
c$$$      CH_FFNS_NC_AB = ( MINIMAX(ICO,IXI+1,X) * DLOG(XI/XIGRID(IXI))
c$$$     1              +   MINIMAX(ICO,IXI,X) * DLOG(XIGRID(IXI+1)/XI) )
c$$$     2              / DLOG(XIGRID(IXI+1)/XIGRID(IXI))
c$$$*
c$$$      RETURN
c$$$      END
c$$$*
c$$$
c$$$
c$$$************************************************************************
