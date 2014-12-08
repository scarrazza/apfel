************************************************************************
*
*     FFCoefficientFunctions.f:
*
*     This file contains all the FF coefficient functions as parametrized
*     in the paper written below.
*
************************************************************************
*
*     Initializes the MINIMAX table with the precomputed coeffients
*     of the Alekhin-Blumlein parametrization of the Mellin space heavy
*     quark coefficient functions.
*
*     Phys.Lett.B594:299-307,2004
*     e-Print: hep-ph/0404034
*
************************************************************************
      subroutine initFFNStables
*
      implicit none
*
      include "../commons/coeffhqmellin.h"
*
      include "../commons/minimax.h"
*
c      open (unit=41,file="minimax.tbl",status="old")
c      read (41,fmt="(10i10)") m_coef
c      read (41,fmt="(6e23.16)") ximin,ximax,xigrid,
c     1    coef_p1,coef_p2,coef
c      close(41)
*
      return
      end
*
************************************************************************
*
*     Linear Iterpolation
*
************************************************************************
      FUNCTION CH_FFNS_NC_AB(ICO,X,Q2,M2)
*
      IMPLICIT NONE
*
      include "../commons/coeffhqmellin.h"
**
*     Input Variables
*
      INTEGER ICO,IXI
      DOUBLE PRECISION X,Q2,M2
**
*     Internal Variables
*
      INTEGER I
      DOUBLE PRECISION XI,EPS,LOGXI
      DOUBLE PRECISION DIFF(NXI),SGN
      DOUBLE PRECISION Y,YY(0:1)
      DOUBLE PRECISION TRANS,TRANS_TMP(0:1)
      DOUBLE PRECISION PI
      PARAMETER(PI = 3.14159265358979D0)
**
*     Output Variables
*
      DOUBLE PRECISION CH_FFNS_NC_AB
*
      XI    = Q2 / M2
      LOGXI = DLOG(XI)
      EPS   = 1D0 / XI
*
*     Find IXI such that XIGRID(IXI) < XI < XIGRID(IXI+1)
*
      DIFF(1) = XI - XIGRID(1)
      DO I=2,NXI
         DIFF(I) = XI - XIGRID(I)
         SGN = DIFF(I-1) * DIFF(I)
         IF (SGN.LT.0.D0) THEN
             IXI = I - 1
         ENDIF
      ENDDO
*
*     Check that the value of XI is within the range of the AB CFs or on the borders
*
      IF(XI.LT.XIMIN.OR.XI.GT.XIMAX)THEN
*
*     Set to zero the O(as2) coefficient function if XI
*     is outside the parametrization range
*
         CH_FFNS_NC_AB = 0D0
         RETURN
      ELSEIF(XI.EQ.XIMAX)THEN
         CALL MINIMAX(ICO,NXI,X,TRANS)
         IF(ICO.EQ.1.OR.ICO.EQ.2)THEN
            CH_FFNS_NC_AB = TRANS / ( PI * EPS ) / X       !NLO (LO)
         ELSE
            CH_FFNS_NC_AB = 16D0 * PI * TRANS / EPS / X    !NNLO (NLO)
         ENDIF
         RETURN
      ELSEIF(XI.EQ.XIMIN)THEN
         CALL MINIMAX(ICO,1,X,TRANS)
         IF(ICO.EQ.1.OR.ICO.EQ.2)THEN
            CH_FFNS_NC_AB = TRANS / ( PI * EPS ) / X       !NLO (LO)
         ELSE
            CH_FFNS_NC_AB = 16D0 * PI * TRANS / EPS / X    !NNLO (NLO)
         ENDIF
         RETURN
      ENDIF
*
*     Linear Interpolation
*
      DO I=0,1
         CALL MINIMAX(ICO,IXI+I,X,TRANS_TMP(I))
         YY(I) = DLOG(XIGRID(IXI+I))
      ENDDO
      Y = DLOG(XI)
*
      TRANS = ( TRANS_TMP(1) - TRANS_TMP(0) ) * ( Y - YY(0) ) 
     1      / ( YY(1) - YY(0) ) + TRANS_TMP(0)
*
*     Different normalization according to the perturbative order
*
      IF(ICO.EQ.1.OR.ICO.EQ.2)THEN
         CH_FFNS_NC_AB = TRANS / ( PI * EPS ) / X       !NLO  (LO)
      ELSE
         CH_FFNS_NC_AB = 16D0 * PI * TRANS / EPS / X    !NNLO (NLO)
      ENDIF
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE MINIMAX(ICOEF,IXI,Z,TRANS)
*
      IMPLICIT NONE
*
      include "../commons/coeffhqmellin.h"
**
*     Input Variables
*
      INTEGER ICOEF,IXI
      DOUBLE PRECISION Z
**
*     Internal Variables
*
      INTEGER I
      DOUBLE PRECISION RHO
      DOUBLE PRECISION P1,P2
**
*     Output Variables
*
      DOUBLE PRECISION TRANS
*
      RHO = 1D0 / ( 4D0 / XIGRID(IXI) + 1D0 )
      P1  = COEF_P1(ICOEF)
      P2  = COEF_P2(ICOEF)
*
      TRANS = 0D0
      IF(Z.GE.RHO)THEN
         RETURN
      ELSE
         DO I=0,M_COEF(ICOEF)-1
           TRANS = TRANS + COEF(ICOEF,IXI,I+1) * Z**I
         ENDDO
         TRANS = Z**(-P1) * ( RHO - Z )**(-P2) * TRANS
      ENDIF
*
      RETURN
      END
*
************************************************************************
*
*     Evaluation of the integral in eq. (97) of ArXiv.1001.2312 needed 
*     for the Adler sum rule to be true
*
************************************************************************
      FUNCTION ADLERSR(Q2,M2)
*
      IMPLICIT NONE
*
      include "../commons/vfns.h"
**
*     Input Variables
*
      DOUBLE PRECISION Q2,M2
**
*     Internal Variables
*
      DOUBLE PRECISION DGAUSS,A,B,EPS
      DOUBLE PRECISION LN,ZETA2,CF,TR
      PARAMETER(TR=1D0/2D0)
      PARAMETER(CF=4D0/3D0)

      DOUBLE PRECISION QQ2,MM2
      COMMON/SCALES/QQ2,MM2

      DOUBLE PRECISION L22Q
      EXTERNAL L22Q
**
*     Output Variables
*
      DOUBLE PRECISION ADLERSR
*
      IF(VFNS.EQ."FFNS")THEN
         QQ2 = Q2
         MM2 = M2
*
         A   = 1D-6
         B   = 1D0
         EPS = 1D-5
         ADLERSR = DGAUSS(L22Q,A,B,EPS)
      ELSEIF(VFNS.EQ."FFN0")THEN
         LN = DLOG(Q2/M2)
         ZETA2 = 1.6449340668D0
*     Eq (4.10) of hep-ph/9601302 (Appendix D)
         ADLERSR = CF * TR * ( 2D0 * LN**2D0 - ( 32D0 * ZETA2 / 3D0 
     1           + 38D0 / 3D0 ) * LN + 268D0 * ZETA2 / 9D0 
     2           + 265D0 / 9D0 )
      ENDIF
*
      RETURN
      END
*
************************************************************************
*
*     Reference: Appendix A of hep-ph/9601302
*
************************************************************************
      function l22q(z)
*
      implicit none
**
*     Input Variables
*
      double precision z
**
*     Internal Variables
*
      double precision xi,zmax
      double precision sq1,sq2,l1,l2,l3,dil1,dil2,dil3,dil4
      double precision ddilog

      double precision qq2,mm2
      common/scales/qq2,mm2
**
*     Output Variables
*
      double precision l22q
*
      xi = qq2 / mm2
      zmax = 1d0/(1d0+4d0/xi)
      if(z.ge.zmax)then
         l22q = 0d0
      else
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
         l22q =2d0/3d0*((4d0/3d0*(1d0+z**2d0)/(1d0-z)
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
      endif
*
      return
      end
*
************************************************************************
      function ll2q(z)
*
      implicit none
**
*     Input Variables
*
      double precision z
**
*     Internal Variables
*
      double precision xi,zmax
      double precision sq1,sq2,l1,l2,l3,dil1,dil2,dil3,dil4
      double precision ddilog

      double precision qq2,mm2
      common/scales/qq2,mm2
**
*     Output Variables
*
      double precision ll2q
*
      xi = qq2 / mm2
      zmax = 1d0/(1d0+4d0/xi)
      if(z.ge.zmax)then
         ll2q = 0d0
      else
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
         ll2q =2d0/3d0*(96d0*z*(z/xi)**2d0*(dlog((1d0-z)/z**2d0)*l1
     1        +l1*l2+2d0*(-dil1+dil2+dil3-dil4))+(z/xi/(1d0-z))**2d0
     2        *(64d0-288d0*z+192d0*z**2d0)*l1+z*(16d0/3d0-416d0*z/3d0/xi
     3        +1408d0*z**2d0/3d0/xi**2d0)*l3/sq2+(16d0/3d0-400d0*z/18d0
     4        +z*(-160d0/3d0+3824d0*z/9d0-992d0*z**2d0/3d0)/(1d0-z)/xi)
     5        *sq1)
      endif
*
      return
      end
