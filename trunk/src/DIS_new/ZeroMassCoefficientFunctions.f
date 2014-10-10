************************************************************************
*
*     ZMCoefficientFunctions.f:
*
*     This file contains all the ZM coefficient functions.
*
************************************************************************
*
*     Order alphas coeficient functions (NLO)
*     Expansion parameter alphas/4*pi
*
************************************************************************
*     F2: quark non-singlet - regular term (A)
************************************************************************
      FUNCTION C2NS1A(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Output Variables
*
      DOUBLE PRECISION C2NS1A
*
      C2NS1A = 2D0 * CF * ( - ( 1D0 + X ) * DLOG( 1D0 - X ) 
     1     - ( 1D0 + X**2D0 ) * DLOG(X) / ( 1D0 - X ) + 3D0 + 2D0 * X )
*
      RETURN
      END
*
************************************************************************
*     F2: quark non-singlet - singular term (B)
************************************************************************
      FUNCTION C2NS1B(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Output Variables
*
      DOUBLE PRECISION C2NS1B
*
      C2NS1B = 2D0 * CF * ( 2D0 * DLOG( 1D0 - X ) - 3D0 / 2D0 ) 
     1       / ( 1D0 - X )
*
      RETURN
      END
*
************************************************************************
*     F2: quark non-singlet - local term (C)
************************************************************************
      FUNCTION C2NS1C(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
      INCLUDE "../commons/consts.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Output Variables
*
      DOUBLE PRECISION C2NS1C
*
      C2NS1C = 2D0 * CF * ( DLOG( 1D0 - X )**2D0 
     1       - 3D0 * DLOG( 1D0 - X ) / 2D0 
     2       - ( 2D0 * ZETA2 + 9D0 / 2D0 ) )
*
      RETURN
      END
*
************************************************************************
*     F2: gluon - regular term (A)
************************************************************************
      FUNCTION C2G1A(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Output Variables
*
      DOUBLE PRECISION C2G1A
*
      C2G1A = 2D0 * TR * ( ( ( 1D0 - X )**2D0 + X**2D0 ) 
     1      * DLOG( ( 1D0 - X ) / X ) - 8D0 * X**2D0 + 8D0 * X - 1D0 )
*
      RETURN
      END
*
************************************************************************
*     FL: quark non-singlet - regular term (A)
************************************************************************
      FUNCTION CLNS1A(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Output Variables
*
      DOUBLE PRECISION CLNS1A
*
      CLNS1A = 4D0 * CF * X
*
      RETURN
      END
*
************************************************************************
*     FL: gluon - regular term (A)
************************************************************************
      FUNCTION CLG1A(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Output Variables
*
      DOUBLE PRECISION CLG1A
*
      CLG1A = 8D0 * TR * X * ( 1D0 - X )
*
      RETURN
      END
*
************************************************************************
*     F3: quark non-singlet - regular term (A)
************************************************************************
      FUNCTION C3NS1A(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Output Variables
*
      DOUBLE PRECISION C3NS1A
*
      C3NS1A = 2D0 * CF * ( - ( 1D0 + X ) * DLOG( 1D0 - X ) 
     1     - ( 1D0 + X**2D0 ) * DLOG(X) / ( 1D0 - X ) + 2D0 + X )
*
      RETURN
      END
*
************************************************************************
*     F3: quark non-singlet - singular term (B)
************************************************************************
      FUNCTION C3NS1B(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Output Variables
*
      DOUBLE PRECISION C3NS1B
*
      C3NS1B = 2D0 * CF * ( 2D0 * DLOG( 1D0 - X ) - 3D0 / 2D0 ) 
     1        / ( 1D0 - X )
*
      RETURN
      END
*
************************************************************************
*     F3: quark non-singlet - local term (C)
************************************************************************
      FUNCTION C3NS1C(X)
*
      IMPLICIT NONE
*
      INCLUDE "../commons/ColorFactors.h"
      INCLUDE "../commons/consts.h"
**
*     Input Variables
*
      DOUBLE PRECISION X
**
*     Output Variables
*
      DOUBLE PRECISION C3NS1C
*
      C3NS1C = 2D0 * CF * ( DLOG( 1D0 - X )**2D0 
     1        - 3D0 * DLOG( 1D0 - X ) / 2D0 
     2        - ( 2D0 * ZETA2 + 9D0 / 2D0 ) )
*
      RETURN
      END
*
c$$$************************************************************************
c$$$*
c$$$*     Order alphas^2 coeficient functions (NNLO)
c$$$*     Expansion parameter alphas/4*pi 
c$$$*     Paramtrization by van Neerven and Vogt:
c$$$*     - hep-ph/9907472
c$$$*     - hep-ph/0006154
c$$$*
c$$$************************************************************************
c$$$*     F2: quark non-singlet plus - regular term (A)
c$$$************************************************************************
c$$$       FUNCTION C2NSP2A (Y, NF)
c$$$       IMPLICIT REAL*8 (A-Z)
c$$$       INTEGER NF
c$$$*
c$$$       DL  = LOG (Y)
c$$$       DL1 = LOG (1.-Y)
c$$$*
c$$$       C2NSP2A = 
c$$$     1          - 69.59 - 1008.* Y
c$$$     2          - 2.835 * DL**3 - 17.08 * DL**2 + 5.986 * DL 
c$$$     3          - 17.19 * DL1**3 + 71.08 * DL1**2 - 660.7 * DL1
c$$$     4          - 174.8 * DL * DL1**2 + 95.09 * DL**2 * DL1
c$$$     5        + NF * ( - 5.691 - 37.91 * Y 
c$$$     6          + 2.244 * DL**2 + 5.770 * DL 
c$$$     7          - 1.707 * DL1**2  + 22.95 * DL1
c$$$     8          + 3.036 * DL**2 * DL1 + 17.97 * DL * DL1 )     
c$$$*
c$$$       RETURN
c$$$       END
c$$$*
c$$$************************************************************************
c$$$*     F2: quark non-singlet minus - regular term (A)
c$$$************************************************************************
c$$$       FUNCTION C2NSM2A (Y, NF)
c$$$       IMPLICIT REAL*8 (A-Z)
c$$$       INTEGER NF
c$$$*
c$$$       DL  = LOG (Y)
c$$$       DL1 = LOG (1.-Y)
c$$$*
c$$$       C2NSM2A = 
c$$$     1          - 84.18 - 1010.* Y
c$$$     2          - 3.748 * DL**3 - 19.56 * DL**2 - 1.235 * DL 
c$$$     3          - 17.19 * DL1**3 + 71.08 * DL1**2 - 663.0 * DL1
c$$$     4          - 192.4 * DL * DL1**2 + 80.41 * DL**2 * DL1
c$$$     5        + NF * ( - 5.691 - 37.91 * Y 
c$$$     6          + 2.244 * DL**2 + 5.770 * DL 
c$$$     7          - 1.707 * DL1**2  + 22.95 * DL1
c$$$     8          + 3.036 * DL**2 * DL1 + 17.97 * DL * DL1 )     
c$$$*
c$$$       RETURN
c$$$       END
c$$$*
c$$$************************************************************************
c$$$*     Derivative with respect to NF of C2NSP2A and C2NSM2A
c$$$************************************************************************
c$$$       FUNCTION C2NS2A_DNF (Y)
c$$$       IMPLICIT REAL*8 (A-Z)
c$$$*
c$$$       DL  = LOG (Y)
c$$$       DL1 = LOG (1.-Y)
c$$$*
c$$$       C2NS2A_DNF = 
c$$$     1          - 5.691 - 37.91 * Y 
c$$$     2          + 2.244 * DL**2 + 5.770 * DL 
c$$$     3          - 1.707 * DL1**2  + 22.95 * DL1
c$$$     4          + 3.036 * DL**2 * DL1 + 17.97 * DL * DL1
c$$$*
c$$$       RETURN
c$$$       END
c$$$*
c$$$************************************************************************
c$$$*     F2: quark non-singlet - singular term (B)
c$$$************************************************************************
c$$$       FUNCTION C2NS2B (Y, NF)
c$$$       IMPLICIT REAL*8 (A-Z)
c$$$       INTEGER NF
c$$$*
c$$$       DL1 = LOG (1.-Y)
c$$$       DM  = 1./(1.-Y)
c$$$*
c$$$       C2NS2B = 
c$$$     1          + 14.2222 * DL1**3 - 61.3333 * DL1**2 - 31.105 * DL1 
c$$$     2          + 188.64 
c$$$     3        + NF * ( 1.77778 * DL1**2 - 8.5926 * DL1 + 6.3489 ) 
c$$$       C2NS2B = DM * C2NS2B
c$$$*
c$$$       RETURN
c$$$       END
c$$$*
c$$$************************************************************************
c$$$*     F2: quark non-singlet plus - local term (C)
c$$$************************************************************************
c$$$       FUNCTION C2NSP2C (Y, NF)
c$$$       IMPLICIT REAL*8 (A-Z)
c$$$       INTEGER NF
c$$$*
c$$$       DL1 = LOG (1.-Y)
c$$$*
c$$$       C2NSP2C = 
c$$$     1          + 3.55555 * DL1**4 - 20.4444 * DL1**3 - 15.5525 * DL1**2
c$$$     2          + 188.64 * DL1 - 338.531 + 0.485 
c$$$     3        + NF * (0.592593 * DL1**3 - 4.2963 * DL1**2 
c$$$     4          + 6.3489 * DL1 + 46.844 - 0.0035)
c$$$*
c$$$       RETURN
c$$$       END
c$$$*
c$$$************************************************************************
c$$$*     F2: quark non-singlet minus - local term (C)
c$$$************************************************************************
c$$$       FUNCTION C2NSM2C (Y, NF)
c$$$       IMPLICIT REAL*8 (A-Z)
c$$$       INTEGER NF
c$$$*
c$$$       DL1 = LOG (1.-Y)
c$$$*
c$$$       C2NSM2C = 
c$$$     1          + 3.55555 * DL1**4 - 20.4444 * DL1**3 - 15.5525 * DL1**2
c$$$     2          + 188.64 * DL1 - 338.531 + 0.537
c$$$     3        + NF * (0.592593 * DL1**3 - 4.2963 * DL1**2 
c$$$     4          + 6.3489 * DL1 + 46.844 - 0.0035)
c$$$*
c$$$       RETURN
c$$$       END
c$$$*
c$$$************************************************************************
c$$$*     F2: quark pure-singlet - regular term (A)
c$$$************************************************************************
c$$$       FUNCTION C2PS2A (Y, NF)
c$$$       IMPLICIT REAL*8 (A-Z)
c$$$       INTEGER NF
c$$$*
c$$$       DL  = LOG (Y)
c$$$       DL1 = LOG (1.-Y)
c$$$*
c$$$       C2PS2A =   NF * ( 5.290 * (1./Y-1.) + 4.310 * DL**3   
c$$$     1         - 2.086 * DL**2 + 39.78 * DL - 0.101 * (1.-Y) * DL1**3 
c$$$     2         - (24.75 - 13.80 * Y) * DL**2 * DL1 + 30.23 * DL * DL1 )
c$$$C       C2PS2A = NF * ( ( 8D0 / 3D0 * DL1**2D0 - 32D0 / 3D0 * DL1 
c$$$C     1       + 9.8937D0 ) * ( 1D0 - Y ) + ( 9.57D0 - 13.41D0 * Y 
c$$$C     2       + 0.08D0 * DL1**3D0 ) * ( 1D0 - Y )**2D0 
c$$$C     3       + 5.667D0 * Y * DL**3D0 - DL**2D0 * DL1 * ( 20.26D0 
c$$$C     4       - 33.93D0 * Y ) + 43.36D0 * ( 1D0 - Y ) * DL 
c$$$C     5       - 1.053D0 * DL**2D0 + 40D0 / 9D0 * DL**3D0 
c$$$C     6       + 5.2903D0 * ( 1D0 - Y )**2D0 / Y )
c$$$*
c$$$       RETURN
c$$$       END
c$$$*
c$$$************************************************************************
c$$$*     F2: gluon - regular term (A)
c$$$************************************************************************
c$$$       FUNCTION C2G2A (Y, NF)
c$$$       IMPLICIT REAL*8 (A-Z)
c$$$       INTEGER NF
c$$$*
c$$$       DL  = LOG (Y)
c$$$       DL1 = LOG (1.-Y)
c$$$*
c$$$       C2G2A =   NF * ( 1./Y * (11.90 + 1494.* DL1) + 5.319 * DL**3  
c$$$     1         - 59.48 * DL**2 - 284.8 * DL + 392.4 - 1483.* DL1
c$$$     2         + (6.445 + 209.4 * (1.-Y)) * DL1**3 - 24.00 * DL1**2
c$$$     3         - 724.1 * DL**2 * DL1 - 871.8 * DL * DL1**2 )
c$$$*
c$$$       RETURN
c$$$       END
c$$$*
c$$$************************************************************************
c$$$*     F2: gluon - local term (C) (Articficial term due to the parametrization)
c$$$************************************************************************
c$$$       FUNCTION C2G2C (Y, NF)
c$$$       IMPLICIT REAL*8 (A-Z)
c$$$       INTEGER NF
c$$$*
c$$$       C2G2C = - NF * 0.28  
c$$$*
c$$$       RETURN
c$$$       END
c$$$*
c$$$************************************************************************
c$$$*     FL: quark non-singlet plus - regular term (A)
c$$$************************************************************************
c$$$       FUNCTION CLNSP2A (Y, NF)
c$$$       IMPLICIT REAL*8 (A-Z)
c$$$       INTEGER NF
c$$$*
c$$$       DL  = LOG (Y)
c$$$       DL1 = LOG (1.-Y)
c$$$*
c$$$       CLNSP2A = 
c$$$     1          - 40.41 + 97.48 * Y
c$$$     2          + (26.56 * Y - 0.031) * DL**2 - 14.85 * DL 
c$$$     3          + 13.62 * DL1**2 - 55.79 * DL1 - 150.5 * DL * DL1 
c$$$     4        + NF * 16./27.D0 * ( 6.* Y*DL1 - 12.* Y*DL - 25.* Y + 6.)
c$$$*
c$$$       RETURN
c$$$       END
c$$$*
c$$$************************************************************************
c$$$*     FL: quark non-singlet minus - regular term (A)
c$$$************************************************************************
c$$$       FUNCTION CLNSM2A (Y, NF)
c$$$       IMPLICIT REAL*8 (A-Z)
c$$$       INTEGER NF
c$$$*
c$$$       DL  = LOG (Y)
c$$$       DL1 = LOG (1.-Y)
c$$$*
c$$$       CLNSM2A = 
c$$$     1          - 52.27 + 100.8 * Y
c$$$     2          + (23.29 * Y - 0.043) * DL**2 - 22.21 * DL 
c$$$     3          + 13.30 * DL1**2 - 59.12 * DL1 - 141.7 * DL * DL1 
c$$$     4        + NF * 16./27.D0 * ( 6.* Y*DL1 - 12.* Y*DL - 25.* Y + 6.)
c$$$*
c$$$       RETURN
c$$$       END
c$$$*
c$$$************************************************************************
c$$$*     Derivative with respect to NF of CLNSP2A and CLNSM2A
c$$$************************************************************************
c$$$       FUNCTION CLNS2A_DNF (Y)
c$$$       IMPLICIT REAL*8 (A-Z)
c$$$*
c$$$       DL  = LOG (Y)
c$$$       DL1 = LOG (1.-Y)
c$$$*
c$$$       CLNS2A_DNF = 
c$$$     1        + 16./27.D0 * ( 6.* Y*DL1 - 12.* Y*DL - 25.* Y + 6.)
c$$$*
c$$$       RETURN
c$$$       END
c$$$*
c$$$************************************************************************
c$$$*     FL: quark non-singlet plus - local term (C)
c$$$************************************************************************
c$$$       FUNCTION CLNSP2C (Y)
c$$$       IMPLICIT REAL*8 (A-Z)
c$$$*
c$$$       CLNSP2C = -0.164
c$$$*
c$$$       RETURN
c$$$       END
c$$$*
c$$$************************************************************************
c$$$*     FL: quark non-singlet minus - local term (C)
c$$$************************************************************************
c$$$       FUNCTION CLNSM2C (Y)
c$$$       IMPLICIT REAL*8 (A-Z)
c$$$*
c$$$       CLNSM2C = -0.150
c$$$*
c$$$       RETURN
c$$$       END
c$$$*
c$$$************************************************************************
c$$$*     FL: quark pure-singlet - regular term (A)
c$$$************************************************************************
c$$$       FUNCTION CLPS2A (Y, NF)
c$$$       IMPLICIT REAL*8 (A-Z)
c$$$       INTEGER NF
c$$$*
c$$$       DL  = LOG (Y)
c$$$       DL1 = LOG (1.-Y)
c$$$*
c$$$       CLPS2A = NF * ( (15.94 - 5.212 * Y) * (1.-Y)**2 * DL1
c$$$     1         + (0.421 + 1.520 * Y) * DL**2 + 28.09 * (1.-Y) * DL
c$$$     2         - (2.370/Y - 19.27) * (1.-Y)**3 )
c$$$*
c$$$       RETURN
c$$$       END
c$$$*
c$$$************************************************************************
c$$$*     FL: gluon - regular term (A)
c$$$************************************************************************
c$$$       FUNCTION CLG2A (Y, NF)
c$$$       IMPLICIT REAL*8 (A-Z)
c$$$       INTEGER NF
c$$$*
c$$$       DL  = LOG (Y)
c$$$       DL1 = LOG (1.-Y)
c$$$*
c$$$       CLG2A = NF * ( (94.74 - 49.20 * Y) * (1.-Y) * DL1**2 
c$$$     1         + 864.8 * (1.-Y) * DL1 + 1161.* Y * DL * DL1 
c$$$     2         + 60.06 * Y * DL**2 + 39.66 * (1.-Y) * DL 
c$$$     3         - 5.333 * (1./Y - 1.) )
c$$$*
c$$$       RETURN
c$$$       END
c$$$*
c$$$************************************************************************
c$$$*     F3: quark non-singlet plus - regular term (A)
c$$$************************************************************************
c$$$       FUNCTION C3NSP2A (Y, NF)
c$$$       IMPLICIT REAL*8 (A-Z)
c$$$       INTEGER NF
c$$$*
c$$$       DL  = LOG (Y)
c$$$       DL1 = LOG (1.-Y)
c$$$*
c$$$       C3NSP2A = 
c$$$     1          - 242.9 - 467.2 * Y
c$$$     2          - 3.049 * DL**3 - 30.14 * DL**2 - 79.14 * DL 
c$$$     3          - 15.20 * DL1**3 + 94.61 * DL1**2 - 396.1 * DL1
c$$$     4          - 92.43 * DL * DL1**2 
c$$$     5        + NF * ( - 6.337 - 14.97 * Y 
c$$$     6          + 2.207 * DL**2 + 8.683 * DL 
c$$$     7          + 0.042 * DL1**3 - 0.808 * DL1**2  + 25.00 * DL1
c$$$     8          + 9.684 * DL * DL1 )     
c$$$*
c$$$       RETURN
c$$$       END
c$$$*
c$$$************************************************************************
c$$$*     F3: quark non-singlet minus - regular term (A)
c$$$************************************************************************
c$$$       FUNCTION C3NSM2A (Y, NF)
c$$$       IMPLICIT REAL*8 (A-Z)
c$$$       INTEGER NF
c$$$*
c$$$       DL  = LOG (Y)
c$$$       DL1 = LOG (1.-Y)
c$$$*
c$$$       C3NSM2A = 
c$$$     1          - 206.1 - 576.8 * Y
c$$$     2          - 3.922 * DL**3 - 33.31 * DL**2 - 67.60 * DL 
c$$$     3          - 15.20 * DL1**3 + 94.61 * DL1**2 - 409.6 * DL1
c$$$     4          - 147.9 * DL * DL1**2 
c$$$     5        + NF * ( - 6.337 - 14.97 * Y 
c$$$     6          + 2.207 * DL**2 + 8.683 * DL 
c$$$     7          + 0.042 * DL1**3 - 0.808 * DL1**2 + 25.00 * DL1
c$$$     8          + 9.684 * DL * DL1 )     
c$$$*
c$$$       RETURN
c$$$       END
c$$$*
c$$$************************************************************************
c$$$*     F3: quark non-singlet - singular term (B)
c$$$************************************************************************
c$$$       FUNCTION C3NS2B (Y, NF)
c$$$       IMPLICIT REAL*8 (A-Z)
c$$$       INTEGER NF
c$$$*
c$$$       DL1 = LOG (1.-Y)
c$$$       DM  = 1./(1.-Y)
c$$$*
c$$$       C3NS2B = 
c$$$     1          + 14.2222 * DL1**3 - 61.3333 * DL1**2 - 31.105 * DL1 
c$$$     2          + 188.64 
c$$$     3        + NF * ( 1.77778 * DL1**2 - 8.5926 * DL1 + 6.3489 ) 
c$$$       C3NS2B = DM * C3NS2B
c$$$*
c$$$       RETURN
c$$$       END
c$$$*
c$$$************************************************************************
c$$$*     F3: quark non-singlet plus - local term (C)
c$$$************************************************************************
c$$$       FUNCTION C3NSP2C (Y, NF)
c$$$       IMPLICIT REAL*8 (A-Z)
c$$$       INTEGER NF
c$$$*
c$$$       DL1 = LOG (1.-Y)
c$$$*
c$$$       C3NSP2C = 
c$$$     1          + 3.55555 * DL1**4 - 20.4444 * DL1**3 - 15.5525 * DL1**2
c$$$     2          + 188.64 * DL1 - 338.531  - 0.152 
c$$$     3        + NF * (0.592593 * DL1**3 - 4.2963 * DL1**2 
c$$$     4          + 6.3489 * DL1 + 46.844 + 0.013)
c$$$*
c$$$       RETURN
c$$$       END
c$$$*
c$$$************************************************************************
c$$$*     F3: quark non-singlet minus - local term (C)
c$$$************************************************************************
c$$$       FUNCTION C3NSM2C (Y, NF)
c$$$       IMPLICIT REAL*8 (A-Z)
c$$$       INTEGER NF
c$$$*
c$$$       DL1 = LOG (1.-Y)
c$$$*
c$$$       C3NSM2C = 
c$$$     1          + 3.55555 * DL1**4 - 20.4444 * DL1**3 - 15.5525 * DL1**2
c$$$     2          + 188.64 * DL1 - 338.531 - 0.104 
c$$$     3        + NF * (0.592593 * DL1**3 - 4.2963 * DL1**2 
c$$$     4          + 6.3489 * DL1 + 46.844 + 0.013)
c$$$*
c$$$       RETURN
c$$$       END
