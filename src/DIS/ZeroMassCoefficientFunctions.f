************************************************************************
*
*     ZMCoefficientFunctions.f:
*
*     This file contains all the ZM coefficient functions.
*
************************************************************************
*
*     Order alphas coefficient functions (NLO)
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
      C2G1A = 4D0 * TR * ( ( ( 1D0 - X )**2D0 + X**2D0 ) 
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
      CLG1A = 16D0 * TR * X * ( 1D0 - X )
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
************************************************************************
*
*     Order alphas^2 coefficient functions (NNLO)
*     Expansion parameter alphas/4*pi 
*     Parametrization by van Neerven and Vogt:
*     - hep-ph/9907472
*     - hep-ph/0006154
*
************************************************************************
*     F2: quark non-singlet plus - regular term (A)
************************************************************************
       FUNCTION C2NSP2A (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       C2NSP2A = 
     1          - 69.59 - 1008.* Y
     2          - 2.835 * DL**3 - 17.08 * DL**2 + 5.986 * DL 
     3          - 17.19 * DL1**3 + 71.08 * DL1**2 - 660.7 * DL1
     4          - 174.8 * DL * DL1**2 + 95.09 * DL**2 * DL1
     5        + NF * ( - 5.691 - 37.91 * Y 
     6          + 2.244 * DL**2 + 5.770 * DL 
     7          - 1.707 * DL1**2  + 22.95 * DL1
     8          + 3.036 * DL**2 * DL1 + 17.97 * DL * DL1 )     
*
       RETURN
       END
*
************************************************************************
*     F2: quark non-singlet minus - regular term (A)
************************************************************************
       FUNCTION C2NSM2A (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       C2NSM2A = 
     1          - 84.18 - 1010.* Y
     2          - 3.748 * DL**3 - 19.56 * DL**2 - 1.235 * DL 
     3          - 17.19 * DL1**3 + 71.08 * DL1**2 - 663.0 * DL1
     4          - 192.4 * DL * DL1**2 + 80.41 * DL**2 * DL1
     5        + NF * ( - 5.691 - 37.91 * Y 
     6          + 2.244 * DL**2 + 5.770 * DL 
     7          - 1.707 * DL1**2  + 22.95 * DL1
     8          + 3.036 * DL**2 * DL1 + 17.97 * DL * DL1 )     
*
       RETURN
       END
*
************************************************************************
*     Derivative with respect to NF of C2NSP2A and C2NSM2A
************************************************************************
       FUNCTION C2NS2A_DNF (Y)
       IMPLICIT REAL*8 (A-Z)
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       C2NS2A_DNF = 
     1          - 5.691 - 37.91 * Y 
     2          + 2.244 * DL**2 + 5.770 * DL 
     3          - 1.707 * DL1**2  + 22.95 * DL1
     4          + 3.036 * DL**2 * DL1 + 17.97 * DL * DL1
*
       RETURN
       END
*
************************************************************************
*     F2: quark non-singlet - singular term (B)
************************************************************************
       FUNCTION C2NS2B (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL1 = LOG (1.-Y)
       DM  = 1./(1.-Y)
*
       C2NS2B = 
     1          + 14.2222 * DL1**3 - 61.3333 * DL1**2 - 31.105 * DL1 
     2          + 188.64 
     3        + NF * ( 1.77778 * DL1**2 - 8.5926 * DL1 + 6.3489 ) 
       C2NS2B = DM * C2NS2B
*
       RETURN
       END
*
************************************************************************
*     F2: quark non-singlet plus - local term (C)
************************************************************************
       FUNCTION C2NSP2C (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL1 = LOG (1.-Y)
*
       C2NSP2C = 
     1          + 3.55555 * DL1**4 - 20.4444 * DL1**3 - 15.5525 * DL1**2
     2          + 188.64 * DL1 - 338.531 + 0.485 
     3        + NF * (0.592593 * DL1**3 - 4.2963 * DL1**2 
     4          + 6.3489 * DL1 + 46.844 - 0.0035)
*
       RETURN
       END
*
************************************************************************
*     F2: quark non-singlet minus - local term (C)
************************************************************************
       FUNCTION C2NSM2C (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL1 = LOG (1.-Y)
*
       C2NSM2C = 
     1          + 3.55555 * DL1**4 - 20.4444 * DL1**3 - 15.5525 * DL1**2
     2          + 188.64 * DL1 - 338.531 + 0.537
     3        + NF * (0.592593 * DL1**3 - 4.2963 * DL1**2 
     4          + 6.3489 * DL1 + 46.844 - 0.0035)
*
       RETURN
       END
*
************************************************************************
*     F2: quark pure-singlet - regular term (A)
************************************************************************
       FUNCTION C2PS2A (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       C2PS2A =   NF * ( 5.290 * (1./Y-1.) + 4.310 * DL**3   
     1         - 2.086 * DL**2 + 39.78 * DL - 0.101 * (1.-Y) * DL1**3 
     2         - (24.75 - 13.80 * Y) * DL**2 * DL1 + 30.23 * DL * DL1 )
C       C2PS2A = NF * ( ( 8D0 / 3D0 * DL1**2D0 - 32D0 / 3D0 * DL1 
C     1       + 9.8937D0 ) * ( 1D0 - Y ) + ( 9.57D0 - 13.41D0 * Y 
C     2       + 0.08D0 * DL1**3D0 ) * ( 1D0 - Y )**2D0 
C     3       + 5.667D0 * Y * DL**3D0 - DL**2D0 * DL1 * ( 20.26D0 
C     4       - 33.93D0 * Y ) + 43.36D0 * ( 1D0 - Y ) * DL 
C     5       - 1.053D0 * DL**2D0 + 40D0 / 9D0 * DL**3D0 
C     6       + 5.2903D0 * ( 1D0 - Y )**2D0 / Y )
*
       RETURN
       END
*
************************************************************************
*     F2: gluon - regular term (A)
************************************************************************
       FUNCTION C2G2A (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       C2G2A =   NF * ( 1./Y * (11.90 + 1494.* DL1) + 5.319 * DL**3  
     1         - 59.48 * DL**2 - 284.8 * DL + 392.4 - 1483.* DL1
     2         + (6.445 + 209.4 * (1.-Y)) * DL1**3 - 24.00 * DL1**2
     3         - 724.1 * DL**2 * DL1 - 871.8 * DL * DL1**2 )
*
       RETURN
       END
*
************************************************************************
*     F2: gluon - local term (C) (Articficial term due to the parametrization)
************************************************************************
       FUNCTION C2G2C (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       C2G2C = - NF * 0.28  
*
       RETURN
       END
*
************************************************************************
*     FL: quark non-singlet plus - regular term (A)
************************************************************************
       FUNCTION CLNSP2A (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       CLNSP2A = 
     1          - 40.41 + 97.48 * Y
     2          + (26.56 * Y - 0.031) * DL**2 - 14.85 * DL 
     3          + 13.62 * DL1**2 - 55.79 * DL1 - 150.5 * DL * DL1 
     4        + NF * 16./27.D0 * ( 6.* Y*DL1 - 12.* Y*DL - 25.* Y + 6.)
*
       RETURN
       END
*
************************************************************************
*     FL: quark non-singlet minus - regular term (A)
************************************************************************
       FUNCTION CLNSM2A (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       CLNSM2A = 
     1          - 52.27 + 100.8 * Y
     2          + (23.29 * Y - 0.043) * DL**2 - 22.21 * DL 
     3          + 13.30 * DL1**2 - 59.12 * DL1 - 141.7 * DL * DL1 
     4        + NF * 16./27.D0 * ( 6.* Y*DL1 - 12.* Y*DL - 25.* Y + 6.)
*
       RETURN
       END
*
************************************************************************
*     Derivative with respect to NF of CLNSP2A and CLNSM2A
************************************************************************
       FUNCTION CLNS2A_DNF (Y)
       IMPLICIT REAL*8 (A-Z)
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       CLNS2A_DNF = 
     1        + 16./27.D0 * ( 6.* Y*DL1 - 12.* Y*DL - 25.* Y + 6.)
*
       RETURN
       END
*
************************************************************************
*     FL: quark non-singlet plus - local term (C)
************************************************************************
       FUNCTION CLNSP2C (Y)
       IMPLICIT REAL*8 (A-Z)
*
       CLNSP2C = -0.164
*
       RETURN
       END
*
************************************************************************
*     FL: quark non-singlet minus - local term (C)
************************************************************************
       FUNCTION CLNSM2C (Y)
       IMPLICIT REAL*8 (A-Z)
*
       CLNSM2C = -0.150
*
       RETURN
       END
*
************************************************************************
*     FL: quark pure-singlet - regular term (A)
************************************************************************
       FUNCTION CLPS2A (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       CLPS2A = NF * ( (15.94 - 5.212 * Y) * (1.-Y)**2 * DL1
     1         + (0.421 + 1.520 * Y) * DL**2 + 28.09 * (1.-Y) * DL
     2         - (2.370/Y - 19.27) * (1.-Y)**3 )
*
       RETURN
       END
*
************************************************************************
*     FL: gluon - regular term (A)
************************************************************************
       FUNCTION CLG2A (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       CLG2A = NF * ( (94.74 - 49.20 * Y) * (1.-Y) * DL1**2 
     1         + 864.8 * (1.-Y) * DL1 + 1161.* Y * DL * DL1 
     2         + 60.06 * Y * DL**2 + 39.66 * (1.-Y) * DL 
     3         - 5.333 * (1./Y - 1.) )
*
       RETURN
       END
*
************************************************************************
*     F3: quark non-singlet plus - regular term (A)
************************************************************************
       FUNCTION C3NSP2A (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       C3NSP2A = 
     1          - 242.9 - 467.2 * Y
     2          - 3.049 * DL**3 - 30.14 * DL**2 - 79.14 * DL 
     3          - 15.20 * DL1**3 + 94.61 * DL1**2 - 396.1 * DL1
     4          - 92.43 * DL * DL1**2 
     5        + NF * ( - 6.337 - 14.97 * Y 
     6          + 2.207 * DL**2 + 8.683 * DL 
     7          + 0.042 * DL1**3 - 0.808 * DL1**2  + 25.00 * DL1
     8          + 9.684 * DL * DL1 )     
*
       RETURN
       END
*
************************************************************************
*     F3: quark non-singlet minus - regular term (A)
************************************************************************
       FUNCTION C3NSM2A (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       C3NSM2A = 
     1          - 206.1 - 576.8 * Y
     2          - 3.922 * DL**3 - 33.31 * DL**2 - 67.60 * DL 
     3          - 15.20 * DL1**3 + 94.61 * DL1**2 - 409.6 * DL1
     4          - 147.9 * DL * DL1**2 
     5        + NF * ( - 6.337 - 14.97 * Y 
     6          + 2.207 * DL**2 + 8.683 * DL 
     7          + 0.042 * DL1**3 - 0.808 * DL1**2 + 25.00 * DL1
     8          + 9.684 * DL * DL1 )     
*
       RETURN
       END
*
************************************************************************
*     F3: quark non-singlet - singular term (B)
************************************************************************
       FUNCTION C3NS2B (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL1 = LOG (1.-Y)
       DM  = 1./(1.-Y)
*
       C3NS2B = 
     1          + 14.2222 * DL1**3 - 61.3333 * DL1**2 - 31.105 * DL1 
     2          + 188.64 
     3        + NF * ( 1.77778 * DL1**2 - 8.5926 * DL1 + 6.3489 ) 
       C3NS2B = DM * C3NS2B
*
       RETURN
       END
*
************************************************************************
*     F3: quark non-singlet plus - local term (C)
************************************************************************
       FUNCTION C3NSP2C (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL1 = LOG (1.-Y)
*
       C3NSP2C = 
     1          + 3.55555 * DL1**4 - 20.4444 * DL1**3 - 15.5525 * DL1**2
     2          + 188.64 * DL1 - 338.531  - 0.152 
     3        + NF * (0.592593 * DL1**3 - 4.2963 * DL1**2 
     4          + 6.3489 * DL1 + 46.844 + 0.013)
*
       RETURN
       END
*
************************************************************************
*     F3: quark non-singlet minus - local term (C)
************************************************************************
       FUNCTION C3NSM2C (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL1 = LOG (1.-Y)
*
       C3NSM2C = 
     1          + 3.55555 * DL1**4 - 20.4444 * DL1**3 - 15.5525 * DL1**2
     2          + 188.64 * DL1 - 338.531 - 0.104 
     3        + NF * (0.592593 * DL1**3 - 4.2963 * DL1**2 
     4          + 6.3489 * DL1 + 46.844 + 0.013)
*
       RETURN
       END
*
************************************************************************
*
*     Order alphas coefficient functions (NLO) for the semi-inclusive
*     e+e- annihilation. Expansion parameter alphas/4*pi
*     Reference: hep-ph/0604160
*
************************************************************************
*     F2: quark non-singlet - regular term (A)
************************************************************************
      FUNCTION C2NS1TA(X)
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
      DOUBLE PRECISION C2NS1TA
*
      C2NS1TA = 2D0 * CF * ( - ( 1D0 + X ) * DLOG( 1D0 - X ) 
     1     + 2D0 * ( 1D0 + X**2D0 ) * DLOG(X) / ( 1D0 - X ) + 5D0 / 2D0 
     2     - 3D0 * X / 2D0 )
*
      RETURN
      END
*
************************************************************************
*     F2: quark non-singlet - singular term (B)
************************************************************************
      FUNCTION C2NS1TB(X)
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
      DOUBLE PRECISION C2NS1TB
*
      C2NS1TB = 2D0 * CF * ( 2D0 * DLOG( 1D0 - X ) - 3D0 / 2D0 ) 
     1        / ( 1D0 - X )
*
      RETURN
      END
*
************************************************************************
*     F2: quark non-singlet - local term (C)
************************************************************************
      FUNCTION C2NS1TC(X)
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
      DOUBLE PRECISION C2NS1TC
*
      C2NS1TC = 2D0 * CF * ( DLOG( 1D0 - X )**2D0 
     1        - 3D0 * DLOG( 1D0 - X ) / 2D0 
     2        + ( 4D0 * ZETA2 - 9D0 / 2D0 ) )
*
      RETURN
      END
*
************************************************************************
*     F2: gluon - regular term (A)
************************************************************************
      FUNCTION C2G1TA(X)
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
      DOUBLE PRECISION C2G1TA
*
      C2G1TA = 4D0 * CF * ( ( 1D0 + ( 1D0 - X )**2D0 ) 
     1       * DLOG( X**2D0 * ( 1D0 - X ) ) / X )
*
      RETURN
      END
*
************************************************************************
*     FL: quark non-singlet - regular term (A)
************************************************************************
      FUNCTION CLNS1TA(X)
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
      DOUBLE PRECISION CLNS1TA
*
      CLNS1TA = 2D0 * CF
*
      RETURN
      END
*
************************************************************************
*     FL: gluon - regular term (A)
************************************************************************
      FUNCTION CLG1TA(X)
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
      DOUBLE PRECISION CLG1TA
*
      CLG1TA = 8D0 * CF * ( 1D0 - X ) / X
*
      RETURN
      END
*
************************************************************************
*     F3: quark non-singlet - regular term (A)
************************************************************************
      FUNCTION C3NS1TA(X)
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
      DOUBLE PRECISION C3NS1TA
*
      C3NS1TA = 2D0 * CF * ( - ( 1D0 + X ) * DLOG( 1D0 - X ) 
     1     + 2d0 * ( 1D0 + X**2D0 ) * DLOG(X) / ( 1D0 - X ) + 1D0 / 2D0
     2     - X / 2D0 )
*
      RETURN
      END
*
************************************************************************
*     F3: quark non-singlet - singular term (B)
************************************************************************
      FUNCTION C3NS1TB(X)
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
      DOUBLE PRECISION C3NS1TB
*
      C3NS1TB = 2D0 * CF * ( 2D0 * DLOG( 1D0 - X ) - 3D0 / 2D0 ) 
     1        / ( 1D0 - X )
*
      RETURN
      END
*
************************************************************************
*     F3: quark non-singlet - local term (C)
************************************************************************
      FUNCTION C3NS1TC(X)
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
      DOUBLE PRECISION C3NS1TC
*
      C3NS1TC = 2D0 * CF * ( DLOG( 1D0 - X )**2D0 
     1        - 3D0 * DLOG( 1D0 - X ) / 2D0 
     2        + ( 4D0 * ZETA2 - 9D0 / 2D0 ) )
*
      RETURN
      END
*
************************************************************************
*
*     Order alphas coefficient functions (NNLO) for the semi-inclusive
*     e+e- annihilation. Expansion parameter alphas/4*pi
*     Reference: hep-ph/0604160
*
************************************************************************
*     F2: quark non-singlet plus - regular term (A)
************************************************************************
      FUNCTION C2NSP2TA (X, NF)
      IMPLICIT REAL*8 (A-Z)
*
      INCLUDE "../commons/ColorFactors.h"
      INCLUDE "../commons/consts.h"
*
      COMPLEX*16 HC1, HC2, HC3, HC4 
      INTEGER NF, NF2, N1, N2, NW, I1, I2, I3, N
      PARAMETER ( N1 = -1, N2 = 1, NW = 3 ) 
      DIMENSION HC1(N1:N2),HC2(N1:N2,N1:N2),HC3(N1:N2,N1:N2,N1:N2), 
     ,          HC4(N1:N2,N1:N2,N1:N2,N1:N2) 
      DIMENSION HR1(N1:N2),HR2(N1:N2,N1:N2),HR3(N1:N2,N1:N2,N1:N2), 
     ,          HR4(N1:N2,N1:N2,N1:N2,N1:N2) 
      DIMENSION HI1(N1:N2),HI2(N1:N2,N1:N2),HI3(N1:N2,N1:N2,N1:N2), 
     ,          HI4(N1:N2,N1:N2,N1:N2,N1:N2) 
*
* ...Some abbreviations
*
      z2 = zeta2
      z3 = zeta3
      DX = 1.D0/X
      DM = 1.D0/(1.D0-X)
      DP = 1.D0/(1.D0+X)
      DL1 = LOG (1.D0-X)
*
* ...The harmonic polylogs up to weight 4 by Gehrmann and Remiddi
*
      CALL HPLOG (X, NW, HC1,HC2,HC3,HC4, HR1,HR2,HR3,HR4,
     ,            HI1,HI2,HI3,HI4, N1, N2) 
*
* ...The coefficient function in terms of the harmonic polylogs
*    (without the delta(1-x) part, but with the soft contribution)
*
      CTeq2 =
     &  + cf*ca * ( 1271.D0/270.D0 + 4829.D0/270.D0*x - 24.D0/5.D0*x**2
     &     + 24.D0/5.D0*dx - 3155.D0/54.D0*dm - 36.D0*z3*x - 28.D0*z3*
     &    dp + 28.D0*z3*dm - 4.D0*z2*x + 24.D0/5.D0*z2*x**3 - 12.D0*
     &    Hr1(-1)*z2 + 20.D0*Hr1(-1)*z2*x + 32.D0*Hr1(-1)*z2*dp - 359.D0
     &    /15.D0*Hr1(0) - 143.D0/5.D0*Hr1(0)*x - 24.D0/5.D0*Hr1(0)*x**2
     &     - 24.D0/5.D0*Hr1(0)*dx + 8.D0*Hr1(0)*dp + 206.D0/3.D0*Hr1(0)
     &    *dm + 12.D0*Hr1(0)*z2 + 20.D0*Hr1(0)*z2*x + 8.D0*Hr1(0)*z2*dp
     &     - 32.D0*Hr1(0)*z2*dm - 25.D0/9.D0*Hr1(1) + 311.D0/9.D0*Hr1(1
     &    )*x - 367.D0/9.D0*Hr1(1)*dm + 8.D0*Hr1(1)*z2 - 8.D0*Hr1(1)*z2
     &    *dm - 4.D0*Hr2(-1,0) - 4.D0*Hr2(-1,0)*x + 24.D0/5.D0*Hr2(-1,0
     &    )*x**3 + 24.D0/5.D0*Hr2(-1,0)*dx**2 + 11.D0/3.D0*Hr2(0,0) +
     &    23.D0/3.D0*Hr2(0,0)*x - 24.D0/5.D0*Hr2(0,0)*x**3 - 22.D0/3.D0
     &    *Hr2(0,0)*dm - 22.D0/3.D0*Hr2(0,1) - 22.D0/3.D0*Hr2(0,1)*x +
     &    44.D0/3.D0*Hr2(0,1)*dm + 22.D0/3.D0*Hr2(1,1) + 22.D0/3.D0*
     &    Hr2(1,1)*x - 44.D0/3.D0*Hr2(1,1)*dm - 8.D0*Hr3(-1,-1,0) + 24.D
     &    0*Hr3(-1,-1,0)*x )
      CTeq2 = CTeq2 + cf*ca * ( 32.D0*Hr3(-1,-1,0)*dp - 16.D0*Hr3(-1,0,
     &    0) + 8.D0*Hr3(-1,0,0)*x + 24.D0*Hr3(-1,0,0)*dp + 8.D0*Hr3(-1,
     &    0,1) - 8.D0*Hr3(-1,0,1)*x - 16.D0*Hr3(-1,0,1)*dp + 16.D0*Hr3(
     &    0,-1,0)*x + 8.D0*Hr3(0,-1,0)*dp - 24.D0*Hr3(0,-1,0)*dm - 36.D0
     &    *Hr3(0,0,0)*x - 36.D0*Hr3(0,0,0)*dp + 36.D0*Hr3(0,0,0)*dm - 4.
     &    D0*Hr3(0,0,1) + 4.D0*Hr3(0,0,1)*x + 8.D0*Hr3(0,0,1)*dp + 4.D0
     &    *Hr3(0,1,0) + 4.D0*Hr3(0,1,0)*x - 8.D0*Hr3(0,1,0)*dm + 4.D0*
     &    Hr3(1,0,0) + 12.D0*Hr3(1,0,0)*x - 16.D0*Hr3(1,0,0)*dm - 4.D0*
     &    Hr3(1,0,1) - 4.D0*Hr3(1,0,1)*x + 8.D0*Hr3(1,0,1)*dm + 4.D0*
     &    Hr3(1,1,0) + 4.D0*Hr3(1,1,0)*x - 8.D0*Hr3(1,1,0)*dm )
      CTeq2 = CTeq2 + cf**2 * ( 279.D0/10.D0 - 279.D0/10.D0*x + 48.D0/5.
     &    D0*x**2 - 48.D0/5.D0*dx + 51.D0/2.D0*dm + 56.D0*z3 + 128.D0*
     &    z3*x + 56.D0*z3*dp - 152.D0*z3*dm - 24.D0*z2 - 48.D0/5.D0*z2*
     &    x**3 + 12.D0*z2*dm + 24.D0*Hr1(-1)*z2 - 40.D0*Hr1(-1)*z2*x -
     &    64.D0*Hr1(-1)*z2*dp + 376.D0/5.D0*Hr1(0) + 166.D0/5.D0*Hr1(0)
     &    *x + 48.D0/5.D0*Hr1(0)*x**2 + 48.D0/5.D0*Hr1(0)*dx - 16.D0*
     &    Hr1(0)*dp - 106.D0*Hr1(0)*dm - 4.D0*Hr1(0)*z2 - 20.D0*Hr1(0)*
     &    z2*x - 16.D0*Hr1(0)*z2*dp + 40.D0*Hr1(0)*z2*dm + 13.D0*Hr1(1)
     &     - 51.D0*Hr1(1)*x + 27.D0*Hr1(1)*dm - 8.D0*Hr1(1)*z2 + 8.D0*
     &    Hr1(1)*z2*x + 8.D0*Hr2(-1,0) + 8.D0*Hr2(-1,0)*x - 48.D0/5.D0*
     &    Hr2(-1,0)*x**3 - 48.D0/5.D0*Hr2(-1,0)*dx**2 - 66.D0*Hr2(0,0)
     &     - 30.D0*Hr2(0,0)*x + 48.D0/5.D0*Hr2(0,0)*x**3 + 66.D0*Hr2(0,
     &    0)*dm + 12.D0*Hr2(0,1) - 4.D0*Hr2(0,1)*x + 12.D0*Hr2(0,1)*dm
     &     - 28.D0*Hr2(1,0) + 28.D0*Hr2(1,0)*x + 24.D0*Hr2(1,0)*dm + 16.
     &    D0*Hr2(1,1) + 8.D0*Hr2(1,1)*x - 36.D0*Hr2(1,1)*dm + 16.D0*
     &    Hr3(-1,-1,0) )
      CTeq2 = CTeq2 + cf**2 * (  - 48.D0*Hr3(-1,-1,0)*x - 64.D0*Hr3(-1,
     &    -1,0)*dp + 32.D0*Hr3(-1,0,0) - 16.D0*Hr3(-1,0,0)*x - 48.D0*
     &    Hr3(-1,0,0)*dp - 16.D0*Hr3(-1,0,1) + 16.D0*Hr3(-1,0,1)*x + 32.
     &    D0*Hr3(-1,0,1)*dp - 32.D0*Hr3(0,-1,0)*x - 16.D0*Hr3(0,-1,0)*
     &    dp + 48.D0*Hr3(0,-1,0)*dm + 66.D0*Hr3(0,0,0) + 138.D0*Hr3(0,0
     &    ,0)*x + 72.D0*Hr3(0,0,0)*dp - 160.D0*Hr3(0,0,0)*dm - 8.D0*
     &    Hr3(0,0,1) - 24.D0*Hr3(0,0,1)*x - 16.D0*Hr3(0,0,1)*dp + 8.D0*
     &    Hr3(0,0,1)*dm + 36.D0*Hr3(0,1,0) + 36.D0*Hr3(0,1,0)*x - 72.D0
     &    *Hr3(0,1,0)*dm - 16.D0*Hr3(0,1,1) - 16.D0*Hr3(0,1,1)*x + 40.D0
     &    *Hr3(0,1,1)*dm - 12.D0*Hr3(1,0,0) - 28.D0*Hr3(1,0,0)*x + 40.D0
     &    *Hr3(1,0,0)*dm - 16.D0*Hr3(1,0,1) - 16.D0*Hr3(1,0,1)*x + 32.D0
     &    *Hr3(1,0,1)*dm - 24.D0*Hr3(1,1,0) - 24.D0*Hr3(1,1,0)*x + 48.D0
     &    *Hr3(1,1,0)*dm + 24.D0*Hr3(1,1,1) + 24.D0*Hr3(1,1,1)*x - 48.D0
     &    *Hr3(1,1,1)*dm )
      CTeq2 = CTeq2 + nf*cf * (  - 59.D0/27.D0 - 17.D0/27.D0*x + 247.D0/
     &    27.D0*dm + 10.D0/3.D0*Hr1(0) + 6.D0*Hr1(0)*x - 32.D0/3.D0*
     &    Hr1(0)*dm - 14.D0/9.D0*Hr1(1) - 26.D0/9.D0*Hr1(1)*x + 58.D0/9.
     &    D0*Hr1(1)*dm - 2.D0/3.D0*Hr2(0,0) - 2.D0/3.D0*Hr2(0,0)*x + 4.D
     &    0/3.D0*Hr2(0,0)*dm + 4.D0/3.D0*Hr2(0,1) + 4.D0/3.D0*Hr2(0,1)*
     &    x - 8.D0/3.D0*Hr2(0,1)*dm - 4.D0/3.D0*Hr2(1,1) - 4.D0/3.D0*
     &    Hr2(1,1)*x + 8.D0/3.D0*Hr2(1,1)*dm )
*
* ...The soft (`+'-distribution) part of the coefficient function
*
      A3 = 
     &     + 8.D0*cf**2
      A2 =
     &     - 22.D0/3.D0*ca*cf
     &     - 18.D0*cf**2
     &     + 4.D0/3.D0*cf*nf
      A1 =
     &     - 8.D0*z2*ca*cf
     &     + 16.D0*z2*cf**2
     &     + 367.D0/9.D0*ca*cf
     &     - 27.D0*cf**2
     &     - 58.D0/9.D0*cf*nf
      A0 =
     &     + 44.D0/3.D0*z2*ca*cf
     &     + 40.D0*z3*ca*cf
     &     - 8.D0*z3*cf**2
     &     - 3155.D0/54.D0*ca*cf
     &     + 51.D0/2.D0*cf**2
     &     + 247.D0/27.D0*cf*nf
     &     - 8.D0/3.D0*z2*cf*nf
*
      CTeq2L = DM * ( DL1**3 * A3 + DL1**2 * A2 + DL1 * A1 + A0) 
*
* ...The regular piece of the coefficient function
*
      C2NSP2TA = CTeq2 - CTeq2L + CLNSP2TA (X, NF)
*
      RETURN
      END
*
************************************************************************
*     F2: quark non-singlet - singular term (B)
************************************************************************
      FUNCTION C2NS2TB (X, NF)
      IMPLICIT REAL*8 (A-Z)
*
      INCLUDE "../commons/ColorFactors.h"
      INCLUDE "../commons/consts.h"
*
      INTEGER NF
*
      z2 = zeta2
      z3 = zeta3
*
      A3 = 
     &     + 8.D0*cf**2
      A2 =
     &     - 22.D0/3.D0*ca*cf
     &     - 18.D0*cf**2
     &     + 4.D0/3.D0*cf*nf
      A1 =
     &     - 8.D0*z2*ca*cf
     &     + 16.D0*z2*cf**2
     &     + 367.D0/9.D0*ca*cf
     &     - 27.D0*cf**2
     &     - 58.D0/9.D0*cf*nf
      A0 =
     &     + 44.D0/3.D0*z2*ca*cf
     &     + 40.D0*z3*ca*cf
     &     - 8.D0*z3*cf**2
     &     - 3155.D0/54.D0*ca*cf
     &     + 51.D0/2.D0*cf**2
     &     + 247.D0/27.D0*cf*nf
     &     - 8.D0/3.D0*z2*cf*nf
*
      DL1 = LOG (1.D0-X)
      DM  = 1.D0/(1.D0-X)
*
      C2NS2TB = DM * ( DL1**3 * A3 + DL1**2 * A2 + DL1 * A1 + A0)
*
      RETURN
      END
*
************************************************************************
*     F2: quark non-singlet plus - local term (C)
************************************************************************
      FUNCTION C2NSP2TC (X, NF)
      IMPLICIT REAL*8 (A-Z)
*
      INCLUDE "../commons/ColorFactors.h"
      INCLUDE "../commons/consts.h"
*
      INTEGER NF, NF2
*
      z2 = zeta2
      z3 = zeta3
*
      A3 = 
     &     + 8.D0*cf**2
      A2 =
     &     - 22.D0/3.D0*ca*cf
     &     - 18.D0*cf**2
     &     + 4.D0/3.D0*cf*nf
      A1 =
     &     - 8.D0*z2*ca*cf
     &     + 16.D0*z2*cf**2
     &     + 367.D0/9.D0*ca*cf
     &     - 27.D0*cf**2
     &     - 58.D0/9.D0*cf*nf
      A0 =
     &     + 44.D0/3.D0*z2*ca*cf
     &     + 40.D0*z3*ca*cf
     &     - 8.D0*z3*cf**2
     &     - 3155.D0/54.D0*ca*cf
     &     + 51.D0/2.D0*cf**2
     &     + 247.D0/27.D0*cf*nf
     &     - 8.D0/3.D0*z2*cf*nf
*
* ...The coefficient of delta(1-x)
*
      C2DELT = 
     ,    + ca*cf * ( - 5465.D0/72.D0 + 140.D0/3.D0*z3 + 215.D0/3.D0*z2
     ,                - 49.D0/5.D0*z2**2 )
     ,    + cf**2 * ( 331.D0/8.D0 - 78.D0*z3 - 39.D0*z2 + 30.D0*z2**2 )
     ,    + cf*nf * ( 457.D0/36.D0 + 4.D0/3.D0*z3 - 38.D0/3.D0*z2 )
*
      DL1 = LOG (1.D0-X)
*
      C2NSP2TC =   DL1**4 * A3/4.D0 + DL1**3 * A2/3.D0 
     ,           + DL1**2 * A1/2.D0 + DL1 * A0 + C2DELT
*
      RETURN
      END
*
************************************************************************
*     F2: quark pure-singlet - regular term (A)
************************************************************************
      FUNCTION C2PS2TA (X, NF)
      IMPLICIT REAL*8 (A-Z)
*
      INCLUDE "../commons/ColorFactors.h"
      INCLUDE "../commons/consts.h"
*
      COMPLEX*16 HC1, HC2, HC3, HC4, HC5
      INTEGER NF, NF2, N1, N2, NW
      PARAMETER ( N1 = -1, N2 = 1, NW = 3 )
      DIMENSION HC1(N1:N2),HC2(N1:N2,N1:N2),HC3(N1:N2,N1:N2,N1:N2),
     ,          HC4(N1:N2,N1:N2,N1:N2,N1:N2) 
      DIMENSION HR1(N1:N2),HR2(N1:N2,N1:N2),HR3(N1:N2,N1:N2,N1:N2),
     ,          HR4(N1:N2,N1:N2,N1:N2,N1:N2) 
      DIMENSION HI1(N1:N2),HI2(N1:N2,N1:N2),HI3(N1:N2,N1:N2,N1:N2),
     ,          HI4(N1:N2,N1:N2,N1:N2,N1:N2) 
*
      z2 = zeta2
      z3 = zeta3
      DX = 1.D0/X
*
* ...Harmonic polylogs (HPLs) up to weight NW=3 by Gehrmann and Remiddi 
*
      CALL HPLOG (X, NW, HC1,HC2,HC3,HC4, HR1,HR2,HR3,HR4,
     ,            HI1,HI2,HI3,HI4, N1, N2)
*
* ...The coefficient function in terms of the HPLs
*
      cTeqps2 =
     &  + nf*cf * (  - 118.D0/3.D0 + 70.D0/3.D0*x + 512.D0/27.D0*x**2
     &     - 80.D0/27.D0*dx + 16*z3 + 16*z3*x - 8*z2 - 32*z2*x - 200.D0/
     &    3.D0*Hr1(0) - 104.D0/3.D0*Hr1(0)*x - 128.D0/9.D0*Hr1(0)*x**2
     &     - 16.D0/3.D0*Hr1(0)*dx + 16*Hr1(0)*z2 + 16*Hr1(0)*z2*x + 92.D
     &    0/3.D0*Hr1(1) - 68.D0/3.D0*Hr1(1)*x - 32.D0/3.D0*Hr1(1)*x**2
     &     + 8.D0/3.D0*Hr1(1)*dx - 16*Hr2(-1,0) - 16*Hr2(-1,0)*x - 16.D0
     &    /3.D0*Hr2(-1,0)*x**2 - 16.D0/3.D0*Hr2(-1,0)*dx - 14*Hr2(0,0)
     &     - 14*Hr2(0,0)*x + 16.D0/3.D0*Hr2(0,0)*x**2 + 64.D0/3.D0*Hr2(
     &    0,0)*dx + 4*Hr2(0,1) + 20*Hr2(0,1)*x + 16.D0/3.D0*Hr2(0,1)*
     &    x**2 - 32.D0/3.D0*Hr2(0,1)*dx + 4*Hr2(1,1) - 4*Hr2(1,1)*x -
     &    16.D0/3.D0*Hr2(1,1)*x**2 + 16.D0/3.D0*Hr2(1,1)*dx + 44*Hr3(0,
     &    0,0) + 44*Hr3(0,0,0)*x - 24*Hr3(0,0,1) - 24*Hr3(0,0,1)*x + 8*
     &    Hr3(0,1,1) + 8*Hr3(0,1,1)*x )
*
      C2PS2TA = cTeqps2 + CLPS2TA (X, NF)
*
      RETURN
      END
*
************************************************************************
*     F2: gluon - regular term (A)
************************************************************************
      FUNCTION C2G2TA (X, NF)
      IMPLICIT REAL*8 (A-Z)
*
      INCLUDE "../commons/ColorFactors.h"
      INCLUDE "../commons/consts.h"
*
      COMPLEX*16 HC1, HC2, HC3, HC4, HC5
      INTEGER NF, NF2, N1, N2, NW
      PARAMETER ( N1 = -1, N2 = 1, NW = 3 )
      DIMENSION HC1(N1:N2),HC2(N1:N2,N1:N2),HC3(N1:N2,N1:N2,N1:N2),
     ,          HC4(N1:N2,N1:N2,N1:N2,N1:N2) 
      DIMENSION HR1(N1:N2),HR2(N1:N2,N1:N2),HR3(N1:N2,N1:N2,N1:N2),
     ,          HR4(N1:N2,N1:N2,N1:N2,N1:N2) 
      DIMENSION HI1(N1:N2),HI2(N1:N2,N1:N2),HI3(N1:N2,N1:N2,N1:N2),
     ,          HI4(N1:N2,N1:N2,N1:N2,N1:N2) 
*
      z2 = zeta2
      z3 = zeta3
      DX = 1.D0/X
*
* ...Harmonic polylogs (HPLs) up to weight NW=3 by Gehrmann and Remiddi 
*
      CALL HPLOG (X, NW, HC1,HC2,HC3,HC4, HR1,HR2,HR3,HR4,
     ,            HI1,HI2,HI3,HI4, N1, N2)
*
* ...The coefficient function in terms of the HPLs
*
      cTeg2 =
     &  + cf*ca * (  - 36 - 106*x - 928.D0/27.D0*x**2 + 4438.D0/27.D0*
     &    dx + 64*z3 - 136*z3*x - 240*z3*dx + 32*z2 + 64*z2*x - 56*z2*
     &    dx + 16*Hr1(-1)*z2 + 8*Hr1(-1)*z2*x + 32*Hr1(-1)*z2*dx + 772.D
     &    0/3.D0*Hr1(0) + 172.D0/3.D0*Hr1(0)*x + 256.D0/9.D0*Hr1(0)*
     &    x**2 + 496.D0/3.D0*Hr1(0)*dx - 128*Hr1(0)*z2 - 16*Hr1(0)*z2*x
     &     + 32*Hr1(0)*z2*dx + 236.D0/3.D0*Hr1(1) + 4.D0/3.D0*Hr1(1)*x
     &     + 32.D0/3.D0*Hr1(1)*x**2 - 356.D0/3.D0*Hr1(1)*dx - 48*Hr1(1)
     &    *z2 + 24*Hr1(1)*z2*x + 32*Hr1(1)*z2*dx + 80*Hr2(-1,0) + 56*
     &    Hr2(-1,0)*x + 32.D0/3.D0*Hr2(-1,0)*x**2 + 80.D0/3.D0*Hr2(-1,0
     &    )*dx - 32*Hr2(0,0) + 4*Hr2(0,0)*x - 32.D0/3.D0*Hr2(0,0)*x**2
     &     - 464.D0/3.D0*Hr2(0,0)*dx - 96*Hr2(0,1) - 16*Hr2(0,1)*x - 32.
     &    D0/3.D0*Hr2(0,1)*x**2 + 496.D0/3.D0*Hr2(0,1)*dx - 64*Hr2(1,0)
     &     + 8*Hr2(1,0)*x + 64*Hr2(1,0)*dx + 96*Hr2(1,1) - 16*Hr2(1,1)*
     &    x + 32.D0/3.D0*Hr2(1,1)*x**2 - 344.D0/3.D0*Hr2(1,1)*dx - 32*
     &    Hr3(-1,-1,0) - 16*Hr3(-1,-1,0)*x + 96*Hr3(-1,0,0) + 48*Hr3(-1
     &    ,0,0)*x + 80*Hr3(-1,0,0)*dx - 32*Hr3(-1,0,1) - 16*Hr3(-1,0,1)
     &    *x - 32*Hr3(-1,0,1)*dx + 64*Hr3(0,-1,0) + 32*Hr3(0,-1,0)*x +
     &    64*Hr3(0,-1,0)*dx - 176*Hr3(0,0,0) - 248*Hr3(0,0,0)*x - 320*
     &    Hr3(0,0,0)*dx + 128*Hr3(0,0,1) + 96*Hr3(0,0,1)*x + 96*Hr3(0,0
     &    ,1)*dx + 96*Hr3(0,1,0) - 48*Hr3(0,1,0)*x - 96*Hr3(0,1,0)*dx
     &     - 48*Hr3(0,1,1) - 24*Hr3(0,1,1)*x - 16*Hr3(0,1,1)*dx + 64*
     &    Hr3(1,0,0) - 32*Hr3(1,0,0)*x - 48*Hr3(1,0,0)*dx - 32*Hr3(1,0,
     &    1) + 16*Hr3(1,0,1)*x + 32*Hr3(1,0,1)*dx - 64*Hr3(1,1,0) + 32*
     &    Hr3(1,1,0)*x + 64*Hr3(1,1,0)*dx + 16*Hr3(1,1,1) - 8*Hr3(1,1,1
     &    )*x - 16*Hr3(1,1,1)*dx )
      cTeg2 = cTeg2 + cf**2 * (  - 604.D0/5.D0 + 154.D0/5.D0*x - 16.D0/
     &    5.D0*x**2 + 316.D0/5.D0*dx + 32*z3 - 80*z3*x - 64*z3*dx - 32*
     &    z2 - 72*z2*x + 16.D0/5.D0*z2*x**3 + 64*Hr1(-1)*z2 + 32*Hr1(-1
     &    )*z2*x + 32*Hr1(-1)*z2*dx + 418.D0/5.D0*Hr1(0) - 262.D0/5.D0*
     &    Hr1(0)*x - 16.D0/5.D0*Hr1(0)*x**2 + 144.D0/5.D0*Hr1(0)*dx +
     &    32*Hr1(0)*z2 - 16*Hr1(0)*z2*x + 24*Hr1(1)*x - 8*Hr1(1)*dx +
     &    80*Hr1(1)*z2 - 40*Hr1(1)*z2*x - 48*Hr1(1)*z2*dx - 64*Hr2(-1,0
     &    ) - 96*Hr2(-1,0)*x + 16.D0/5.D0*Hr2(-1,0)*x**3 - 64.D0/5.D0*
     &    Hr2(-1,0)*dx**2 - 64*Hr2(0,0) + 166*Hr2(0,0)*x - 16.D0/5.D0*
     &    Hr2(0,0)*x**3 - 80*Hr2(0,1) + 12*Hr2(0,1)*x + 96*Hr2(0,1)*dx
     &     - 16*Hr2(1,0)*x + 112*Hr2(1,1) - 28*Hr2(1,1)*x - 96*Hr2(1,1)
     &    *dx + 128*Hr3(-1,-1,0) + 64*Hr3(-1,-1,0)*x + 64*Hr3(-1,-1,0)*
     &    dx - 64*Hr3(-1,0,0) - 32*Hr3(-1,0,0)*x - 32*Hr3(-1,0,0)*dx -
     &    128*Hr3(0,-1,0) + 88*Hr3(0,0,0) - 44*Hr3(0,0,0)*x + 16*Hr3(0,
     &    0,1) - 8*Hr3(0,0,1)*x - 64*Hr3(0,0,1)*dx + 64*Hr3(0,1,0) - 32
     &    *Hr3(0,1,0)*x - 64*Hr3(0,1,0)*dx - 80*Hr3(0,1,1) + 40*Hr3(0,1
     &    ,1)*x + 96*Hr3(0,1,1)*dx - 128*Hr3(1,0,0) + 64*Hr3(1,0,0)*x
     &     + 96*Hr3(1,0,0)*dx - 48*Hr3(1,0,1) + 24*Hr3(1,0,1)*x + 48*
     &    Hr3(1,0,1)*dx - 16*Hr3(1,1,0) + 8*Hr3(1,1,0)*x + 16*Hr3(1,1,0
     &    )*dx + 80*Hr3(1,1,1) - 40*Hr3(1,1,1)*x - 80*Hr3(1,1,1)*dx )
*
      C2G2TA = cTeg2 + CLG2TA (X, NF)
*
      RETURN
      END
*
************************************************************************
*     FL: quark non-singlet plus - regular term (A)
************************************************************************
      FUNCTION CLNSP2TA (X, NF)
      IMPLICIT REAL*8 (A-Z)
*
      INCLUDE "../commons/ColorFactors.h"
      INCLUDE "../commons/consts.h"
*
      COMPLEX*16 HC1, HC2, HC3, HC4, HC5
      INTEGER NF, NF2, N1, N2, NW
      PARAMETER ( N1 = -1, N2 = 1, NW = 3 )
      DIMENSION HC1(N1:N2),HC2(N1:N2,N1:N2),HC3(N1:N2,N1:N2,N1:N2),
     ,          HC4(N1:N2,N1:N2,N1:N2,N1:N2) 
      DIMENSION HR1(N1:N2),HR2(N1:N2,N1:N2),HR3(N1:N2,N1:N2,N1:N2),
     ,          HR4(N1:N2,N1:N2,N1:N2,N1:N2) 
      DIMENSION HI1(N1:N2),HI2(N1:N2,N1:N2),HI3(N1:N2,N1:N2,N1:N2),
     ,          HI4(N1:N2,N1:N2,N1:N2,N1:N2) 
*
      z2 = zeta2
      DX = 1.D0/X
*
* ...Harmonic polylogs (HPLs) up to weight NW=3 by Gehrmann and Remiddi 
*
      CALL HPLOG (X, NW, HC1,HC2,HC3,HC4, HR1,HR2,HR3,HR4,
     ,            HI1,HI2,HI3,HI4, N1, N2)
*
* ...The coefficient function in terms of the HPLs
*
      CLeq2 =
     &  + cf*ca * ( 1729.D0/45.D0 - 98.D0/15.D0*x - 16.D0/5.D0*x**2 -
     &    24.D0/5.D0*dx + 16.D0*z2*x + 16.D0/5.D0*z2*x**3 - 8.D0*Hr1(-1
     &    )*z2 - 146.D0/15.D0*Hr1(0) + 8.D0/5.D0*Hr1(0)*x - 16.D0/5.D0*
     &    Hr1(0)*x**2 + 24.D0/5.D0*Hr1(0)*dx + 46.D0/3.D0*Hr1(1) - 8.D0
     &    *Hr1(1)*z2 + 8.D0*Hr2(-1,0) + 16.D0*Hr2(-1,0)*x + 16.D0/5.D0*
     &    Hr2(-1,0)*x**3 - 24.D0/5.D0*Hr2(-1,0)*dx**2 - 16.D0*Hr2(0,0)*
     &    x - 16.D0/5.D0*Hr2(0,0)*x**3 - 16.D0*Hr3(-1,-1,0) + 8.D0*Hr3(
     &    -1,0,0) + 16.D0*Hr3(0,-1,0) + 8.D0*Hr3(1,0,0) )
      CLeq2 = CLeq2 + cf**2 * (  - 147.D0/5.D0 - 18.D0/5.D0*x + 32.D0/5.
     &    D0*x**2 + 48.D0/5.D0*dx + 4.D0*z2 - 32.D0*z2*x - 32.D0/5.D0*
     &    z2*x**3 + 16.D0*Hr1(-1)*z2 + 34.D0/5.D0*Hr1(0) + 24.D0/5.D0*
     &    Hr1(0)*x + 32.D0/5.D0*Hr1(0)*x**2 - 48.D0/5.D0*Hr1(0)*dx - 14.
     &    D0*Hr1(1) - 4.D0*Hr1(1)*x + 16.D0*Hr1(1)*z2 - 16.D0*Hr2(-1,0)
     &     - 32.D0*Hr2(-1,0)*x - 32.D0/5.D0*Hr2(-1,0)*x**3 + 48.D0/5.D0
     &    *Hr2(-1,0)*dx**2 - 12.D0*Hr2(0,0) + 32.D0*Hr2(0,0)*x + 32.D0/
     &    5.D0*Hr2(0,0)*x**3 - 4.D0*Hr2(0,1) - 16.D0*Hr2(1,0) + 8.D0*
     &    Hr2(1,1) + 32.D0*Hr3(-1,-1,0) - 16.D0*Hr3(-1,0,0) - 32.D0*
     &    Hr3(0,-1,0) - 16.D0*Hr3(1,0,0) )
      CLeq2 = CLeq2 + nf*cf * (  - 50.D0/9.D0 + 4.D0/3.D0*x + 4.D0/3.D0
     &    *Hr1(0) - 4.D0/3.D0*Hr1(1) )
*
      CLNSP2TA = CLeq2
*
      RETURN
      END
*
************************************************************************
*     FL: quark pure-singlet - regular term (A)
************************************************************************
      FUNCTION CLPS2TA (X, NF)
      IMPLICIT REAL*8 (A-Z)
*
      INCLUDE "../commons/ColorFactors.h"
      INCLUDE "../commons/consts.h"
*
      COMPLEX*16 HC1, HC2, HC3, HC4, HC5
      INTEGER NF, NF2, N1, N2, NW
      PARAMETER ( N1 = -1, N2 = 1, NW = 2 )
      DIMENSION HC1(N1:N2),HC2(N1:N2,N1:N2),HC3(N1:N2,N1:N2,N1:N2),
     ,          HC4(N1:N2,N1:N2,N1:N2,N1:N2) 
      DIMENSION HR1(N1:N2),HR2(N1:N2,N1:N2),HR3(N1:N2,N1:N2,N1:N2),
     ,          HR4(N1:N2,N1:N2,N1:N2,N1:N2) 
      DIMENSION HI1(N1:N2),HI2(N1:N2,N1:N2),HI3(N1:N2,N1:N2,N1:N2),
     ,          HI4(N1:N2,N1:N2,N1:N2,N1:N2) 
*
      z2 = zeta2
      DX = 1.D0/X
*
* ...Harmonic polylogs (HPLs) up to weight NW=2 by Gehrmann and Remiddi 
*
      CALL HPLOG (X, NW, HC1,HC2,HC3,HC4, HR1,HR2,HR3,HR4,
     ,            HI1,HI2,HI3,HI4, N1, N2)
*
* ...The coefficient function in terms of the HPLs
*
      cLeqps2 =
     &  + nf*cf * (  - 56.D0/3.D0 + 104.D0/3.D0*x - 8*x**2 - 8*dx + 8*
     &    z2 - 16*Hr1(0) - 16*Hr1(0)*x + 8.D0/3.D0*Hr1(0)*x**2 + 32.D0/
     &    3.D0*Hr1(0)*dx + 8*Hr1(1)*x - 8.D0/3.D0*Hr1(1)*x**2 - 16.D0/3.
     &    D0*Hr1(1)*dx + 24*Hr2(0,0) - 8*Hr2(0,1) )
*
      CLPS2TA = cLeqps2
*
      RETURN
      END
*
************************************************************************
*     FL: gluon - regular term (A)
************************************************************************
      FUNCTION CLG2TA (X, NF)
      IMPLICIT REAL*8 (A-Z)
*
      INCLUDE "../commons/ColorFactors.h"
      INCLUDE "../commons/consts.h"
*
      COMPLEX*16 HC1, HC2, HC3, HC4, HC5
      INTEGER NF, NF2, N1, N2, NW
      PARAMETER ( N1 = -1, N2 = 1, NW = 2 )
      DIMENSION HC1(N1:N2),HC2(N1:N2,N1:N2),HC3(N1:N2,N1:N2,N1:N2),
     ,          HC4(N1:N2,N1:N2,N1:N2,N1:N2) 
      DIMENSION HR1(N1:N2),HR2(N1:N2,N1:N2),HR3(N1:N2,N1:N2,N1:N2),
     ,          HR4(N1:N2,N1:N2,N1:N2,N1:N2) 
      DIMENSION HI1(N1:N2),HI2(N1:N2,N1:N2),HI3(N1:N2,N1:N2,N1:N2),
     ,          HI4(N1:N2,N1:N2,N1:N2,N1:N2) 
*
      z2 = zeta2
      DX = 1.D0/X
*
* ...Harmonic polylogs (HPLs) up to weight NW=2 by Gehrmann and Remiddi 
*
      CALL HPLOG (X, NW, HC1,HC2,HC3,HC4, HR1,HR2,HR3,HR4,
     ,            HI1,HI2,HI3,HI4, N1, N2)
*
* ...The coefficient function in terms of the HPLs
*
      cLeg2 =
     &  + cf*ca * (  - 320.D0/3.D0 - 160.D0/3.D0*x + 32.D0/3.D0*x**2 +
     &    448.D0/3.D0*dx - 64*z2 + 32*z2*dx + 112*Hr1(0) + 32*Hr1(0)*x
     &     - 16.D0/3.D0*Hr1(0)*x**2 - 352.D0/3.D0*Hr1(0)*dx - 144*Hr1(1
     &    ) - 16*Hr1(1)*x + 16.D0/3.D0*Hr1(1)*x**2 + 464.D0/3.D0*Hr1(1)
     &    *dx + 32*Hr2(-1,0) + 32*Hr2(-1,0)*dx - 96*Hr2(0,0) - 128*Hr2(
     &    0,0)*dx + 64*Hr2(0,1) + 64*Hr2(1,0) - 64*Hr2(1,0)*dx - 32*
     &    Hr2(1,1) + 32*Hr2(1,1)*dx )
      cLeg2 = cLeg2 + cf**2 * ( 24.D0/5.D0 + 248.D0/15.D0*x - 32.D0/15.D
     &    0*x**2 - 96.D0/5.D0*dx + 16*z2 + 32.D0/15.D0*z2*x**3 - 8.D0/5.
     &    D0*Hr1(0) - 224.D0/15.D0*Hr1(0)*x - 32.D0/15.D0*Hr1(0)*x**2
     &     + 96.D0/5.D0*Hr1(0)*dx + 24*Hr1(1) + 8*Hr1(1)*x - 32*Hr1(1)*
     &    dx - 32.D0/3.D0*Hr2(-1,0) + 32.D0/15.D0*Hr2(-1,0)*x**3 + 64.D0
     &    /5.D0*Hr2(-1,0)*dx**2 + 48*Hr2(0,0) - 32.D0/15.D0*Hr2(0,0)*
     &    x**3 - 16*Hr2(0,1) )
*
      CLG2TA = cLeg2
*
      RETURN
      END
*
************************************************************************
*     F3: quark non-singlet plus - regular term (A)
************************************************************************
      FUNCTION C3NSP2TA (X, NF)
      IMPLICIT REAL*8 (A - Z)
*
      INCLUDE "../commons/ColorFactors.h"
      INCLUDE "../commons/consts.h"
*
      COMPLEX*16 HC1, HC2, HC3, HC4 
      INTEGER NF, NF2, N1, N2, NW, I1, I2, I3, N
      PARAMETER ( N1 = -1, N2 = 1, NW = 3 ) 
      DIMENSION HC1(N1:N2),HC2(N1:N2,N1:N2),HC3(N1:N2,N1:N2,N1:N2), 
     ,          HC4(N1:N2,N1:N2,N1:N2,N1:N2) 
      DIMENSION HR1(N1:N2),HR2(N1:N2,N1:N2),HR3(N1:N2,N1:N2,N1:N2), 
     ,          HR4(N1:N2,N1:N2,N1:N2,N1:N2) 
      DIMENSION HI1(N1:N2),HI2(N1:N2,N1:N2),HI3(N1:N2,N1:N2,N1:N2), 
     ,          HI4(N1:N2,N1:N2,N1:N2,N1:N2) 
*
* ...Some abbreviations
*
      z2 = zeta2
      z3 = zeta3
      DX = 1.D0/X
      DM = 1.D0/(1.D0-X)
      DP = 1.D0/(1.D0+X)
      DL1 = LOG (1.D0-X)
*
* ...The harmonic polylogs up to weight 4 by Gehrmann and Remiddi
*
      CALL HPLOG (X, NW, HC1,HC2,HC3,HC4, HR1,HR2,HR3,HR4,
     ,            HI1,HI2,HI3,HI4, N1, N2) 
*
* ...The coefficient function in terms of the harmonic polylogs
*    (without the delta(1-x) part, but with the soft contribution)
*
      CAeq2 =
     &  + cf*ca * ( 325.D0/54.D0 + 895.D0/54.D0*x - 3155.D0/54.D0*dm - 
     &    36.D0*z3 + 28.D0*z3*dp + 28.D0*z3*dm + 12.D0*z2 + 8.D0*z2*x
     &     + 8.D0*z2*x**2 + 20.D0*Hr1(-1)*z2 - 12.D0*Hr1(-1)*z2*x - 32.D
     &    0*Hr1(-1)*z2*dp + 27.D0*Hr1(0) - 193.D0/3.D0*Hr1(0)*x - 8.D0*
     &    Hr1(0)*dp + 206.D0/3.D0*Hr1(0)*dm + 20.D0*Hr1(0)*z2 + 12.D0*
     &    Hr1(0)*z2*x - 8.D0*Hr1(0)*z2*dp - 32.D0*Hr1(0)*z2*dm - 19.D0/
     &    9.D0*Hr1(1) + 305.D0/9.D0*Hr1(1)*x - 367.D0/9.D0*Hr1(1)*dm + 
     &    8.D0*Hr1(1)*z2*x - 8.D0*Hr1(1)*z2*dm + 4.D0*Hr2(-1,0) + 4.D0*
     &    Hr2(-1,0)*x + 8.D0*Hr2(-1,0)*x**2 + 8.D0*Hr2(-1,0)*dx + 59.D0/
     &    3.D0*Hr2(0,0) + 71.D0/3.D0*Hr2(0,0)*x - 8.D0*Hr2(0,0)*x**2 - 
     &    22.D0/3.D0*Hr2(0,0)*dm - 46.D0/3.D0*Hr2(0,1) - 46.D0/3.D0*
     &    Hr2(0,1)*x + 44.D0/3.D0*Hr2(0,1)*dm + 22.D0/3.D0*Hr2(1,1) + 
     &    22.D0/3.D0*Hr2(1,1)*x - 44.D0/3.D0*Hr2(1,1)*dm + 24.D0*Hr3(-1
     &    ,-1,0) - 8.D0*Hr3(-1,-1,0)*x - 32.D0*Hr3(-1,-1,0)*dp + 8.D0*
     &    Hr3(-1,0,0) - 16.D0*Hr3(-1,0,0)*x - 24.D0*Hr3(-1,0,0)*dp - 8.D
     &    0*Hr3(-1,0,1) )
      CAeq2 = CAeq2 + cf*ca * ( 8.D0*Hr3(-1,0,1)*x + 16.D0*Hr3(-1,0,1)*
     &    dp + 16.D0*Hr3(0,-1,0) - 8.D0*Hr3(0,-1,0)*dp - 24.D0*Hr3(0,-1
     &    ,0)*dm - 36.D0*Hr3(0,0,0) + 36.D0*Hr3(0,0,0)*dp + 36.D0*Hr3(0
     &    ,0,0)*dm + 4.D0*Hr3(0,0,1) - 4.D0*Hr3(0,0,1)*x - 8.D0*Hr3(0,0
     &    ,1)*dp + 4.D0*Hr3(0,1,0) + 4.D0*Hr3(0,1,0)*x - 8.D0*Hr3(0,1,0
     &    )*dm + 12.D0*Hr3(1,0,0) + 4.D0*Hr3(1,0,0)*x - 16.D0*Hr3(1,0,0
     &    )*dm - 4.D0*Hr3(1,0,1) - 4.D0*Hr3(1,0,1)*x + 8.D0*Hr3(1,0,1)*
     &    dm + 4.D0*Hr3(1,1,0) + 4.D0*Hr3(1,1,0)*x - 8.D0*Hr3(1,1,0)*dm
     &     )
      CAeq2 = CAeq2 + cf**2 * (  - 19.D0/2.D0 + 19.D0/2.D0*x + 51.D0/2.D
     &    0*dm + 128.D0*z3 + 56.D0*z3*x - 56.D0*z3*dp - 152.D0*z3*dm - 
     &    52.D0*z2 - 20.D0*z2*x - 16.D0*z2*x**2 + 12.D0*z2*dm - 40.D0*
     &    Hr1(-1)*z2 + 24.D0*Hr1(-1)*z2*x + 64.D0*Hr1(-1)*z2*dp - 2.D0*
     &    Hr1(0) + 92.D0*Hr1(0)*x + 16.D0*Hr1(0)*dp - 106.D0*Hr1(0)*dm
     &     - 20.D0*Hr1(0)*z2 - 4.D0*Hr1(0)*z2*x + 16.D0*Hr1(0)*z2*dp + 
     &    40.D0*Hr1(0)*z2*dm - 5.D0*Hr1(1) - 33.D0*Hr1(1)*x + 27.D0*
     &    Hr1(1)*dm + 8.D0*Hr1(1)*z2 - 8.D0*Hr1(1)*z2*x - 8.D0*Hr2(-1,0
     &    ) - 8.D0*Hr2(-1,0)*x - 16.D0*Hr2(-1,0)*x**2 - 16.D0*Hr2(-1,0)
     &    *dx - 86.D0*Hr2(0,0) - 74.D0*Hr2(0,0)*x + 16.D0*Hr2(0,0)*x**2
     &     + 66.D0*Hr2(0,0)*dm + 32.D0*Hr2(0,1) + 8.D0*Hr2(0,1)*x + 12.D
     &    0*Hr2(0,1)*dm - 12.D0*Hr2(1,0) + 12.D0*Hr2(1,0)*x + 24.D0*
     &    Hr2(1,0)*dm + 8.D0*Hr2(1,1) + 16.D0*Hr2(1,1)*x - 36.D0*Hr2(1,
     &    1)*dm - 48.D0*Hr3(-1,-1,0) + 16.D0*Hr3(-1,-1,0)*x + 64.D0*
     &    Hr3(-1,-1,0)*dp - 16.D0*Hr3(-1,0,0) + 32.D0*Hr3(-1,0,0)*x + 
     &    48.D0*Hr3(-1,0,0)*dp )
      CAeq2 = CAeq2 + cf**2 * ( 16.D0*Hr3(-1,0,1) - 16.D0*Hr3(-1,0,1)*x
     &     - 32.D0*Hr3(-1,0,1)*dp - 32.D0*Hr3(0,-1,0) + 16.D0*Hr3(0,-1,
     &    0)*dp + 48.D0*Hr3(0,-1,0)*dm + 138.D0*Hr3(0,0,0) + 66.D0*Hr3(
     &    0,0,0)*x - 72.D0*Hr3(0,0,0)*dp - 160.D0*Hr3(0,0,0)*dm - 24.D0
     &    *Hr3(0,0,1) - 8.D0*Hr3(0,0,1)*x + 16.D0*Hr3(0,0,1)*dp + 8.D0*
     &    Hr3(0,0,1)*dm + 36.D0*Hr3(0,1,0) + 36.D0*Hr3(0,1,0)*x - 72.D0
     &    *Hr3(0,1,0)*dm - 16.D0*Hr3(0,1,1) - 16.D0*Hr3(0,1,1)*x + 40.D0
     &    *Hr3(0,1,1)*dm - 28.D0*Hr3(1,0,0) - 12.D0*Hr3(1,0,0)*x + 40.D0
     &    *Hr3(1,0,0)*dm - 16.D0*Hr3(1,0,1) - 16.D0*Hr3(1,0,1)*x + 32.D0
     &    *Hr3(1,0,1)*dm - 24.D0*Hr3(1,1,0) - 24.D0*Hr3(1,1,0)*x + 48.D0
     &    *Hr3(1,1,0)*dm + 24.D0*Hr3(1,1,1) + 24.D0*Hr3(1,1,1)*x - 48.D0
     &    *Hr3(1,1,1)*dm )
      CAeq2 = CAeq2 + nf*cf * ( 55.D0/27.D0 - 131.D0/27.D0*x + 247.D0/
     &    27.D0*dm + 2.D0*Hr1(0) + 22.D0/3.D0*Hr1(0)*x - 32.D0/3.D0*
     &    Hr1(0)*dm - 2.D0/9.D0*Hr1(1) - 38.D0/9.D0*Hr1(1)*x + 58.D0/9.D
     &    0*Hr1(1)*dm - 2.D0/3.D0*Hr2(0,0) - 2.D0/3.D0*Hr2(0,0)*x + 4.D0
     &    /3.D0*Hr2(0,0)*dm + 4.D0/3.D0*Hr2(0,1) + 4.D0/3.D0*Hr2(0,1)*x
     &     - 8.D0/3.D0*Hr2(0,1)*dm - 4.D0/3.D0*Hr2(1,1) - 4.D0/3.D0*
     &    Hr2(1,1)*x + 8.D0/3.D0*Hr2(1,1)*dm )
*
* ...The soft (`+'-distribution) part of the coefficient function
*
      A3 = 
     &     + 8.D0*cf**2
      A2 =
     &     - 22.D0/3.D0*ca*cf
     &     - 18.D0*cf**2
     &     + 4.D0/3.D0*cf*nf
      A1 =
     &     - 8.D0*z2*ca*cf
     &     + 16.D0*z2*cf**2
     &     + 367.D0/9.D0*ca*cf
     &     - 27.D0*cf**2
     &     - 58.D0/9.D0*cf*nf
      A0 =
     &     + 44.D0/3.D0*z2*ca*cf
     &     + 40.D0*z3*ca*cf
     &     - 8.D0*z3*cf**2
     &     - 3155.D0/54.D0*ca*cf
     &     + 51.D0/2.D0*cf**2
     &     + 247.D0/27.D0*cf*nf
     &     - 8.D0/3.D0*z2*cf*nf
*
      CAeq2L = DM * ( DL1**3 * A3 + DL1**2 * A2 + DL1 * A1 + A0 ) 
*
* ...The regular piece of the coefficient function
*
      C3NSP2TA = CAeq2 - CAeq2L
*
      RETURN
      END
*
************************************************************************
*     F3: quark non-singlet - singular term (B)
************************************************************************
      FUNCTION C3NS2TB (X, NF)
      IMPLICIT REAL*8 (A - Z)
*
      INCLUDE "../commons/ColorFactors.h"
      INCLUDE "../commons/consts.h"
*
      INTEGER NF
*
      z2 = zeta2
      z3 = zeta3
*
      A3 = 
     &     + 8.D0*cf**2
      A2 =
     &     - 22.D0/3.D0*ca*cf
     &     - 18.D0*cf**2
     &     + 4.D0/3.D0*cf*nf
      A1 =
     &     - 8.D0*z2*ca*cf
     &     + 16.D0*z2*cf**2
     &     + 367.D0/9.D0*ca*cf
     &     - 27.D0*cf**2
     &     - 58.D0/9.D0*cf*nf
      A0 =
     &     + 44.D0/3.D0*z2*ca*cf
     &     + 40.D0*z3*ca*cf
     &     - 8.D0*z3*cf**2
     &     - 3155.D0/54.D0*ca*cf
     &     + 51.D0/2.D0*cf**2
     &     + 247.D0/27.D0*cf*nf
     &     - 8.D0/3.D0*z2*cf*nf
*
       DL1 = LOG (1.D0-X)
       DM  = 1.D0/(1.D0-X)
*
       C3NS2TB = DM * ( DL1**3 * A3 + DL1**2 * A2 + DL1 * A1 + A0)
*
       RETURN
       END
*
************************************************************************
*     F3: quark non-singlet plus - local term (C)
************************************************************************
      FUNCTION C3NSP2TC (X, NF)
      IMPLICIT REAL*8 (A - Z)
*
      INCLUDE "../commons/ColorFactors.h"
      INCLUDE "../commons/consts.h"
*
      INTEGER NF
*
      z2 = zeta2
      z3 = zeta3
*
      A3 = 
     &     + 8.D0*cf**2
      A2 =
     &     - 22.D0/3.D0*ca*cf
     &     - 18.D0*cf**2
     &     + 4.D0/3.D0*cf*nf
      A1 =
     &     - 8.D0*z2*ca*cf
     &     + 16.D0*z2*cf**2
     &     + 367.D0/9.D0*ca*cf
     &     - 27.D0*cf**2
     &     - 58.D0/9.D0*cf*nf
      A0 =
     &     + 44.D0/3.D0*z2*ca*cf
     &     + 40.D0*z3*ca*cf
     &     - 8.D0*z3*cf**2
     &     - 3155.D0/54.D0*ca*cf
     &     + 51.D0/2.D0*cf**2
     &     + 247.D0/27.D0*cf*nf
     &     - 8.D0/3.D0*z2*cf*nf
*
* ...The coefficient of delta(1-x)
*
      C2DELT = 
     ,    + ca*cf * ( - 5465.D0/72.D0 + 140.D0/3.D0*z3 + 215.D0/3.D0*z2
     ,                - 49.D0/5.D0*z2**2 )
     ,    + cf**2 * ( 331.D0/8.D0 - 78.D0*z3 - 39.D0*z2 + 30.D0*z2**2 )
     ,    + cf*nf * ( 457.D0/36.D0 + 4.D0/3.D0*z3 - 38.D0/3.D0*z2 )
*
      DL1 = LOG (1.D0-X)
*
      C3NSP2TC =   DL1**4 * A3/4.D0 + DL1**3 * A2/3.D0 
     ,           + DL1**2 * A1/2.D0 + DL1 * A0 + C2DELT
*
       RETURN
       END
*
c$$$************************************************************************
c$$$*
c$$$*     Analytical terms neededed for the scale variations
c$$$*
c$$$************************************************************************
c$$$      FUNCTION P0P0NSA(X)
c$$$*
c$$$      IMPLICIT NONE
c$$$*
c$$$      INCLUDE "../commons/ColorFactors.h"
c$$$**
c$$$*     Input Variables
c$$$*
c$$$      DOUBLE PRECISION X
c$$$**
c$$$*     Output Variables
c$$$*
c$$$      DOUBLE PRECISION P0P0NSA
c$$$*
c$$$      P0P0NSA = 4D0 * CF**2D0 * ( - 4D0 * DLOG(X) / ( 1D0 - X )
c$$$     1     - 4D0 * ( 1D0 + X ) * DLOG( 1D0 - X )
c$$$     2     + 3D0 * ( 1D0 + X ) * DLOG(X) - ( X + 5D0 ) )
c$$$*
c$$$      RETURN
c$$$      END
c$$$*
c$$$************************************************************************
c$$$      FUNCTION P0P0NSB(X)
c$$$*
c$$$      IMPLICIT NONE
c$$$*
c$$$      INCLUDE "../commons/ColorFactors.h"
c$$$**
c$$$*     Input Variables
c$$$*
c$$$      DOUBLE PRECISION X
c$$$**
c$$$*     Output Variables
c$$$*
c$$$      DOUBLE PRECISION P0P0NSB
c$$$*
c$$$      P0P0NSB = 4D0 * CF**2D0 * ( 8D0 * DLOG( 1D0 - X ) + 6D0 ) 
c$$$     1     / ( 1D0 - X )
c$$$*
c$$$      RETURN
c$$$      END
c$$$*
c$$$************************************************************************
c$$$      FUNCTION P0P0NSC(X)
c$$$*
c$$$      IMPLICIT NONE
c$$$*
c$$$      INCLUDE "../commons/ColorFactors.h"
c$$$      INCLUDE "../commons/consts.h"
c$$$**
c$$$*     Input Variables
c$$$*
c$$$      DOUBLE PRECISION X
c$$$**
c$$$*     Output Variables
c$$$*
c$$$      DOUBLE PRECISION P0P0NSC
c$$$*
c$$$      P0P0NSC = 4D0 * CF**2D0 * ( 4D0 * DLOG( 1D0 - X )**2D0 
c$$$     1     + 6D0 * DLOG( 1D0 - X ) + ( 9D0 / 4D0 - 4D0 * ZETA2 ) )
c$$$*
c$$$      RETURN
c$$$      END
c$$$*
c$$$************************************************************************
c$$$      FUNCTION P0P0PSA(X,NF)
c$$$*
c$$$      IMPLICIT NONE
c$$$*
c$$$      INCLUDE "../commons/ColorFactors.h"
c$$$**
c$$$*     Input Variables
c$$$*
c$$$      INTEGER NF
c$$$      DOUBLE PRECISION X
c$$$**
c$$$*     Output Variables
c$$$*
c$$$      DOUBLE PRECISION P0P0PSA
c$$$*
c$$$      P0P0PSA = NF * CF * TR * ( 8D0 * ( 3D0 + 4D0 / X - 3D0 * X 
c$$$     1     - 4D0 * X**2D0 ) / 3D0 + 16D0 * ( 1D0 + X ) * DLOG(X) ) 
c$$$*
c$$$      RETURN
c$$$      END
c$$$*
c$$$************************************************************************
c$$$      FUNCTION P0P0GA(X,NF)
c$$$*
c$$$      IMPLICIT NONE
c$$$*
c$$$      INCLUDE "../commons/ColorFactors.h"
c$$$**
c$$$*     Input Variables
c$$$*
c$$$      INTEGER NF
c$$$      DOUBLE PRECISION X
c$$$**
c$$$*     Output Variables
c$$$*
c$$$      DOUBLE PRECISION P0P0GA
c$$$*
c$$$      P0P0GA = CA * TR * NF * ( 16D0 * ( 2D0 * X**2D0 - 2D0 * X + 1D0 )
c$$$     1     * DLOG( 1D0 - X ) + 16D0 * ( 4D0 * X + 1 D0 ) * DLOG(X)
c$$$     2     + 4D0 * ( - 40D0 * X**2D0 + 26D0 * X + 17D0 + 8D0 / X )
c$$$     3     / 3D0 )
c$$$     4     + CF * TR * NF * ( ( 2D0 * X**2D0 - 2D0 * X + 1D0 )
c$$$     5     * DLOG( 1D0 - X ) - 8D0 * ( 4D0 *X**2D0 - 2D0 * X + 1D0 )
c$$$     6     * DLOG(X) + 4D0 * ( 4D0 * X - 1D0 ) )
c$$$     7     + NF**2D0 * TR**2D0 * ( - 16D0 * ( 2D0 * X**2D0 - 2D0 * X
c$$$     8     + 1D0 ) )
c$$$*
c$$$      RETURN
c$$$      END
c$$$*
c$$$************************************************************************
c$$$      FUNCTION CL1P0NSA(X)
c$$$*
c$$$      IMPLICIT NONE
c$$$*
c$$$      INCLUDE "../commons/ColorFactors.h"
c$$$**
c$$$*     Input Variables
c$$$*
c$$$      DOUBLE PRECISION X
c$$$**
c$$$*     Output Variables
c$$$*
c$$$      DOUBLE PRECISION CL1P0NSA
c$$$*
c$$$      CL1P0NSA = 4D0 * CF**2D0 * ( ( X + 2D0 ) + 4D0 * X
c$$$     1     * DLOG( 1D0 - X ) - 2D0 * X * DLOG(X) )
c$$$*
c$$$      RETURN
c$$$      END
c$$$*
c$$$************************************************************************
c$$$      FUNCTION CL1P0PSA(X)
c$$$*
c$$$      IMPLICIT NONE
c$$$*
c$$$      INCLUDE "../commons/ColorFactors.h"
c$$$**
c$$$*     Input Variables
c$$$*
c$$$      DOUBLE PRECISION X
c$$$**
c$$$*     Output Variables
c$$$*
c$$$      DOUBLE PRECISION CL1P0PSA
c$$$*
c$$$      CL1P0PSA = CF * TR * ( 32D0 * ( 1D0 / X - 3D0 + 2D0 * X**2D0 ) 
c$$$     1     / 3D0 - 32D0 * X * DLOG(X) )
c$$$*
c$$$      RETURN
c$$$      END
c$$$*
c$$$************************************************************************
c$$$      FUNCTION CL1P0GA(X,NF)
c$$$*
c$$$      IMPLICIT NONE
c$$$*
c$$$      INCLUDE "../commons/ColorFactors.h"
c$$$**
c$$$*     Input Variables
c$$$*
c$$$      INTEGER NF
c$$$      DOUBLE PRECISION X
c$$$**
c$$$*     Output Variables
c$$$*
c$$$      DOUBLE PRECISION CL1P0GA
c$$$*
c$$$      CL1P0GA = CA * TR * ( 64D0 * X * ( 1D0 - X ) * DLOG( 1D0 - X )
c$$$     1     - 128D0 * X * DLOG(X) + 16D0 * ( 23D0 * X**2D0 - 19D0 * X
c$$$     2     - 6D0 + 2D0 / X ) / 3D0 )
c$$$     3     + CF * TR * ( 16D0 * X * DLOG(X) / 3D0 - 8D0 * ( 2D0 * X**2D0 
c$$$     4     - X - 1D0 ) / 3D0 )
c$$$     5     + NF * TR**2D0 * ( - 64D0 * X * ( 1D0 - X ) / 3D0 )
c$$$*
c$$$      RETURN
c$$$      END
