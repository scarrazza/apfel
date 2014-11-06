************************************************************************
*
*     ZMCoefficientFunctions.f:
*
*     This file contains all the ZM coefficient functions as parametrized
*     in the papers written below.
*
************************************************************************
*
* ..File: xc2sg2p.f    F2_S
*
*
* ..Calculation of the 2-loop x-space MS(bar) coefficient functions 
*    for F2 via compact parametrizations involving only logarithms.
*    Singlet, mu_r = mu_f = Q. Expansion parameter: alpha_s/(4 pi).
*
*  ..The distributions (in the mathematical sense) are given as in eq.
*    (B.26) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.
*    The name-endings A, B, and C of the functions below correspond to 
*    the kernel superscripts [2], [3], and [1] in that equation.
*
*  ..The relative accuracy of the coefficient functions, as well as of 
*    the convolution results, amounts to a few thousandth.
*
*  ..Reference: W.L. van Neerven and A. Vogt, 
*               hep-ph/0006154 = Nucl. Phys. B588 (2000) 345
*  ..The user should also cite the original calculations,
*     E.B. Zijlstra and W.L. van Neerven, Phys. Lett. B272 (1991) 127
*     (pure singlet) and Phys. Lett. B273 (1991) 476 (gluon).
*
* 
* =====================================================================
*
*
* ..This is the pure singlet piece, denoted by C2S in WvN's program. 
*    Seven numerical coefficients (all but the one of 1/y, which is 
*    exact up to truncation) are fitted to his results, using x values
*    between 10^-6 and 1-10^-6.
*
       FUNCTION C2S2A (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       C2S2A =   NF * ( 5.290 * (1./Y-1.) + 4.310 * DL**3   
     1         - 2.086 * DL**2 + 39.78 * DL - 0.101 * (1.-Y) * DL1**3 
     2         - (24.75 - 13.80 * Y) * DL**2 * DL1 + 30.23 * DL * DL1 )
C       C2S2A = NF * ( ( 8D0 / 3D0 * DL1**2D0 - 32D0 / 3D0 * DL1 
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
* ---------------------------------------------------------------------
*
*
* ..This is the regular gluon piece, denoted by C2G2 in WvN's program. 
*    Nine numerical coefficients are fitted as above, the ones of 1/y, 
*    ln^3(1-y), and ln^2(1-y) are exact up to truncation.
*
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
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' gluon piece, which has no counterpart in WvN's
*    program, as it does not exist in the exact expressions. Here it 
*    is, however, relevant for achieving the highest accuracy of the 
*    convolution, as are the adjustments of the constant in the non-
*    singlet quark coefficient functions. The value is fixed from the 
*    lowest moments. 
*
       FUNCTION C2G2C (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       C2G2C = - NF * 0.28  
*
       RETURN
       END
*
* =================================================================av==
*
* ..File: xc2ns2p.f    F2_NS
*
*
* ..Calculation of the 2-loop x-space MS(bar) coefficient functions 
*    for F2 via compact parametrizations involving only logarithms.
*    Non-singlet, mu_r = mu_f = Q. Expansion parameter: alpha_s/(4 pi).
*
*  ..The distributions (in the mathematical sense) are given as in eq.
*    (B.26) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.
*    The name-endings A, B, and C of the functions below correspond to 
*    the kernel superscripts [2], [3], and [1] in that equation.
*
*  ..The relative accuracy of the coefficient functions, as well as of 
*    the convolution results, amounts to a few thousandth.
*    
*  ..Reference: W.L. van Neerven and A. Vogt, 
*               hep-ph/9907472 = Nucl. Phys. B568 (2000) 263
*  ..The user should also cite the original calculation, 
*     E.B. Zijlstra and W.L. van Neerven, Phys. Lett. B272 (1991) 127.
*
* 
* =====================================================================
*
*
* ..This is the regular non-singlet piece for the electromagnetic F2, 
*    corresponding to C2NSP+C2NSN in W. van Neerven's program. The 
*    (10+8) numerical coefficients are fitted to his results, using x 
*    values between 10^-6 and 1-10^-6. 
*
       FUNCTION C2NN2A (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       C2NN2A = 
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
* =====================================================================
*
* ..Derivative with respect to NF of C2NN2A
*
       FUNCTION C2NN2A_DNF (Y)
       IMPLICIT REAL*8 (A-Z)
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       C2NN2A_DNF = 
     1          - 5.691 - 37.91 * Y 
     2          + 2.244 * DL**2 + 5.770 * DL 
     3          - 1.707 * DL1**2  + 22.95 * DL1
     4          + 3.036 * DL**2 * DL1 + 17.97 * DL * DL1
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the regular non-singlet piece for the odd-moment (CC) F2, 
*    corresponding to C2NSP-C2NSN in WvN's program. For the NF^0 piece
*    8 numerical coefficients are fitted to his results, the ones of
*    ln^3(1-y) and ln^2(1-y) are taken over from C3NN2A. The NF piece
*    is also the same as in C3NN2A.
*
       FUNCTION C2NC2A (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       C2NC2A = 
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
* ---------------------------------------------------------------------
*
*
* ..This is the singular NS piece, denoted by SOFT2 in WvN's program. 
*    It is the same for all F2 and F3 cases. The numerical coefficients 
*    are exact, but truncated.
*
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
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' NS piece for the e.m. F2, denoted by COR2 in 
*    WvN's program. The numerical coefficients of the logs are exact,
*    but truncated, the constant one (from the delta-function) is 
*    slightly adjusted (+ 0.485 - 0.0035 NF) using the lowest moments.
*
       FUNCTION C2NN2C (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL1 = LOG (1.-Y)
*
       C2NN2C = 
     1          + 3.55555 * DL1**4 - 20.4444 * DL1**3 - 15.5525 * DL1**2
     2          + 188.64 * DL1 - 338.531 + 0.485 
     3        + NF * (0.592593 * DL1**3 - 4.2963 * DL1**2 
     4          + 6.3489 * DL1 + 46.844 - 0.0035)
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' NS piece for the CC F2, also given by COR2 in 
*    WvN's program. The numerical coefficients of the logs are exact,
*    but truncated, the constant one is adjusted (- 0.2652 - 0.0035 NF) 
*    using the lowest moments.
*
       FUNCTION C2NC2C (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL1 = LOG (1.-Y)
*
       C2NC2C = 
     1          + 3.55555 * DL1**4 - 20.4444 * DL1**3 - 15.5525 * DL1**2
     2          + 188.64 * DL1 - 338.531 + 0.537
     3        + NF * (0.592593 * DL1**3 - 4.2963 * DL1**2 
     4          + 6.3489 * DL1 + 46.844 - 0.0035)
*
       RETURN
       END
*
* =================================================================av==
*
* ..File: xc3ns2p.f    F3_NS
*
*
* ..Calculation of the 2-loop x-space MS(bar) coefficient functions 
*    for xF3 via compact parametrizations involving only logarithms.
*    mu_r = mu_f = Q. Expansion parameter: alpha_s/(4 pi).
*
*  ..The distributions (in the mathematical sense) are given as in eq.
*    (B.26) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.
*    The name-endings A, B, and C of the functions below correspond to 
*    the kernel superscripts [2], [3], and [1] in that equation.
*
*  ..The relative accuracy of the coefficient functions, as well as of 
*    the convolution results, amounts to a few thousandth.
*    
*  ..Reference: W.L. van Neerven and A. Vogt,
*               hep-ph/9907472 = Nucl. Phys. B568 (2000) 263
*  ..The user should also cite the original calculation,
*     E.B. Zijlstra and W.L. van Neerven, Phys. Lett. B297 (1992) 377.
*
* 
* =====================================================================
*
*                                                               __
* ..This is the regular non-singlet piece for the sum F3(nu)+F3(nu), 
*    corresponding to C3NSP-C3NSN in W. van Neerven's program. The 
*    (9+8) numerical coefficients are fitted to his results, using x 
*    values between 10^-6 and 1-10^-6. 
*
       FUNCTION C3NM2A (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       C3NM2A = 
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
* ---------------------------------------------------------------------
*
*                                                                 __
* ..This is the regular non-singlet piece for the diff. F3(nu)-F3(nu), 
*    corresponding to C3NSP+C3NSN in WvN's program. For the NF^0 piece
*    7 numerical coefficients are fitted to his results, the ones of
*    ln^3(1-y) and ln^2(1-y) are taken over from C3NM2A. The NF piece
*    is also the same as in C3NM2A.
*
       FUNCTION C3NP2A (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       C3NP2A = 
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
* ---------------------------------------------------------------------
*
*
* ..This is the singular NS piece, denoted by SOFT2 in WvN's program. 
*    It is the same for all F2 and F3 cases. The numerical coefficients 
*    are exact, but truncated.
*
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
* ---------------------------------------------------------------------
*
*                                        __
* ..This is the 'local' NS piece for the nu+nu F3, denoted by COR2 in 
*    WvN's program. The numerical coefficients of the logs are exact,
*    but truncated, the constant one (from the delta-function) is 
*    slightly adjusted (- 0.104 + 0.013 NF) using the lowest moments.
*
       FUNCTION C3NM2C (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL1 = LOG (1.-Y)
*
       C3NM2C = 
     1          + 3.55555 * DL1**4 - 20.4444 * DL1**3 - 15.5525 * DL1**2
     2          + 188.64 * DL1 - 338.531 - 0.104 
     3        + NF * (0.592593 * DL1**3 - 4.2963 * DL1**2 
     4          + 6.3489 * DL1 + 46.844 + 0.013)
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*                                        __
* ..This is the 'local' NS piece for the nu-nu F3, also given by COR2 in
*    WvN's program. The numerical coefficients of the logs are exact,
*    but truncated, the constant one (from the delta-function) is 
*    slightly adjusted (- 0.152 + 0.013 NF) using the lowest moments.
*
       FUNCTION C3NP2C (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL1 = LOG (1.-Y)
*
       C3NP2C = 
     1          + 3.55555 * DL1**4 - 20.4444 * DL1**3 - 15.5525 * DL1**2
     2          + 188.64 * DL1 - 338.531  - 0.152 
     3        + NF * (0.592593 * DL1**3 - 4.2963 * DL1**2 
     4          + 6.3489 * DL1 + 46.844 + 0.013)
*
       RETURN
       END
*
* =================================================================av==
*
* ..File: xclsg2p.f    FL_S
*
*
* ..Calculation of the 2-loop x-space MS(bar) coefficient functions 
*    for F2 via compact parametrizations involving only logarithms.
*    Singlet, mu_r = mu_f = Q. Expansion parameter: alpha_s/(4 pi).
*
*  ..The relative accuracy of the coefficient functions, as well as of 
*    the convolution results, amounts to a few thousandth.
*
*  ..Reference: W.L. van Neerven and A. Vogt, 
*               hep-ph/0006154 = Nucl. Phys. B588 (2000) 345
*  ..The user should also cite the original calculations,
*      J. Sanchez Guillen et al, Nucl. Phys. B353 (1991) 337  and 
*      E.B. Zijlstra and W.L. van Neerven, Phys. Lett. B273 (1991) 476 
*     
* 
* =====================================================================
*
*
* ..This is the pure singlet piece, denoted by C2S in WvN's program. 
*    Six numerical coefficients (all but the one of 1/y, which is 
*    exact up to truncation) are fitted to his results, using x values 
*    between 10^-6 and 1-10^-6.
*
       FUNCTION CLS2A (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       CLS2A = NF * ( (15.94 - 5.212 * Y) * (1.-Y)**2 * DL1
     1         + (0.421 + 1.520 * Y) * DL**2 + 28.09 * (1.-Y) * DL
     2         - (2.370/Y - 19.27) * (1.-Y)**3 )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the gluon contribution, denoted by C2G2 in WvN's program. 
*    Six numerical coefficients are fitted as above, the one of 1/y is
*    exact up to truncation.
*
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
* =================================================================av==
*
* ..File: xclns2p.f    FL_NS
*
*
* ..Calculation of the 2-loop x-space MS(bar) coefficient functions 
*    for FL via compact parametrizations involving only logarithms.
*    Non-singlet, mu_r = mu_f = Q. Expansion parameter: alpha_s/(4 pi).
*
*  ..The distributions (in the mathematical sense) are given as in eq.
*    (B.26) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.
*    The name-endings A, B, and C of the functions below correspond to 
*    the kernel superscripts [2], [3], and [1] in that equation.
*
*  ..The relative accuracy of the coefficient functions, as well as of 
*    the convolution results, amounts to a few thousandth.
*
*  ..Reference: W.L. van Neerven and A. Vogt, 
*               hep-ph/9907472 = Nucl. Phys. B568 (2000) 263
*  ..The user should also cite the original calculation,
*      J. Sanchez Guillen et al, Nucl. Phys. B353 (1991) 337.
*    
* 
* =====================================================================
*
*
* ..This is the regular non-singlet piece for the electromagnetic F2, 
*    corresponding to CLNSP+C2NSM in W. van Neerven's program. The 
*    8 numerical coefficients are fitted to his results, using x values 
*    between 10^-6 and 1-10^-6. 
*
       FUNCTION CLNN2A (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       CLNN2A = 
     1          - 40.41 + 97.48 * Y
     2          + (26.56 * Y - 0.031) * DL**2 - 14.85 * DL 
     3          + 13.62 * DL1**2 - 55.79 * DL1 - 150.5 * DL * DL1 
     4        + NF * 16./27.D0 * ( 6.* Y*DL1 - 12.* Y*DL - 25.* Y + 6.)
*
       RETURN
       END
* =====================================================================
*
* ..Derivartive with respect to NF of CLNN2A
*
       FUNCTION CLNN2A_DNF (Y)
       IMPLICIT REAL*8 (A-Z)
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       CLNN2A_DNF = 
     1        + 16./27.D0 * ( 6.* Y*DL1 - 12.* Y*DL - 25.* Y + 6.)
*
       RETURN
       END
* ---------------------------------------------------------------------
*
*
* ..This is the regular non-singlet piece for the odd-moment (CC) FL, 
*    corresponding to CLNSP-CLNSN in WvN's program. The 8 numerical 
*    coefficients the NF^0 piece are fitted to his results. The NF part
*    is the same as in CLNN2A.
*
       FUNCTION CLNC2A (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       CLNC2A = 
     1          - 52.27 + 100.8 * Y
     2          + (23.29 * Y - 0.043) * DL**2 - 22.21 * DL 
     3          + 13.30 * DL1**2 - 59.12 * DL1 - 141.7 * DL * DL1 
     4        + NF * 16./27.D0 * ( 6.* Y*DL1 - 12.* Y*DL - 25.* Y + 6.)
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' NS piece for the e.m. FL, with no counterpart 
*    in WvN's program, as it does not exist in the exact expressions.
*    The value is fixed from the lowest integer moments.
*
       FUNCTION CLNN2C (Y)
       IMPLICIT REAL*8 (A-Z)
*
       CLNN2C = -0.164
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' NS piece for the CC FL, with no counterpart 
*    in WvN's program, as it does not exist in the exact expressions.  
*    The value is fixed from the lowest integer moments.
*
       FUNCTION CLNC2C (Y)
       IMPLICIT REAL*8 (A-Z)
*
       CLNC2C = -0.150
*
       RETURN
       END
*
* =================================================================av==
