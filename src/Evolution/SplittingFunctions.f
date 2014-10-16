************************************************************************
*
*     SplittingFunctions.f:
*
*     Collections of the QCD splitting fuctions up to NNLO.
*
************************************************************************
*                           __
* ..The parametrized 3-loop MS non-singlet splitting functions P^(2)
*    for the evolution of unpolarized partons densities, mu_r = mu_f.
*    The expansion parameter is alpha_s/(4 pi).
*
* ..The distributions (in the mathematical sense) are given as in eq.
*    (B.26) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.
*    The name-endings A, B, and C of the functions below correspond to
*    the kernel superscripts [2], [3], and [1] in that equation.
*
*  ..The relative accuracy of these parametrizations, as well as of
*    the convolution results, is better than one part in thousand.
*
* ..References: S. Moch, J. Vermaseren and A. Vogt,
*               hep-ph/0209100 = Nucl. Phys. B646 (2002) 181,
*               hep-ph/0403192 = Nucl. Phys. B688 (2004) 101
*
* =====================================================================
*
*
* ..This is the regular piece of P2_NS+.  The rational coefficients are
*    exact, the rest has been fitted for x between 10^-6 and 1 - 10^-6.
*    The N_f^2 part is exact and was first determined in N-space by 
*    J.A. Gracey in Phys. Lett. B322 (1994) 141.
*
       FUNCTION P2NSPA (Y, NF)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
       D81 = 1./81.D0
*
       P2NSPA =   1641.1 - 3135.* Y + 243.6 * Y**2 - 522.1 * Y**3
     ,            + 128.*D81 * DL**4 + 2400.*D81 * DL**3
     ,            + 294.9 * DL**2 + 1258.* DL
     ,            + 714.1 * DL1 + DL*DL1 * (563.9 + 256.8 * DL)
     ,        + NF * ( -197.0 + 381.1 * Y + 72.94 * Y**2 + 44.79 * Y**3
     ,            - 192.*D81 * DL**3  - 2608.*D81 * DL**2 - 152.6 * DL
     ,            - 5120.*D81 * DL1 - 56.66 * DL*DL1 - 1.497 * Y*DL**3 )
     ,        + NF**2 * ( 32.* Y*DL/(1.-Y) * (3.* DL + 10.) + 64.
     ,            + (48.* DL**2 + 352.* DL + 384.) * (1.-Y) ) * D81
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the regular piece of P2_NS-.  The rational coefficients are 
*    exact, the rest has been fitted for x between 10^-6 and 1 - 10^-6.
*    The N_f^2 part is exact (and identical to that of P2_NS+). 
*
       FUNCTION P2NSMA (Y, NF)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
       D81 = 1./81.D0
*
       P2NSMA =   1860.2 - 3505.* Y + 297.0 * Y**2 - 433.2 * Y**3
     ,            + 116.*D81 * DL**4 + 2880.*D81 * DL**3 
     ,            + 399.2 * DL**2 + 1465.2 * DL
     ,            + 714.1 * DL1 + DL*DL1 * (684.0 + 251.2 * DL)
     ,        + NF * ( -216.62 + 406.5 * Y + 77.89 * Y**2 + 34.76 * Y**3
     ,            - 256.*D81 * DL**3  - 3216.*D81 * DL**2 - 172.69 * DL 
     ,            - 5120.*D81 * DL1 - 65.43 * DL*DL1 - 1.136 * Y*DL**3 )
     ,        + NF**2 * ( 32.* Y*DL/(1.-Y) * (3.* DL + 10.) + 64.
     ,            + (48.* DL**2 + 352.* DL + 384.) * (1.-Y) ) * D81
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the singular piece of both P2_NS+ and P2_NS-. It is exact 
*    up to the truncation of the irrational coefficients.
*
       FUNCTION P2NSB (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       P2NSB = ( 1174.898 - NF * 183.187 - NF**2 * 64./81.D0 ) / (1.-Y)
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' piece of P2_NS+. The coefficients of delta(1-x)
*    have been partly shifted relative to the exact (truncated) values.
*
       FUNCTION P2NSPC (Y, NF)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
*
       DL1 = LOG (1.-Y)
*
       P2NSPC =       1174.898 * DL1 + 1295.624 - 0.24
     ,        - NF * ( 183.187 * DL1 + 173.938 - 0.011 )
     ,        + NF**2 * ( - 64./81.D0 * DL1 + 1.13067 )
*
       RETURN
       END
*
*
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' piece of P2_NS-. The coefficients of delta(1-x) 
*    have been partly shifted relative to the exact (truncated) values.
*
       FUNCTION P2NSMC (Y, NF)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
*
       DL1 = LOG (1.-Y)
*
       P2NSMC =       1174.898 * DL1 + 1295.624 - 0.154
     ,        - NF * ( 183.187 * DL1 + 173.938  - 0.005 )
     ,        + NF**2 * ( - 64./81.D0 * DL1 + 1.13067 )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is P2_NSS, the difference of P2_NSV and P2_NS-.
*
       FUNCTION P2NSSA (Y, NF)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       D27 = 1./27.D0
       DL  = LOG (Y)
       Y1  = 1.- Y
       DL1 = LOG (Y1)
*
       P2NSSA = Y1* ( 151.49 + 44.51 * Y - 43.12 * Y**2 + 4.820 * Y**3 )
     1          + 40.*D27 * DL**4 - 80.*D27 * DL**3 + 6.892 * DL**2 
     2          + 178.04 * DL + DL*DL1 * ( - 173.1 + 46.18 * DL )
     4          + Y1*DL1 * ( - 163.9 / Y - 7.208 * Y ) 
*
       P2NSSA  = NF * P2NSSA
*
       RETURN
       END
*
* =================================================================av==
*
* ..The 1- and 2-loop MS(bar) non-singlet splitting functions P_NS^(2) 
*    for the evolution of unpolarized partons densities, mu_r = mu_f.
*    The expansion parameter is alpha_s/(4 pi).
* 
* ..The distributions (in the mathematical sense) are given as in eq.
*    (B.26) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.
*    The name-endings A, B, and C of the functions below correspond to 
*    the kernel superscripts [2], [3], and [1] in that equation.
*
* ..The code uses the package of Gehrmann and Remiddi for the harmonic
*    polylogarithms published in hep-ph/0107173 = CPC 141 (2001) 296.
*
* =====================================================================
*
*
* ..This is the regular 2-loop piece for P_NS^+. 
*
      FUNCTION X1NSPA (X, NF)
*
      IMPLICIT REAL*8 (A - Z)
*
      include "../commons/consts.h"
*
      INTEGER NF
*
      PARAMETER ( Z2 = 1.6449 34066 84822 64365 D0,
     1            Z3 = 1.2020 56903 15959 42854 D0 )
*
* ..The soft coefficient for use in X2NSB and X2NSC
*
      COMMON / P1SOFT / A2
*
* ...Colour factors and some abbreviations
*
      CF  = 4D0/3D0
      CA  = 3D0
      lnx = dlog(x)
      ln1mx = dlog(1d0 - x)
      pqq   = 2d0 / ( 1d0 - x ) - 1d0 - x
      pqqmx = 2d0 / ( 1d0 + x ) - 1d0 + x
      S2x   = S2(x)
      DM = 1D0/(1D0-X)
*
      gqq1 = 2d0 * CF * NF * ( ( - 1.1111111111111112d0 - ( 2d0 * lnx ) 
     1     / 3d0 ) * pqq - ( 4d0 * ( 1d0 - x ) ) / 3d0 ) 
     2     + 4d0 * CA * CF * ( ( 3.7222222222222223d0 + ( 11d0 * lnx ) 
     3     / 6d0 + lnx**2d0 / 2d0 - pi**2d0 / 6d0 ) * pqq
     4     + ( 20d0 * ( 1d0 - x ) ) / 3d0 + lnx * ( 1d0 + x) )
     5     + 4d0 * CF**2d0 * ( ( ( - 3d0 * lnx ) / 2d0 - 2d0 * ln1mx 
     6     * lnx ) * pqq - 5d0 * ( 1d0 - x ) - ( lnx**2d0 * ( 1d0 
     7     + x ) ) / 2d0 - lnx * ( 1.5d0 + ( 7d0 * x ) / 2d0 ) )
     8     + 4d0 * CF * ( CF - CA / 2d0 ) * ( 2d0 * pqqmx * S2x 
     9     + 4d0 * ( 1d0 - x ) + 2d0 * lnx * ( 1d0 + x ) )
*
* ...The soft (`+'-distribution) part of the splitting function
*
       A2 = - 40D0/9D0*cf*nf + 268D0/9D0*ca*cf - 8D0*z2*ca*cf
*
       GQQ1L = DM * A2
*
* ...The regular piece of the coefficient function
*
       X1NSPA = GQQ1 - GQQ1L
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the regular 2-loop piece for P_NS^-. 
*
      FUNCTION X1NSMA (X, NF)
*
      IMPLICIT REAL*8 (A - Z)
*
      include "../commons/consts.h"
*
      INTEGER NF
*
      PARAMETER ( Z2 = 1.6449 34066 84822 64365 D0,
     1            Z3 = 1.2020 56903 15959 42854 D0 )
*
* ..The soft coefficient for use in X2NSB and X2NSC
*
      COMMON / P1SOFT / A2
*
* ...Colour factors and some abbreviations
*
      CF  = 4D0/3D0
      CA  = 3D0
      lnx = dlog(x)
      ln1mx = dlog(1d0 - x)
      pqq   = 2d0 / ( 1d0 - x ) - 1d0 - x
      pqqmx = 2d0 / ( 1d0 + x ) - 1d0 + x
      S2x   = S2(x)
      DM = 1D0/(1D0-X)
*
      gqq1 = 2d0 * CF * NF * ( ( - 1.1111111111111112d0 - ( 2d0 * lnx ) 
     1     / 3d0 ) * pqq - ( 4d0 * ( 1d0 - x ) ) / 3d0 ) 
     2     + 4d0 * CA * CF * ( ( 3.7222222222222223d0 + ( 11d0 * lnx ) 
     3     / 6d0 + lnx**2d0 / 2d0 - pi**2d0 / 6d0 ) * pqq
     4     + ( 20d0 * ( 1d0 - x ) ) / 3d0 + lnx * ( 1d0 + x) )
     5     + 4d0 * CF**2d0 * ( ( ( - 3d0 * lnx ) / 2d0 - 2d0 * ln1mx 
     6     * lnx ) * pqq - 5d0 * ( 1d0 - x ) - ( lnx**2d0 * ( 1d0 
     7     + x ) ) / 2d0 - lnx * ( 1.5d0 + ( 7d0 * x ) / 2d0 ) )
     8     - 4d0 * CF * ( CF - CA / 2d0 ) * ( 2d0 * pqqmx * S2x 
     9     + 4d0 * ( 1d0 - x ) + 2d0 * lnx * ( 1d0 + x ) )
*
* ...The soft (`+'-distribution) part of the splitting function
*
       A2 = - 40D0/9D0*cf*nf + 268D0/9D0*ca*cf - 8D0*z2*ca*cf
*
       GQQ1L = DM * A2
*
* ...The regular piece of the coefficient function
*
       X1NSMA = GQQ1 - GQQ1L
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the singular (soft) piece.
*
       FUNCTION X1NSB (Y)
       IMPLICIT REAL*8 (A - Z)
*
       COMMON / P1SOFT / A2
*
       X1NSB  = A2/(1D0-Y)
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' piece.
*
       FUNCTION X1NSC (Y, NF)
*
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
       PARAMETER ( Z2 = 1.6449 34066 84822 64365 D0,
     ,             Z3 = 1.2020 56903 15959 42854 D0 )
c*
c       COMMON / P1SOFT / A2
*
* ...Colour factors
*
       CF  = 4D0/3D0
       CA  = 3D0
*
* ...The coefficient of delta(1-x)
*
       P1DELT = 
     &     - 1D0/3D0*cf*nf
     &     + 3D0/2D0*cf**2
     &     + 17D0/6D0*ca*cf
     &     + 24D0*z3*cf**2
     &     - 12D0*z3*ca*cf
     &     - 8D0/3D0*z2*cf*nf
     &     - 12D0*z2*cf**2
     &     + 44D0/3D0*z2*ca*cf
*
* ...The soft (`+'-distribution) part of the splitting function
*
       A2 = - 40D0/9D0*cf*nf + 268D0/9D0*ca*cf - 8D0*z2*ca*cf
*
       X1NSC = LOG (1D0-Y) * A2 + P1DELT
*
       RETURN
       END
*
* =====================================================================
*
*
* ..This is the regular 1-loop piece. 
*
       FUNCTION X0NSA (X)
       IMPLICIT REAL*8 (A - Z)
*
       CF = 4D0/3D0
       X0NSA = - 2D0*CF * (1.+ X)
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the singular (soft) piece.
*
       FUNCTION X0NSB (Y)
       IMPLICIT REAL*8 (A - Z)
*
       CF = 4D0/3D0
       X0NSB = 4D0*CF/(1D0-Y)
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' piece.
*
       FUNCTION X0NSC (Y)
       IMPLICIT REAL*8 (A - Z)
*
       CF = 4D0/3D0
       X0NSC = 4D0*CF * LOG (1D0-Y) + 3D0*CF
*
       RETURN
       END
*
* =================================================================av==
*
* ..The parametrized 3-loop MS singlet splitting functions P^(2) for 
*    the evolution of unpol. singlet parton densities at mu_r = mu_f.
*    The expansion parameter is alpha_s/(4 pi).
*
* ..The distributions (in the mathematical sense) are given as in eq.
*   (B.27) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.
*   The name-endings A, B, and C of the functions below correspond to 
*   the kernel superscripts [2], [3], and [1] in that equation.
*
* ..The relative accuracy of these parametrisations, as well as of
*    the convolution results, is better than one part in thousand.

* ..The coefficients of 1/(1-x)_+, (ln x)/x and 1/x are exact (up
*    to a truncation of irrational coefficients).  Furthermore all
*    coefficients written as fractions (e.g., 160./27.D0) are exact.
*    The other terms at x < 1 have fitted to the exact results for x 
*    between 10^-6 and 1 - 10^-6.  The coefficient of delta(1-x) of
*    P_gg^(2) have been slightly adjusted using the second moments.
*
* ..References: S. Moch, J. Vermaseren and A. Vogt,
*               hep-ph/0403192 = Nucl. Phys. B688 (2004) 101
*               A. Vogt, S. Moch and J. Vermaseren,
*               hep-ph/0404111 = Nucl. Phys. B691 (2004) 129
* 
* =====================================================================
*
*
* ..The (regular) pure-singlet splitting functions P_ps^(2).
*    P_qq^(2) is obtained by adding the non-singlet quantity P_NS^(2)+.
*    A parametrization of the latter is provided in the file  xpns2p.f.

       FUNCTION P2PSA (Y, NF)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       P2PS1 = - 3584./(27.D0*Y) * DL - 506.0/ Y + 160./27.D0 * DL**4
     ,         - 400./9.D0 * DL**3 + 131.4 * DL**2 - 661.6 * DL
     ,         - 5.926  * DL1**3 - 9.751 * DL1**2 - 72.11 * DL1
     ,         + 177.4 + 392.9 * Y - 101.4 * Y**2 - 57.04 * DL*DL1
       P2PS2 =   256./(81.*Y) + 32./27.D0 * DL**3 + 17.89 * DL**2
     ,         + 61.75 * DL + 1.778 * DL1**2 + 5.944 * DL1 + 100.1
     ,         - 125.2 * Y + 49.26 * Y**2 - 12.59 * Y**3 
     ,         - 1.889 * DL*DL1 
*
       P2PSA = (1.-Y) * NF * ( P2PS1 + NF * P2PS2 )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The gluon->quark splitting functions P_qg^(2).
*
       FUNCTION P2QGA (Y, NF)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER  NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       P2QG1 = - 896./(3.D0*Y) * DL - 1268.3 / Y + 536./27.D0 * DL**4 
     ,         - 44./3.D0 * DL**3 + 881.5 * DL**2 + 424.9 * DL 
     ,         + 100./27.D0 * DL1**4 - 70./9.D0 * DL1**3 
     ,         - 120.5 * DL1**2 + 104.42 * DL1
     ,         + 2522. - 3316.* Y + 2126.* Y**2
     ,         + DL*DL1 * (1823. - 25.22 * DL) - 252.5 * Y*DL**3  
       P2QG2 =   1112./(243.D0*Y) - 16./9.D0 * DL**4 
     ,         - 376./27.D0 * DL**3 - 90.8 * DL**2 - 254.0 * DL   
     ,         + 20./27.D0 * DL1**3 + 200./27.D0 * DL1**2 - 5.496 * DL1
     ,         - 252.0  + 158.0 * Y + 145.4 * Y**2 - 139.28 * Y**3
     ,         - DL*DL1 * ( 53.09  + 80.616 * DL) - 98.07 * Y*DL**2
     ,         + 11.70 * Y*DL**3
* 
       P2QGA = NF * ( P2QG1 + NF * P2QG2 )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The quark->gluon splitting functions P_gq^(2).  P2GQ2 is exact.
*
       FUNCTION P2GQA (Y, NF)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       P2GQ0 =   1189.3 * DL/Y + 6163.1 / Y - 4288./81.D0 * DL**4
     ,         + 1568./9.D0 * DL**3 - 1794. * DL**2 + 4033. * DL
     ,         + 400./81.D0 * DL1**4 + 2200./27.D0 * DL1**3
     ,         + 606.3 * DL1**2 + 2193.* DL1 
     ,         - 4307. + 489.3 * Y + 1452.* Y**2 + 146.0 * Y**3
     ,         - 447.3 * DL**2*DL1 - 972.9 * Y*DL**2
       P2GQ1 =   71.082 * DL/Y  - 46.41 / Y + 128./27.D0 * DL**4
     ,         + 704/81.D0 * DL**3 + 20.39  * DL**2 + 174.8 * DL
     ,         - 400./81.D0 * DL1**3 - 68.069 * DL1**2 - 296.7 * DL1
     ,         - 183.8 + 33.35 * Y - 277.9 * Y**2 + 108.6 * Y*DL**2
     ,         - 49.68 * DL*DL1
       P2GQ2 = ( 64. * ( - 1./Y + 1. + 2.* Y)
     ,         + 320.* DL1 * ( 1./Y - 1. + 0.8 * Y)
     ,         + 96.* DL1**2 * ( 1./Y - 1. + 0.5 * Y) ) / 27.D0
*
       P2GQA = ( P2GQ0 + NF * (P2GQ1 + NF * P2GQ2) )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The regular piece of the gluon-gluon splitting function P_gg^(2).
*
       FUNCTION P2GGA (Y, NF)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       P2GGA0 = 2675.8 * DL/Y + 14214./ Y - 144. * DL**4 + 72. * DL**3
     1          - 7471. * DL**2 + 274.4 * DL + 3589. * DL1 - 20852. 
     2          + 3968.* Y - 3363. * Y**2 + 4848. * Y**3 
     3          + DL*DL1 * ( 7305. + 8757. * DL )
       P2GGA1 = 157.27 * DL/Y + 182.96 / Y + 512./27.D0 * DL**4
     1          + 832./9.D0 * DL**3 + 491.3 * DL**2 + 1541. * DL
     2          - 320.0 * DL1 - 350.2 + 755.7 * Y - 713.8 * Y**2 
     3          + 559.3 * Y**3 + DL*DL1 * ( 26.15 - 808.7 * DL )
       P2GGA2 = - 680./(243.D0 * Y) - 32./27.D0 * DL**3 + 9.680 * DL**2
     1          - 3.422 * DL - 13.878 + 153.4 * Y - 187.7 * Y**2 
     2          + 52.75 * Y**3 - DL*DL1 * (115.6 - 85.25* Y + 63.23* DL)
*
       P2GGA = P2GGA0 + NF * ( P2GGA1 + NF * P2GGA2 )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The singular piece of the gluon-gluon splitting function P_gg^(2).
*
       FUNCTION P2GGB (Y, NF)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       P2GGB = ( 2643.521 - NF * 412.172 - NF**2 * 16./9.D0 ) / ( 1.-Y) 
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The 'local' piece of the gluon-gluon splitting function P_gg^(2).  
*
       FUNCTION P2GGC (Y, NF)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL1 = LOG (1.-Y)
*
       P2GGC =       2643.521 * DL1 + 4425.448 + 0.446
     ,       - NF * ( 412.172 * DL1 +  528.720 + 0.003 )
     ,       + NF**2 * ( - 16./9.D0 * DL1 + 6.4630)
*
       RETURN
       END
*
* =================================================================av==
*
* ..The 1- and 2-loop MS(bar) singlet splitting functions  P_ij^(2) 
*    for the evolution of unpolarized partons densities, mu_r = mu_f.
*    The expansion parameter is alpha_s/(4 pi).
* 
* ..The distributions (in the mathematical sense) are given as in eq.
*    (B.26) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.
*    The name-endings A, B, and C of the functions below correspond to 
*    the kernel superscripts [2], [3], and [1] in that equation.
*
* ..The code uses the package of Gehrmann and Remiddi for the harmonic
*    polylogarithms published in hep-ph/0107173 = CPC 141 (2001) 296.
*
* =====================================================================
*
*
* ..The 2-loop pure-singlet splitting functions P_ps^(1)
*
      FUNCTION X1PSA (X, NF)
*
      IMPLICIT REAL*8 (A - Z)
      INTEGER NF
*
* ...Colour factors and abbreviation
*
      CF = 4D0/3D0
      DX = 1D0/X
      LNX = DLOG(X)
      HR200 = LNX * LNX / 2D0
*
* ...The splitting function in terms of the harmonic polylogs
*
      X1PSA =
     &  + nf*cf * (  - 8D0 + 24D0*x - 224D0/9D0*x**2 + 80D0/9D0*
     &    dx + 4D0*LNX + 20D0*LNX*x + 32D0/3D0*LNX*x**2 -
     &    8D0*HR200 - 8D0*HR200*x )
*
      RETURN
      END
*
* ---------------------------------------------------------------------
*
*
* ..The 2-loop gluon->quark splitting functions P_qg^(1)
*
      FUNCTION X1QGA (X, NF)
*
      IMPLICIT REAL*8 (A - Z)
*
      include "../commons/consts.h"
*
      INTEGER NF
*
      CF = 4D0/3D0
      CA = 3D0
      lnx = dlog(x)
      ln1mx = dlog(1d0 - x)
      pqg   = x**2d0 + ( 1d0 - x )**2d0
      pqgmx = x**2d0 + ( 1d0 + x )**2d0
      S2x   = S2(x)
*
      X1QGA = 2d0 * CF * NF * ( 4d0  + 4d0 * ln1mx + ( 10d0 - 4d0 
     1      * ( ln1mx - lnx ) + 2d0 * ( - ln1mx + lnx )**2d0 
     2      - 2d0 * pi**2d0 / 3d0 ) * pqg - lnx * ( 1d0 - 4d0 * x )
     3      - lnx**2d0 * ( 1d0  - 2d0 * x ) - 9d0 * x )
     4      + 2d0 * CA * NF * ( 20.22222222222222d0 - 4d0 * ln1mx
     5      + ( - 24.22222222222222d0 + 4d0 * ln1mx - 2d0 * ln1mx**2d0
     6      + ( 44d0 * lnx ) /3d0 - lnx**2d0 + pi**2d0 / 3d0 ) * pqg 
     7      + 2d0 * pqgmx * S2x + 40d0 / ( 9d0 * x ) 
     8      + ( 14d0 * x ) / 9d0 - lnx**2d0 * ( 2d0 + 8d0 * x ) 
     9      + lnx * ( - 12.666666666666666d0 + ( 136d0 * x ) / 3d0 ) )
*
      RETURN
      END
*
* ---------------------------------------------------------------------
*
*
* ..The 2-loop quark->gluon splitting functions P_gq^(1)
*
      FUNCTION X1GQA (X, NF)
*
      IMPLICIT REAL*8 (A - Z)
*
      include "../commons/consts.h"
*
      INTEGER NF
*
* ...Colour factors and abbreviation
*
      CF = 4D0/3D0
      CA = 3D0
      lnx = dlog(x)
      ln1mx = dlog(1d0 - x)
      pgq   = ( 1d0 + ( 1d0 - x )**2d0 ) / x
      pgqmx = - ( 1d0 + ( 1d0 + x )**2d0 ) / x
      S2x   = S2(x)
*
      X1GQA = 2d0 * CF * NF * ( - ( ( 2.2222222222222223d0 
     1      + ( 4d0 * ln1mx ) / 3d0 ) * pgq ) - ( 4d0 * x ) / 3d0 )
     2      + 4d0 * CF**2d0 * ( - 2.5d0 - ( 3d0 * ln1mx + ln1mx**2d0 )
     3      * pgq - lnx**2d0 * ( 1d0 - x / 2d0 ) - ( 7d0 * x ) / 2d0 
     4      - 2d0 * ln1mx * x + lnx * ( 2d0 + ( 7d0 * x ) / 2d0 ) ) 
     5      + 4d0 * CA * CF * ( 3.111111111111111d0 + pgq * ( 0.5 d0 
     6      + ( 11d0 * ln1mx ) / 3d0 + ln1mx**2d0 - 2d0 * ln1mx * lnx
     7      + lnx**2d0 / 2d0 - pi**2d0 / 6d0 ) + pgqmx * S2x 
     8      + ( 65d0 * x ) / 18d0 + 2d0 * ln1mx * x + ( 44d0 * x**2d0 ) 
     9      / 9d0 + lnx**2d0 * ( 4d0 + x ) - lnx * ( 12d0 + 5d0 * x 
     1      + ( 8d0 * x**2d0 ) / 3d0 ) )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The regular piece of the 2-loop gg splitting function P_gg^(1) 
*
      FUNCTION X1GGA (X, NF)
*
      IMPLICIT REAL*8 (A - Z)
*
      include "../commons/consts.h"
*
      INTEGER NF
      PARAMETER ( Z2 = 1.6449 34066 84822 64365 D0,
     1            Z3 = 1.2020 56903 15959 42854 D0 )
*
* ..The soft coefficient for use in X1GGB and X1GGC
*
       COMMON / P1GSOFT / A2G
*
* ...Colour factors and abbreviation
*
      CF = 4D0/3D0
      CA = 3D0
      lnx = dlog(x)
      ln1mx = dlog(1d0 - x)
      pgg   = ( 1d0 / ( 1d0 - x ) +  1d0 / x - 2d0 + x * ( 1d0 - x ) )
      pggmx = ( 1d0 / ( 1d0 + x ) -  1d0 / x - 2d0 - x * ( 1d0 + x ) )
      S2x   = S2(x)
      DM    = 1D0/(1D0-X)
*
      ggg1  = 2d0 * CF * NF * ( - 16d0 + 4d0 / ( 3d0 * x ) + 8d0 * x 
     1      + ( 20d0 * x**2d0 ) / 3d0 - lnx**2d0 * ( 2d0 + 2d0 * x ) 
     2      - lnx * ( 6d0 + 10d0 * x ) )
     3      + 2d0 * CA * NF * ( 2d0 - ( 20d0 * pgg ) / 9d0 - 2d0 * x 
     4      - ( 4d0 * lnx * ( 1d0 + x ) ) / 3d0 + ( 26d0 * ( 
     5      - ( 1d0 / x ) + x**2d0 ) ) / 9d0 ) 
     6      + 4d0 * CA**2d0 * ( pgg * ( 7.444444444444445d0 
     7      - 4d0 * ln1mx * lnx + lnx**2d0 - pi**2d0 / 3d0 ) 
     8      + 2d0 * pggmx * S2x + ( 27d0 * ( 1d0 - x ) ) / 2d0 
     9      + 4d0 * lnx**2d0 * ( 1d0 + x ) + ( 67d0 * ( - ( 1d0 / x ) 
     1      + x**2d0 ) ) / 9d0 - lnx * ( 8.333333333333334d0 
     2      - ( 11d0 * x ) / 3d0 + ( 44d0 * x**2d0 ) / 3d0 ) )
*
* ...The soft (`+'-distribution) part of the splitting function
*
       A2G = - 40D0/9D0*ca*nf + 268D0/9D0*ca**2 - 8D0*z2*ca**2
*
       GGG1L = DM * A2G
*
* ...The regular piece of the coefficient function
*
       X1GGA = GGG1 - GGG1L
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the singular (soft) piece.
*
       FUNCTION X1GGB (Y)
       IMPLICIT REAL*8 (A - Z)
*
       COMMON / P1GSOFT / A2G
*
       X1GGB  = A2G/(1D0-Y)
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' piece.
*
       FUNCTION X1GGC (Y, NF)
*
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
       PARAMETER ( Z2 = 1.6449 34066 84822 64365 D0,
     ,             Z3 = 1.2020 56903 15959 42854 D0 )
c*
c       COMMON / P1GSOFT / A2G
*
* ...Colour factors
*
       CF  = 4D0/3D0
       CA  = 3D0
*
* ...The coefficient of delta(1-x)
*
       P1DELT = 
     ,    - 2D0*cf*nf
     ,    - 8D0/3D0*ca*nf
     ,    + 32D0/3D0*ca**2
     ,    + 12D0*z3*ca**2
*
* ...The soft (`+'-distribution) part of the splitting function                                                                                                                   *
       A2G = - 40D0/9D0*ca*nf + 268D0/9D0*ca**2 - 8D0*z2*ca**2
*
       X1GGC = DLOG (1D0-Y) * A2G + P1DELT
*
       RETURN
       END
*
* =====================================================================
*
*
* ..The 1-loop gluon->quark splitting functions P_qg^(0)
*
       FUNCTION X0QGA (X, NF)
*
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
*
       X0QGA = 2D0* NF * ( 1. - 2. * X + 2. * X**2 )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The 1-loop quark->gluon splitting functions P_gq^(0)
*
       FUNCTION X0GQA (X)
*
       IMPLICIT REAL*8 (A - Z)
*
       CF = 4D0/3D0
       X0GQA = 4D0*CF * ( - 1. + 0.5 * X + 1D0/X )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The regular piece of the 1-loop gg splitting function P_gg^(0)
*
       FUNCTION X0GGA (X)
       IMPLICIT REAL*8 (A - Z)
*
       CA = 3D0
       X0GGA = 4D0*CA * ( - 2. + X - X**2 + 1D0/X )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the singular (soft) piece.
*
       FUNCTION X0GGB (X)
       IMPLICIT REAL*8 (A - Z)
*
       CA = 3D0
       X0GGB = 4D0*CA / (1D0-X)
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' piece.
*
       FUNCTION X0GGC (X,NF)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
*
       CA = 3D0
       X0GGC = 4D0*CA * LOG (1D0-X) - 2D0/3D0 * NF + 11D0/3D0 * CA
*
       RETURN
       END
*
* =================================================================av==
      function S2(x)
*
      implicit none
*
      include "../commons/consts.h"
**
*     Input Variables
*
      double precision x
**
*     Internal Variables
*
      double precision lnx
      double precision ddilog!,wgplg
**
*     Output Variables
*
      double precision S2
*
      lnx = dlog(x)
*
      S2 = - 2d0 * ddilog(-x) + lnx**2d0 / 2d0 
     1     - 2d0 * lnx * log(1d0+x) - pi**2 / 6d0
c      S2 = - 2d0 * wgplg(1,1,-x) + lnx**2d0 / 2d0 
c     1     - 2d0 * lnx * log(1d0+x) - pi**2 / 6d0
*
      return
      end
*
************************************************************************
*
*     Time-like splitting functions for the fragmentation function 
*     evolution.
*     Notice that the one-loop time-like splitting functions coincide
*     with the space-like ones.
*
***********************************************************************
*
* ..This is the regular 2-loop piece for P_NS^+. 
*
      FUNCTION X1NSPTA (X, NF)
*
      IMPLICIT REAL*8 (A - Z)
*
      include "../commons/consts.h"
*
      INTEGER NF
*
      PARAMETER ( Z2 = 1.6449 34066 84822 64365 D0,
     1            Z3 = 1.2020 56903 15959 42854 D0 )
*
* ..The soft coefficient for use in X2NSBT and X2NSCT
*
      COMMON / P1SOFTT / A2T
*
* ...Colour factors and some abbreviations
*
      CF  = 4D0/3D0
      CA  = 3D0
      lnx = dlog(x)
      ln1mx = dlog(1d0 - x)
      pqq   = 2d0 / ( 1d0 - x ) - 1d0 - x
      pqqmx = 2d0 / ( 1d0 + x ) - 1d0 + x
      S2x   = S2(x)
      DM = 1D0/(1D0-X)
*
      gqq1 = 2d0 * CF * NF * ( ( - 1.1111111111111112d0 - ( 2d0 * lnx ) 
     1     / 3d0 ) * pqq - ( 4d0 * ( 1d0 - x ) ) / 3d0 ) 
     2     + 4d0 * CA * CF * ( ( 3.7222222222222223d0 + ( 11d0 * lnx ) 
     3     / 6d0 + lnx**2d0 / 2d0 - pi**2d0 / 6d0 ) * pqq
     4     + ( 20d0 * ( 1d0 - x ) ) / 3d0 + lnx * ( 1d0 + x ) )
     5     + 4d0 * CF**2d0 * ( ( ( 3d0 * lnx ) / 2d0 + 2d0 * ln1mx 
     6     * lnx - 2d0 * lnx**2d0 ) * pqq - 5d0 * ( 1d0 - x ) 
     7     + ( lnx**2d0 * ( 1d0 + x ) ) / 2d0 
     8     - lnx * ( 3.5d0 + ( 3d0 * x ) / 2d0 ) )
     9     + 4d0 * CF * ( CF - CA / 2d0 ) * ( 2d0 * pqqmx * S2x 
     1     + 4d0 * ( 1d0 - x ) + 2d0 * lnx * ( 1d0 + x ) )
*
* ...The soft (`+'-distribution) part of the splitting function
*
       A2T = - 40D0/9D0*cf*nf + 268D0/9D0*ca*cf - 8D0*z2*ca*cf
*
       GQQ1L = DM * A2T
*
* ...The regular piece of the coefficient function
*
       X1NSPTA = GQQ1 - GQQ1L
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the regular 2-loop piece for P_NS^-. 
*
      FUNCTION X1NSMTA (X, NF)
*
      IMPLICIT REAL*8 (A - Z)
*
      include "../commons/consts.h"
*
      INTEGER NF
*
      PARAMETER ( Z2 = 1.6449 34066 84822 64365 D0,
     1            Z3 = 1.2020 56903 15959 42854 D0 )
*
* ..The soft coefficient for use in X2NSBT and X2NSCT
*
      COMMON / P1SOFTT / A2T
*
* ...Colour factors and some abbreviations
*
      CF  = 4D0/3D0
      CA  = 3D0
      lnx = dlog(x)
      ln1mx = dlog(1d0 - x)
      pqq   = 2d0 / ( 1d0 - x ) - 1d0 - x
      pqqmx = 2d0 / ( 1d0 + x ) - 1d0 + x
      S2x   = S2(x)
      DM = 1D0/(1D0-X)
*
      gqq1 = 2d0 * CF * NF * ( ( - 1.1111111111111112d0 - ( 2d0 * lnx ) 
     1     / 3d0 ) * pqq - ( 4d0 * ( 1d0 - x ) ) / 3d0 ) 
     2     + 4d0 * CA * CF * ( ( 3.7222222222222223d0 + ( 11d0 * lnx ) 
     3     / 6d0 + lnx**2d0 / 2d0 - pi**2d0 / 6d0 ) * pqq
     4     + ( 20d0 * ( 1d0 - x ) ) / 3d0 + lnx * ( 1d0 + x) )
     5     + 4d0 * CF**2d0 * ( ( ( 3d0 * lnx ) / 2d0 + 2d0 * ln1mx 
     6     * lnx - 2d0 * lnx**2d0 ) * pqq - 5d0 * ( 1d0 - x ) 
     7     + ( lnx**2d0 * ( 1d0 + x ) ) / 2d0 
     8     - lnx * ( 3.5d0 + ( 3d0 * x ) / 2d0 ) )
     9     - 4d0 * CF * ( CF - CA / 2d0 ) * ( 2d0 * pqqmx * S2x 
     1     + 4d0 * ( 1d0 - x ) + 2d0 * lnx * ( 1d0 + x ) )
*
* ...The soft (`+'-distribution) part of the splitting function
*
       A2T = - 40D0/9D0*cf*nf + 268D0/9D0*ca*cf - 8D0*z2*ca*cf
*
       GQQ1L = DM * A2T
*
* ...The regular piece of the coefficient function
*
       X1NSMTA = GQQ1 - GQQ1L
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the singular (soft) piece.
*
       FUNCTION X1NSTB (Y)
       IMPLICIT REAL*8 (A - Z)
*
       COMMON / P1SOFTT / A2T
*
       X1NSTB  = A2T/(1D0-Y)
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' piece.
*
       FUNCTION X1NSTC (Y, NF)
*
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
       PARAMETER ( Z2 = 1.6449 34066 84822 64365 D0,
     ,             Z3 = 1.2020 56903 15959 42854 D0 )
c*
c       COMMON / P1SOFTT / A2T
*
* ...Colour factors
*
       CF  = 4D0/3D0
       CA  = 3D0
*
* ...The coefficient of delta(1-x)
*
       P1DELT = 
     &     - 1D0/3D0*cf*nf
     &     + 3D0/2D0*cf**2
     &     + 17D0/6D0*ca*cf
     &     + 24D0*z3*cf**2
     &     - 12D0*z3*ca*cf
     &     - 8D0/3D0*z2*cf*nf
     &     - 12D0*z2*cf**2
     &     + 44D0/3D0*z2*ca*cf
*
* ...The soft (`+'-distribution) part of the splitting function
*
       A2T = - 40D0/9D0*cf*nf + 268D0/9D0*ca*cf - 8D0*z2*ca*cf
*
       X1NSTC = LOG (1D0-Y) * A2T + P1DELT
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The 2-loop pure-singlet time-like splitting functions P_ps^(1)
*
      FUNCTION X1PSTA (X, NF)
*
      IMPLICIT REAL*8 (A - Z)
      INTEGER NF
*
* ...Colour factors and abbreviation
*
      CF = 4D0/3D0
      DX = 1D0/X
      LNX = DLOG(X)
      HR200 = LNX * LNX / 2D0
*
* ...The splitting function in terms of the harmonic polylogs
*
      X1PSTA =
     &  + nf*cf * ( - 32D0 + 16D0*x + 224D0/9D0*x**2 - 80D0/9D0*
     &    dx - 20D0*LNX - 36D0*LNX*x - 32D0/3D0*LNX*x**2 +
     &    8D0*HR200 + 8D0*HR200*x )
*
      RETURN
      END
*
* ---------------------------------------------------------------------
*
*
* ..The 2-loop gluon->quark time-like splitting functions P_qg^(1)
*
      FUNCTION X1QGTA (X, NF)
*
      IMPLICIT REAL*8 (A - Z)
*
      include "../commons/consts.h"
*
      INTEGER NF
*
      CF = 4D0/3D0
      CA = 3D0
      lnx = dlog(x)
      ln1mx = dlog(1d0 - x)
      pqg   = x**2d0 + ( 1d0 - x )**2d0
      pqgmx = x**2d0 + ( 1d0 + x )**2d0
      S1x   = - ddilog(1d0-x)
      S2x   = S2(x)
*
      X1QGTA = ( NF**2d0 * ( - 8d0 / 3d0 - ( 16d0 / 9d0 
     1       + 8d0 * lnx / 3d0 + 8d0 * ln1mx / 3d0 ) * pqg )
     2       + 2d0 * CF * NF * ( - 2d0 + 3d0 * x 
     3       + ( - 7d0 + 8d0 * x ) * lnx - 4d0 * ln1mx 
     4       + ( 1d0 - 2d0 * x ) * lnx**2d0 
     5       + ( - 2d0 * ( lnx + ln1mx )**2d0 - 2d0 * ( ln1mx - lnx )
     6       + 16d0 * S1x + 2d0 * pi**2d0 - 10d0 ) * pqg )
     7       + 2d0 * CA * NF * ( - 152d0 / 9d0 + 166d0 * x / 9d0 
     8       - 40d0 / 9d0 / x + ( - 4d0 / 3d0 - 76d0 * x / 3d0 ) * lnx 
     9       + 4d0 * ln1mx + ( 2d0 + 8d0 * x ) * lnx**2d0 
     1       + ( 8d0 * lnx * ln1mx - lnx**2d0 - 4d0 * lnx / 3d0 
     2       + 10d0 * ln1mx / 3d0 + 2d0 * ln1mx**2d0 - 16d0 * S1x 
     3       - 7d0 * pi**2d0 / 3d0 + 178d0 / 9d0 ) * pqg 
     4       + 2d0 * pqgmx * S2x ) ) / 2d0 / NF
*
      RETURN
      END
*
* ---------------------------------------------------------------------
*
*
* ..The 2-loop quark->gluon time-like splitting functions P_gq^(1)
*
      FUNCTION X1GQTA (X, NF)
*
      IMPLICIT REAL*8 (A - Z)
*
      include "../commons/consts.h"
*
      INTEGER NF
*
* ...Colour factors and abbreviation
*
      CF = 4D0/3D0
      CA = 3D0
      lnx = dlog(x)
      ln1mx = dlog(1d0 - x)
      pgq   = ( 1d0 + ( 1d0 - x )**2d0 ) / x
      pgqmx = - ( 1d0 + ( 1d0 + x )**2d0 ) / x
      S1x   = - ddilog(1d0-x)
      S2x   = S2(x)
*
      X1GQTA = 2d0 * NF * ( 4d0 * CF**2d0 * ( - 1d0 / 2d0 + 9d0 * x 
     1       / 2d0 + ( - 8d0 + x / 2d0 ) * lnx + 2d0 * x * ln1mx 
     2       + ( 1d0 - x / 2d0 ) * lnx**2d0 + ( ln1mx**2d0 
     3       + 4d0 * lnx * ln1mx - 8d0 * S1x 
     4       - 4d0 * pi**2d0 / 3d0 ) * pgq )
     5       + 4d0 * CF * CA * ( 62d0 / 9d0 - 35d0 * x / 18d0 
     6       - 44d0 * x**2d0 / 9d0 
     7       + ( 2d0 + 12d0 * x + 8d0 * x**2d0 / 3d0 ) * lnx 
     8       - 2d0 * x * ln1mx - ( 4d0 + x ) * lnx**2d0 + pgqmx * S2x 
     9       + ( - 2d0 * lnx * ln1mx - 3d0 * lnx - 3d0 * lnx**2d0 / 2d0 
     1       - ln1mx**2d0 + 8d0 * S1x + 7d0 * pi**2d0 / 6d0 
     2       + 17d0 / 18d0 ) * pgq ) )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The regular piece of the 2-loop gg time-like splitting function P_gg^(1) 
*
      FUNCTION X1GGTA (X, NF)
*
      IMPLICIT REAL*8 (A - Z)
*
      INTEGER NF
      PARAMETER ( Z2 = 1.6449 34066 84822 64365 D0,
     1            Z3 = 1.2020 56903 15959 42854 D0 )
*
* ..The soft coefficient for use in X1GGB and X1GGC
*
       COMMON / P1GSOFTT / A2GT
*
* ...Colour factors and abbreviation
*
      CF = 4D0/3D0
      CA = 3D0
      lnx = dlog(x)
      ln1mx = dlog(1d0 - x)
      pgg   = ( 1d0 / ( 1d0 - x ) +  1d0 / x - 2d0 + x * ( 1d0 - x ) )
      pggmx = ( 1d0 / ( 1d0 + x ) -  1d0 / x - 2d0 - x * ( 1d0 + x ) )
      S2x   = S2(x)
      DM    = 1D0/(1D0-X)
*
      ggg1  = 2d0 * CF * NF * ( - 4d0 + 12d0 * x - 164d0 * x**2d0 / 9d0 
     1      + ( 10d0 + 14d0 * x + 16d0 * x**2d0 / 3d0 + 16d0 / 3d0 / x )
     2      * lnx + 92d0 / 9d0 / x + 2d0 * ( 1d0 + x ) * lnx**2 )
     3      + 2d0 * CA * NF * ( 2d0 - 2d0 * x 
     4      + 26d0 * ( x**2d0 - 1d0 / x ) / 9d0 
     5      - 4d0 * ( 1d0 + x ) * lnx / 3d0 
     6      - ( 20d0 / 9d0 + 8d0 * lnx / 3d0 ) * pgg )
     7      + 4d0 * CA * CA * ( 27d0 * ( 1d0 - x ) / 2d0 
     8      + 67d0 * ( x**2d0 - 1d0 / x ) / 9d0
     9      + ( 11d0 / 3d0 - 25d0 * x / 3d0 - 44d0 / 3d0 / x ) * lnx
     1      - 4d0 * ( 1d0 + x ) * lnx**2d0
     2      + ( 4d0 * lnx * ln1mx - 3d0 * lnx**2d0 + 22d0 * lnx / 3d0
     3      - 2d0 * z2 + 67d0 / 9d0 ) * pgg + 2d0 * pggmx * S2x )
*
* ...The soft (`+'-distribution) part of the splitting function
*
       A2GT = - 40D0/9D0*ca*nf + 268D0/9D0*ca**2 - 8D0*z2*ca**2
*
       GGG1L = DM * A2GT
*
* ...The regular piece of the coefficient function
*
       X1GGTA = GGG1 - GGG1L
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the singular (soft) piece.
*
       FUNCTION X1GGTB (Y)
       IMPLICIT REAL*8 (A - Z)
*
       COMMON / P1GSOFTT / A2GT
*
       X1GGTB  = A2GT/(1D0-Y)
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' piece.
*
       FUNCTION X1GGTC (Y, NF)
*
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
       PARAMETER ( Z2 = 1.6449 34066 84822 64365 D0,
     ,             Z3 = 1.2020 56903 15959 42854 D0 )
c*
c       COMMON / P1GSOFTT / A2GT
*
* ...Colour factors
*
       CF  = 4D0/3D0
       CA  = 3D0
*
* ...The coefficient of delta(1-x)
*
       P1DELT = 
     ,    - 2D0*cf*nf
     ,    - 8D0/3D0*ca*nf
     ,    + 32D0/3D0*ca**2
     ,    + 12D0*z3*ca**2
*
* ...The soft (`+'-distribution) part of the splitting function
*
       A2GT = - 40D0/9D0*ca*nf + 268D0/9D0*ca**2 - 8D0*z2*ca**2
*
       X1GGTC = DLOG (1D0-Y) * A2GT + P1DELT
*
       RETURN
       END
*
* =====================================================================
*
*
* ..File: xpij2pt.f 
*
*                           __
* ..The parametrized 3-loop MS singlet splitting functions  P^(2)T 
*    for the evolution of unpolarized singlet fragmentation densities 
*    at mu_r = mu_f.  The expansion parameter is alpha_s/(4 pi).

* ..The distributions (in the mathematical sense) are given as in eq.
*   (B.27) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.
*   The name-endings A, B, and C of the functions below correspond to
*   the kernel superscripts [2], [3], and [1] in that equation.
*
* ..The relative accuracy of these parametrisations, as well as of
*    the convolution results, is better than one part in thousand.

* ..The coefficients of 1/(1-x)_+, (ln x)/x and 1/x are exact (up
*    to a truncation of irrational coefficients).  Furthermore all
*    coefficients written as fractions (e.g., 160./27.D0) are exact.
*    The other terms at x < 1 have fitted to the exact results for x
*    between 10^-6 and 1 - 10^-6.  The coefficient of delta(1-x) of
*    P_gg^(2) have been very slightly adjusted using the low moments.
*
* ..References: S. Moch and A. Vogt,
*               Phys. Lett. B659 (2008) 290, arXiv:0708.3899  
*               A. Almasy, S. Moch and A. Vogt,
*               arXiv:1107.nnnn                            
*                           
* =====================================================================
*
* ..The (regular) pure-singlet splitting functions P_ps^(2)T.
*    P_qq^(2)T is obtained by adding the non-singlet quantity P_NS+^(2)T
*    A parametrization of the latter is provided in the file  xpns2pt.f.
*
       FUNCTION P2PSTA (Y, NF)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       P2PST1 = - 256./(9.D0*Y)* DL**3 - 128./(9.D0*Y)* DL**2 
     &          + 324.07/Y* DL + 479.87/Y
     &          - 5.926* DL1**3 - 9.751* DL1**2 - 8.65* DL1 - 106.65
     &          - 848.97* Y + 368.79* Y**2 - 61.284* Y**3
     &          + 96.171* DL*DL1 + 656.49* DL + 425.14* DL**2 
     &          + 47.322* DL**3 + 9.072* DL**4
       P2PST2 = - 128./(81.D0*Y) + 1.778* DL1**2 + 16.611* DL1 + 87.795
     &          - 57.688* Y - 41.827* Y**2 + 25.628* Y**3 - 7.9934* Y**4
     &          - 2.1031* DL*DL1 + 57.713* DL + 9.1682* DL**2 
     &          - 1.9* DL**3 + 0.019122* DL**4
     &          + 26.294* Y*DL - 7.8645* Y*DL**3
*
       P2PSTA = (1.-Y) * NF * ( P2PST1 + NF * P2PST2 )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The quark-gluon splitting functions P_qg^(2)T.
*   The nf^3 part is exact up to a truncation of zeta_2
*
       FUNCTION P2QGTA (Y, NF)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER  NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)

*
       P2QG1 = - 64./Y* DL**3 - 64./Y* DL**2 + 675.83/Y* DL + 1141.7/Y
     &         + 100./27.D0* DL1**4 + 350./9.D0* DL1**3
     &         + 263.07* DL1**2 + 693.84* DL1 + 603.71
     &         - 882.48* Y + 4723.2* Y**2 - 4745.8* Y**3 - 175.28* Y**4
     &         + 1864.* DL + 1512.* DL**2 + 361.28* DL**3 
     &         + 42.328* DL**4 - 1809.4* DL*DL1 - 107.59* Y*DL*DL1 
     &         - 885.5* Y*DL**4
       P2QG2 = - 32./(27.D0*Y)* DL**2 - 3.1752/Y* DL - 2.8986/Y
     &         - 100./27.D0* DL1**3 - 35.446* DL1**2 - 103.609* DL1
     &         - 113.81 + 341.26* Y - 853.35* Y**2 + 492.1* Y**3 
     &         + 14.803* Y**4 + 619.75* DL + 255.62* DL**2 
     &         + 21.569* DL**3 + 966.96* DL*DL1 - 1.593*DL*DL1**2 
     &         - 333.8* Y*DL**3 - 709.1* Y*DL*DL1 
       P2QG3 = 4./9.D0* (4. + 6.* (DL + DL1)
     &         + (1. - 2.*Y + 2.*Y**2) * ( 3.8696 + 4.* (DL + DL1) 
     &         + 3.* (DL + DL1)**2 ) )
* 
       P2QGTA = NF * ( P2QG1 + NF * P2QG2 + NF**2 * P2QG3 )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The gluon-quark splitting functions P_gq^(2)T
*
       FUNCTION P2GQTA (Y, NF)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       P2GQ0 =  400./81.D0* DL1**4 + 520./27.D0* DL1**3 
     &        - 220.13* DL1**2 - 152.60* DL1 + 272.85 - 7188.7* Y 
     &        + 5693.2* Y**2 + 146.98* Y**3 + 128.19* Y**4
     &        - 30.062* DL**4 - 126.38* DL**3 - 0.71252* DL**2 
     &        + 4.4136* DL - 1300.6* DL*DL1 - 71.23* DL*DL1**2 
     &        + 543.8* Y*DL**3 
     &        + 256./Y* DL**4 + 3712./(3.D0*Y)* DL**3
     &        + 1001.89/Y* DL**2 + 4776.5/Y* DL + 5803.7/Y
       P2GQ1 =  80./81.D0* DL1**3 + 1040./81.D0* DL1**2 - 16.914* DL1
     &        - 871.3 + 790.13* Y - 241.23* Y**2 + 43.252* Y**3
     &        - 48.600 * DL**3 - 343.1* DL**2 - 492.* DL
     &        + 55.048* DL*DL1 - 4.3465* Y*DL**3
     &        + 6.0041/Y + 141.93/Y* DL 
     &        + 2912./(27.D0*Y)* DL**2 + 1280./(81.D0*Y)*DL**3
*
       P2GQTA = ( P2GQ0 + NF * P2GQ1 )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The regular piece of the gluon-gluon splitting function P_gg^(2)
*
       FUNCTION P2GGTA (Y, NF)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       P2GGTA0 =   576./Y* DL**4 + 3168./Y* DL**3 + 3651.1/Y* DL**2 
     &           + 10233./Y* DL + 14214.4/Y - 3590.1* DL1 - 28489.
     &           + 7469.* Y + 30421.* Y**2 - 53017.* Y**3 + 19556.* Y**4
     &           + 191.99* DL**4 + 3281.7* DL**3 + 13528.* DL**2 
     &           + 12258.* DL - 186.4* DL*DL1 - 21328.* DL**2*DL1 
     &           + 5685.8* Y*DL**3
       P2GGTA1 = + 448./(9.D0*Y)* DL**3 + 2368./(9.D0*Y)* DL**2 
     &           - 5.47/Y* DL - 804.13/Y + 248.95 + 319.97* DL1
     &           + 260.6* Y + 272.79* Y**2 + 2133.2* Y**3 - 926.87* Y**4
     &           + 4.9934* DL + 482.94* DL**2 + 155.10* DL**3 
     &           + 18.085* DL**4 + 485.18* Y*DL**3 + 1266.5* DL*DL1 
     &           - 29.709* DL**2*DL1 + 87.771* DL*DL1**2
       P2GGTA2 =   32./(27.D0*Y)* DL**2 + 368./(81.D0*Y)* DL 
     &           + 472./(243.D0*Y) 
     &           - 77.190 + 153.27* Y - 106.03* Y**2 + 11.995* Y**3
     &           - 5.0372* DL**3 - 44.8* DL**2 - 69.712* DL
     &           - 115.01* DL*DL1 + 96.522* Y*DL*DL1 - 62.908* DL**2*DL1
*
       P2GGTA = P2GGTA0 + NF * ( P2GGTA1 + NF * P2GGTA2 )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The singular piece of the gluon-gluon splitting function P_gg^(2)T
*    (identical to the spacelike case)
*
       FUNCTION P2GGTB (Y, NF)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       P2GGTB = (2643.521 - NF * 412.172 - NF**2 * 16./9.D0) / ( 1.-Y)
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The 'local' piece of the gluon-gluon splitting function P_gg^(2)
*   (as in the spacelike case, up to the smaller delta(1-x) adjustment)
*
       FUNCTION P2GGTC (Y, NF)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL1 = LOG (1.-Y)
*
       P2GGTC =       2643.521 * DL1 + 4425.448 + 0.003
     ,       - NF * ( 412.172 * DL1 +  528.720 - 0.001 )
     ,       + NF**2 * ( - 16./9.D0 * DL1 + 6.4630 - 0.0002)
*
       RETURN
       END
*
*
* =================================================================aaa=
* ..File: xpns2tp.f 
*
*                           __
* ..The parametrized 3-loop MS non-singlet splitting functions P^(2)T
*    for the evolution of unpolarized fragmentation densities at 
*    mu_r = mu_f.  The expansion parameter is alpha_s/(4 pi).
*
* ..The distributions (in the mathematical sense) are given as in eq.
*    (B.26) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.
*    The name-endings A, B, and C of the functions below correspond to
*    the kernel superscripts [2], [3], and [1] in that equation.
*
*  ..The relative accuracy of these parametrizations, as well as of
*    the convolution results, is better than one part in thousand.
*
* ..References: A. Mitov, S. Moch and A. Vogt,
*               Phys. Lett. B638 (2006) 61, hep-ph/0604053 (exact res.)
*               A. Almasy, S. Moch and A. Vogt,
*               arXiv:1107.nnnn                            (this code)
*
* ..Some (parts) of these functions are identical to the spacelike
*    results of S. Moch, J. Vermaseren and A. Vogt,
*               Nucl. Phys. B688 (2004) 101, hep-ph/040319
*      
* =====================================================================
*
*
* ..This is the regular piece of P2_NS+.  The rational coefficients are
*    exact, the rest has been fitted for x between 10^-6 and 1 - 10^-6.
*    The N_f^2 part is exact is identical to the spacelike case.
*
       FUNCTION P2NSPTA (Y, NF)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
       D81 = 1./81.D0
*
       P2NSPTA =   1658.7 - 707.67* DL1 + 1327.5* DL - 56.907* DL*DL1 
     &           - 189.37* DL**2 - 519.37* DL1*DL**2 - 352./9.D0* DL**3 
     &           + 128./81.D0* DL**4 - 4249.4* Y - 559.1* DL1*DL*Y 
     &           - 1075.3* Y**2 + 593.9* Y**3
     &      + NF * (64./27.D0* DL**3 - 176./81.D0* DL**2 - 168.89* DL
     &           - 198.10 + 466.29* Y + 181.18* Y**2 - 31.84* Y**3
     &           + 5120./81.D0* DL1 - 50.758* DL*DL1 + 28.551* DL**2*DL1
     &           - 39.113* Y*DL + 85.72* Y*DL*DL1 - 23.102* Y*DL**2*DL1)
     &      + NF**2 * ( 32.* Y*DL/(1.-Y) * (3.* DL + 10.) + 64.
     &           + (48.* DL**2 + 352.* DL + 384.) * (1.-Y) ) * D81
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the regular piece of P2_NS-.  The rational coefficients are 
*    exact, the rest has been fitted for x between 10^-6 and 1 - 10^-6.
*    The N_f^2 part is exact (and identical to that of P2_NS+). 
*
       FUNCTION P2NSMTA (Y, NF)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
       D81 = 1./81.D0
*
       P2NSMTA =  - 140./81.D0* DL**4 - 1024./27.D0* DL**3 
     &            - 38.298* DL**2 + 1625.5* DL - 707.94* DL1 + 1981.3
     &            - 4885.7* Y - 577.42* Y**2 + 407.89* Y**3
     &            + 1905.4* DL**2*DL1 + 1969.5* Y*DL**2*DL1 
     &            + 4563.2* DL*DL1 - 34.683* Y*DL**4 
     &            - 5140.6* Y*DL*DL1 - 437.03* Y*DL**3
     &        + NF * ( 128./81.D0* DL**3 - 784./81.D0* DL**2 
     &            - 188.99* DL - 217.84 + 511.92* Y + 209.19* Y**2 
     &            - 85.786* Y**3 + 5120./81.D0* DL1 + 71.428* DL*DL1 
     &            + 30.554* DL**2*DL1 + 92.453* Y*DL - 23.722* Y*DL*DL1 
     &            - 18.975* Y*DL**2*DL1 )
     ,        + NF**2 * ( 32.* Y*DL/(1.-Y) * (3.* DL + 10.) + 64.
     ,            + (48.* DL**2 + 352.* DL + 384.) * (1.-Y) ) * D81
*
       RETURN
       END
*
*
* ---------------------------------------------------------------------
*
*
* ..This is the singular piece of both P2_NS+ and P2_NS-. It is exact
*    up to the truncation of the irrational coefficients, and identical
*    to the spacelike result.
*
       FUNCTION P2NSTB (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       P2NSTB = ( 1174.898 - NF * 183.187 - NF**2 * 64./81.D0 ) / (1.-Y)
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' piece of P2_NS+. Coefficients of delta(1-x) have 
*    been shifted (minimally) relative to the exact (truncated) values.
*
       FUNCTION P2NSPTC (Y, NF)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
*
       DL1 = LOG (1.-Y)
*
       P2NSPTC =       1174.898 * DL1 + 1295.624 + 0.001
     ,         - NF * ( 183.187 * DL1 + 173.938 - 0.003)
     ,         + NF**2 * ( - 64./81.D0 * DL1 + 1.13067 )
*
       RETURN
       END
*
*
* ---------------------------------------------------------------------
*
* ..This is the 'local' piece of P2_NS-. Coefficients of delta(1-x) have
*    been shifted (minimally) relative to the exact (truncated) values.
*
       FUNCTION P2NSMTC (Y, NF)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
*
       DL1 = LOG (1.-Y)
*
       P2NSMTC =       1174.898 * DL1 + 1295.624 - 0.002
     ,        - NF * ( 183.187 * DL1 + 173.938 - 0.0004 )
     ,        + NF**2 * ( - 64./81.D0 * DL1 + 1.13067 )
*
       RETURN
       END
*
*
* ---------------------------------------------------------------------
*
*
* ..This is P2_NSS, the difference of P2_NSV and P2_NS-. It is identical
*    to the spacelike result.
*
       FUNCTION P2NSSTA (Y, NF)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       D27 = 1./27.D0
       DL  = LOG (Y)
       Y1  = 1.- Y
       DL1 = LOG (Y1)
*
       P2NSSA = Y1* ( 151.49 + 44.51 * Y - 43.12 * Y**2 + 4.820 * Y**3 )
     1          + 40.*D27 * DL**4 - 80.*D27 * DL**3 + 6.892 * DL**2
     2          + 178.04 * DL + DL*DL1 * ( - 173.1 + 46.18 * DL )
     4          + Y1*DL1 * ( - 163.9 / Y - 7.208 * Y )
*
       P2NSSTA  = NF * P2NSSA
*
       RETURN
       END
*
* =================================================================aaa=
c$$$*
c$$$* ..File: xpij2te.f 
c$$$*
c$$$*
c$$$* ..The exact 3-loop MS(bar) singlet splitting functions  P_ij^(2)T
c$$$*    for the evolution of unpolarized fragmentation densities at
c$$$*    mu_r = mu_f.  The expansion parameter is alpha_s/(4 pi).
c$$$*
c$$$* ..The function  X2QGAT  includes a parameter IMOD for estimating, 
c$$$*    via  IMOD = 1 and IMOD = -1, a residual uncertainty discussed in 
c$$$*    the paper.  IMOD = 0  provides the main result.
c$$$* 
c$$$* ..The distributions (in the mathematical sense) are given as in eq.
c$$$*    (B.26) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.
c$$$*    The name-endings A, B, and C of the functions below correspond to 
c$$$*    the kernel superscripts [2], [3], and [1] in that equation.
c$$$*
c$$$* ..The code uses the package of Gehrmann and Remiddi for the harmonic
c$$$*    polylogarithms published in hep-ph/0107173 = CPC 141 (2001) 296.
c$$$*
c$$$* ..References: S. Moch and A. Vogt, 
c$$$*               Phys. Lett. B659 (2008) 290, arXiv:0708.3899   (ps, gg)
c$$$*               A. Almasy, S. Moch and A. Vogt,
c$$$*               arXiv:1107.nnnn                                (qg, gq)
c$$$*
c$$$* =====================================================================
c$$$*
c$$$*
c$$$* ..The (regular) pure-singlet splitting functions P_ps^(2)T.
c$$$*    P_qq^(2) is obtained by adding the non-singlet quantity P_NS^(2)T+
c$$$*
c$$$       FUNCTION X2PSTA (X, NF)
c$$$*
c$$$       IMPLICIT REAL*8 (A - Z)
c$$$       COMPLEX*16 HC1, HC2, HC3, HC4
c$$$       INTEGER NF, NF2, N1, N2, NW, I1, I2, I3, N, IMOD
c$$$       PARAMETER ( N1 = -1, N2 = 1, NW = 4 )
c$$$       DIMENSION HC1(N1:N2),HC2(N1:N2,N1:N2),HC3(N1:N2,N1:N2,N1:N2),
c$$$     ,           HC4(N1:N2,N1:N2,N1:N2,N1:N2)
c$$$       DIMENSION HR1(N1:N2),HR2(N1:N2,N1:N2),HR3(N1:N2,N1:N2,N1:N2),
c$$$     ,           HR4(N1:N2,N1:N2,N1:N2,N1:N2)
c$$$       DIMENSION HI1(N1:N2),HI2(N1:N2,N1:N2),HI3(N1:N2,N1:N2,N1:N2),
c$$$     ,           HI4(N1:N2,N1:N2,N1:N2,N1:N2)
c$$$       PARAMETER ( Z2 = 1.6449 34066 84822 64365 D0,
c$$$     ,             Z3 = 1.2020 56903 15959 42854 D0 )
c$$$*
c$$$* ...Colour factors and an abbreviation
c$$$*
c$$$       CF  = 4./3.D0
c$$$       CA  = 3.D0
c$$$       NF2 = NF*NF
c$$$*
c$$$       DX = 1.D0/X
c$$$*
c$$$* ...The harmonic polylogs up to weight 4 by Gehrmann and Remiddi
c$$$*
c$$$       CALL HPLOG (X, NW, HC1,HC2,HC3,HC4, HR1,HR2,HR3,HR4,
c$$$     ,            HI1,HI2,HI3,HI4, N1, N2)
c$$$*
c$$$* ...The splitting function in terms of the harmonic polylogs
c$$$*
c$$$      PTqqps2 =
c$$$     &  + nf*cf*ca * ( 820.D0/3.D0 - 1172.D0/3.D0*x + 7496.D0/81.D0*
c$$$     &    x**2 + 2008.D0/81.D0*dx + 296.D0/3.D0*z3 - 1312.D0/3.D0*z3*x
c$$$     &     + 256.D0/3.D0*z3*x**2 + 2212.D0/9.D0*z2 + 4912.D0/9.D0*z2*x
c$$$     &     + 640.D0/9.D0*z2*x**2 + 512.D0/9.D0*z2*dx - 852.D0/5.D0*
c$$$     &    z2**2 - 812.D0/5.D0*z2**2*x + 40.D0*Hr1(-1)*z2 + 40.D0*Hr1(-1
c$$$     &    )*z2*x + 32.D0/3.D0*Hr1(-1)*z2*x**2 + 32.D0/3.D0*Hr1(-1)*z2*
c$$$     &    dx - 2380.D0/27.D0*Hr1(0) + 23540.D0/27.D0*Hr1(0)*x + 9172.D0/
c$$$     &    27.D0*Hr1(0)*x**2 + 1240.D0/27.D0*Hr1(0)*dx - 160.D0*Hr1(0)*
c$$$     &    z3 - 32.D0*Hr1(0)*z3*x + 196.D0/3.D0*Hr1(0)*z2 - 404.D0/3.D0*
c$$$     &    Hr1(0)*z2*x + 32.D0/3.D0*Hr1(0)*z2*x**2 + 64.D0/3.D0*Hr1(0)*
c$$$     &    z2*dx - 2596.D0/9.D0*Hr1(1) + 592.D0/9.D0*Hr1(1)*x + 7868.D0/
c$$$     &    27.D0*Hr1(1)*x**2 - 1856.D0/27.D0*Hr1(1)*dx - 72.D0*Hr1(1)*z2
c$$$     &     + 72.D0*Hr1(1)*z2*x + 96.D0*Hr1(1)*z2*x**2 - 96.D0*Hr1(1)*z2
c$$$     &    *dx + 688.D0/3.D0*Hr2(-1,0) + 400.D0/3.D0*Hr2(-1,0)*x - 992.D0
c$$$     &    /9.D0*Hr2(-1,0)*x**2 - 128.D0/9.D0*Hr2(-1,0)*dx + 16.D0*Hr2(0
c$$$     &    ,-1)*z2 )
c$$$      PTqqps2 = PTqqps2 + nf*cf*ca * (  - 16.D0*Hr2(0,-1)*z2*x - 2440.D0
c$$$     &    /9.D0*Hr2(0,0) - 10468.D0/9.D0*Hr2(0,0)*x - 2312.D0/9.D0*Hr2(
c$$$     &    0,0)*x**2 - 64.D0/9.D0*Hr2(0,0)*dx - 8.D0*Hr2(0,0)*z2 + 104.D0
c$$$     &    *Hr2(0,0)*z2*x - 1540.D0/9.D0*Hr2(0,1) - 3520.D0/9.D0*Hr2(0,1
c$$$     &    )*x - 2032.D0/9.D0*Hr2(0,1)*x**2 - 112.D0/9.D0*Hr2(0,1)*dx - 
c$$$     &    144.D0*Hr2(0,1)*z2 - 144.D0*Hr2(0,1)*z2*x + 4.D0/3.D0*Hr2(1,0
c$$$     &    ) + 140.D0/3.D0*Hr2(1,0)*x - 72.D0*Hr2(1,0)*x**2 + 24.D0*Hr2(
c$$$     &    1,0)*dx - 12.D0*Hr2(1,1) + 12.D0*Hr2(1,1)*x + 16.D0/9.D0*Hr2(
c$$$     &    1,1)*x**2 - 16.D0/9.D0*Hr2(1,1)*dx + 16.D0*Hr3(-1,-1,0) + 16.D
c$$$     &    0*Hr3(-1,-1,0)*x - 64.D0/3.D0*Hr3(-1,-1,0)*x**2 - 64.D0/3.D0*
c$$$     &    Hr3(-1,-1,0)*dx - 56.D0*Hr3(-1,0,0) - 56.D0*Hr3(-1,0,0)*x + 
c$$$     &    128.D0/3.D0*Hr3(-1,0,0)*x**2 + 128.D0/3.D0*Hr3(-1,0,0)*dx - 
c$$$     &    32.D0*Hr3(-1,0,1) - 32.D0*Hr3(-1,0,1)*x - 64.D0/3.D0*Hr3(-1,0
c$$$     &    ,1)*x**2 - 64.D0/3.D0*Hr3(-1,0,1)*dx + 120.D0*Hr3(0,-1,0) - 
c$$$     &    232.D0*Hr3(0,-1,0)*x + 256.D0/3.D0*Hr3(0,-1,0)*x**2 + 64.D0/3.
c$$$     &    D0*Hr3(0,-1,0)*dx )
c$$$      PTqqps2 = PTqqps2 + nf*cf*ca * (  - 112.D0/3.D0*Hr3(0,0,0) + 2288.
c$$$     &    D0/3.D0*Hr3(0,0,0)*x - 128.D0/3.D0*Hr3(0,0,0)*dx + 236.D0/3.D0
c$$$     &    *Hr3(0,0,1) + 380.D0/3.D0*Hr3(0,0,1)*x + 224.D0/3.D0*Hr3(0,0,
c$$$     &    1)*x**2 + 64.D0/3.D0*Hr3(0,0,1)*dx - 64.D0*Hr3(0,1,0) + 24.D0
c$$$     &    *Hr3(0,1,0)*x + 224.D0/3.D0*Hr3(0,1,0)*x**2 - 64.D0/3.D0*Hr3(
c$$$     &    0,1,0)*dx - 136.D0/3.D0*Hr3(0,1,1) - 88.D0/3.D0*Hr3(0,1,1)*x
c$$$     &     + 32.D0/3.D0*Hr3(0,1,1)*x**2 - 32.D0/3.D0*Hr3(0,1,1)*dx + 60.
c$$$     &    D0*Hr3(1,0,0) - 60.D0*Hr3(1,0,0)*x - 96.D0*Hr3(1,0,0)*x**2 + 
c$$$     &    96.D0*Hr3(1,0,0)*dx + 48.D0*Hr3(1,0,1) - 48.D0*Hr3(1,0,1)*x
c$$$     &     - 64.D0*Hr3(1,0,1)*x**2 + 64.D0*Hr3(1,0,1)*dx + 32.D0*Hr3(1,
c$$$     &    1,0) - 32.D0*Hr3(1,1,0)*x - 128.D0/3.D0*Hr3(1,1,0)*x**2 + 128.
c$$$     &    D0/3.D0*Hr3(1,1,0)*dx + 16.D0*Hr3(1,1,1) - 16.D0*Hr3(1,1,1)*x
c$$$     &     - 64.D0/3.D0*Hr3(1,1,1)*x**2 + 64.D0/3.D0*Hr3(1,1,1)*dx + 32.
c$$$     &    D0*Hr4(0,-1,-1,0) - 32.D0*Hr4(0,-1,-1,0)*x - 80.D0*Hr4(0,-1,0
c$$$     &    ,0) + 80.D0*Hr4(0,-1,0,0)*x - 48.D0*Hr4(0,0,-1,0) + 112.D0*
c$$$     &    Hr4(0,0,-1,0)*x )
c$$$      PTqqps2 = PTqqps2 + nf*cf*ca * ( 64.D0*Hr4(0,0,0,0) - 304.D0*Hr4(
c$$$     &    0,0,0,0)*x + 8.D0*Hr4(0,0,0,1) + 8.D0*Hr4(0,0,0,1)*x - 144.D0
c$$$     &    *Hr4(0,0,1,0) - 144.D0*Hr4(0,0,1,0)*x - 48.D0*Hr4(0,0,1,1) - 
c$$$     &    48.D0*Hr4(0,0,1,1)*x + 136.D0*Hr4(0,1,0,0) + 136.D0*Hr4(0,1,0
c$$$     &    ,0)*x + 96.D0*Hr4(0,1,0,1) + 96.D0*Hr4(0,1,0,1)*x + 64.D0*
c$$$     &    Hr4(0,1,1,0) + 64.D0*Hr4(0,1,1,0)*x + 32.D0*Hr4(0,1,1,1) + 32.
c$$$     &    D0*Hr4(0,1,1,1)*x )
c$$$      PTqqps2 = PTqqps2 + nf*cf**2 * (  - 74.D0/27.D0 + 6266.D0/27.D0*x
c$$$     &     - 156.D0*x**2 - 220.D0/3.D0*dx + 288.D0*z3 + 352.D0*z3*x + 
c$$$     &    256.D0/3.D0*z3*x**2 + 64.D0*z3*dx - 200.D0/3.D0*z2 - 1444.D0/
c$$$     &    3.D0*z2*x - 2624.D0/9.D0*z2*x**2 + 984.D0/5.D0*z2**2 + 984.D0/
c$$$     &    5.D0*z2**2*x - 598.D0/9.D0*Hr1(0) - 4406.D0/9.D0*Hr1(0)*x - 
c$$$     &    9836.D0/27.D0*Hr1(0)*x**2 + 160.D0*Hr1(0)*z3 + 160.D0*Hr1(0)*
c$$$     &    z3*x + 80.D0*Hr1(0)*z2*x + 128.D0/3.D0*Hr1(0)*z2*x**2 + 856.D0
c$$$     &    /3.D0*Hr1(1) - 140.D0/3.D0*Hr1(1)*x - 8972.D0/27.D0*Hr1(1)*
c$$$     &    x**2 + 2528.D0/27.D0*Hr1(1)*dx + 80.D0*Hr1(1)*z2 - 80.D0*Hr1(
c$$$     &    1)*z2*x - 320.D0/3.D0*Hr1(1)*z2*x**2 + 320.D0/3.D0*Hr1(1)*z2*
c$$$     &    dx + 976.D0/3.D0*Hr2(0,0) + 1184.D0/3.D0*Hr2(0,0)*x + 608.D0/
c$$$     &    9.D0*Hr2(0,0)*x**2 - 16.D0*Hr2(0,0)*z2 - 16.D0*Hr2(0,0)*z2*x
c$$$     &     + 296.D0*Hr2(0,1) + 348.D0*Hr2(0,1)*x + 512.D0/3.D0*Hr2(0,1)
c$$$     &    *x**2 + 224.D0/9.D0*Hr2(0,1)*dx + 160.D0*Hr2(0,1)*z2 + 160.D0
c$$$     &    *Hr2(0,1)*z2*x + 428.D0/3.D0*Hr2(1,0) - 284.D0/3.D0*Hr2(1,0)*
c$$$     &    x )
c$$$      PTqqps2 = PTqqps2 + nf*cf**2 * (  - 1120.D0/9.D0*Hr2(1,0)*x**2 + 
c$$$     &    688.D0/9.D0*Hr2(1,0)*dx + 28.D0/3.D0*Hr2(1,1) - 28.D0/3.D0*
c$$$     &    Hr2(1,1)*x - 64.D0/3.D0*Hr2(1,1)*x**2 + 64.D0/3.D0*Hr2(1,1)*
c$$$     &    dx + 72.D0*Hr3(0,0,0) + 184.D0*Hr3(0,0,0)*x + 160.D0/3.D0*
c$$$     &    Hr3(0,0,0)*x**2 + 32.D0*Hr3(0,0,1) + 80.D0*Hr3(0,0,1)*x - 128.
c$$$     &    D0/3.D0*Hr3(0,0,1)*dx + 208.D0*Hr3(0,1,0) + 224.D0*Hr3(0,1,0)
c$$$     &    *x + 64.D0/3.D0*Hr3(0,1,0)*x**2 + 128.D0/3.D0*Hr3(0,1,0)*dx
c$$$     &     + 64.D0*Hr3(0,1,1) + 48.D0*Hr3(0,1,1)*x + 64.D0/3.D0*Hr3(0,1
c$$$     &    ,1)*dx - 88.D0*Hr3(1,0,0) + 88.D0*Hr3(1,0,0)*x + 352.D0/3.D0*
c$$$     &    Hr3(1,0,0)*x**2 - 352.D0/3.D0*Hr3(1,0,0)*dx - 48.D0*Hr3(1,0,1
c$$$     &    ) + 48.D0*Hr3(1,0,1)*x + 64.D0*Hr3(1,0,1)*x**2 - 64.D0*Hr3(1,
c$$$     &    0,1)*dx - 32.D0*Hr3(1,1,0) + 32.D0*Hr3(1,1,0)*x + 128.D0/3.D0
c$$$     &    *Hr3(1,1,0)*x**2 - 128.D0/3.D0*Hr3(1,1,0)*dx - 16.D0*Hr3(1,1,
c$$$     &    1) + 16.D0*Hr3(1,1,1)*x + 64.D0/3.D0*Hr3(1,1,1)*x**2 - 64.D0/
c$$$     &    3.D0*Hr3(1,1,1)*dx - 64.D0*Hr4(0,0,0,0) - 64.D0*Hr4(0,0,0,0)*
c$$$     &    x )
c$$$      PTqqps2 = PTqqps2 + nf*cf**2 * (  - 112.D0*Hr4(0,0,0,1) - 112.D0*
c$$$     &    Hr4(0,0,0,1)*x + 80.D0*Hr4(0,0,1,0) + 80.D0*Hr4(0,0,1,0)*x + 
c$$$     &    48.D0*Hr4(0,0,1,1) + 48.D0*Hr4(0,0,1,1)*x - 176.D0*Hr4(0,1,0,
c$$$     &    0) - 176.D0*Hr4(0,1,0,0)*x - 96.D0*Hr4(0,1,0,1) - 96.D0*Hr4(0
c$$$     &    ,1,0,1)*x - 64.D0*Hr4(0,1,1,0) - 64.D0*Hr4(0,1,1,0)*x - 32.D0
c$$$     &    *Hr4(0,1,1,1) - 32.D0*Hr4(0,1,1,1)*x )
c$$$      PTqqps2 = PTqqps2 + nf2*cf * (  - 344.D0/27.D0 - 376.D0/27.D0*x
c$$$     &     + 752.D0/27.D0*x**2 - 32.D0/27.D0*dx + 64.D0/3.D0*z3 + 64.D0/
c$$$     &    3.D0*z3*x + 296.D0/9.D0*z2 + 440.D0/9.D0*z2*x + 32.D0/3.D0*z2
c$$$     &    *x**2 + 736.D0/27.D0*Hr1(0) + 304.D0/27.D0*Hr1(0)*x - 16.D0*
c$$$     &    Hr1(0)*x**2 - 16.D0/3.D0*Hr1(0)*z2 - 16.D0/3.D0*Hr1(0)*z2*x
c$$$     &     - 128.D0/3.D0*Hr1(1) + 80.D0/3.D0*Hr1(1)*x + 752.D0/27.D0*
c$$$     &    Hr1(1)*x**2 - 320.D0/27.D0*Hr1(1)*dx + 64.D0/9.D0*Hr2(0,0) + 
c$$$     &    208.D0/9.D0*Hr2(0,0)*x + 32.D0/3.D0*Hr2(0,0)*x**2 - 296.D0/9.D
c$$$     &    0*Hr2(0,1) - 440.D0/9.D0*Hr2(0,1)*x - 32.D0/3.D0*Hr2(0,1)*
c$$$     &    x**2 + 8.D0*Hr2(1,0) - 8.D0*Hr2(1,0)*x - 32.D0/3.D0*Hr2(1,0)*
c$$$     &    x**2 + 32.D0/3.D0*Hr2(1,0)*dx + 8.D0/3.D0*Hr2(1,1) - 8.D0/3.D0
c$$$     &    *Hr2(1,1)*x - 32.D0/9.D0*Hr2(1,1)*x**2 + 32.D0/9.D0*Hr2(1,1)*
c$$$     &    dx - 32.D0/3.D0*Hr3(0,0,0) - 32.D0/3.D0*Hr3(0,0,0)*x + 16.D0/
c$$$     &    3.D0*Hr3(0,0,1) + 16.D0/3.D0*Hr3(0,0,1)*x + 16.D0*Hr3(0,1,0)
c$$$     &     + 16.D0*Hr3(0,1,0)*x + 16.D0/3.D0*Hr3(0,1,1) + 16.D0/3.D0*
c$$$     &    Hr3(0,1,1)*x )
c$$$
c$$$*
c$$$       X2PSTA = PTqqps2
c$$$*
c$$$       RETURN
c$$$       END
c$$$*
c$$$* ---------------------------------------------------------------------
c$$$*
c$$$*
c$$$* ..The quark-gluon splitting functions P_qg^(2)T
c$$$*
c$$$       FUNCTION X2QGTA (X, NF, IMOD)
c$$$*
c$$$       IMPLICIT REAL*8 (A - Z)
c$$$       COMPLEX*16 HC1, HC2, HC3, HC4 
c$$$       INTEGER NF, NF2, N1, N2, NW, I1, I2, I3, N, CORR, IMOD
c$$$       PARAMETER ( N1 = -1, N2 = 1, NW = 4 ) 
c$$$       DIMENSION HC1(N1:N2),HC2(N1:N2,N1:N2),HC3(N1:N2,N1:N2,N1:N2), 
c$$$     ,           HC4(N1:N2,N1:N2,N1:N2,N1:N2) 
c$$$       DIMENSION HR1(N1:N2),HR2(N1:N2,N1:N2),HR3(N1:N2,N1:N2,N1:N2), 
c$$$     ,           HR4(N1:N2,N1:N2,N1:N2,N1:N2) 
c$$$       DIMENSION HI1(N1:N2),HI2(N1:N2,N1:N2),HI3(N1:N2,N1:N2,N1:N2), 
c$$$     ,           HI4(N1:N2,N1:N2,N1:N2,N1:N2) 
c$$$       PARAMETER ( Z2 = 1.6449 34066 84822 64365 D0,
c$$$     ,             Z3 = 1.2020 56903 15959 42854 D0 )
c$$$*
c$$$* ...Colour factors and an abbreviation
c$$$*
c$$$       CF  = 4./3.D0
c$$$       CA  = 3.D0
c$$$       NF2 = NF*NF
c$$$*
c$$$       DX = 1.D0/X
c$$$*
c$$$* ...The harmonic polylogs up to weight 4 by Gehrmann and Remiddi
c$$$*
c$$$       CALL HPLOG (X, NW, HC1,HC2,HC3,HC4, HR1,HR2,HR3,HR4,
c$$$     ,            HI1,HI2,HI3,HI4, N1, N2) 
c$$$*
c$$$* ...The splitting function in terms of the harmonic polylogs
c$$$*
c$$$      PTqg2 =
c$$$     &  + nf*ca**2 * ( 2057.D0/9.D0 - 672.D0*x + 53753.D0/81.D0*x**2 + 
c$$$     &    2092.D0/81.D0*dx - 196.D0*z3 - 2416.D0/3.D0*z3*x + 328.D0/3.D0
c$$$     &    *z3*x**2 - 310.D0/3.D0*z2 + 308.D0*z2*x - 536.D0/3.D0*z2*x**2
c$$$     &     + 544.D0/9.D0*z2*dx + 34.D0*z2**2 + 252.D0/5.D0*z2**2*x + 
c$$$     &    704.D0/5.D0*z2**2*x**2 + 136.D0*Hr1(-1)*z3 + 272.D0*Hr1(-1)*
c$$$     &    z3*x + 272.D0*Hr1(-1)*z3*x**2 + 512.D0/3.D0*Hr1(-1)*z2 + 760.D
c$$$     &    0/3.D0*Hr1(-1)*z2*x - 24.D0*Hr1(-1)*z2*x**2 + 32.D0/3.D0*Hr1(
c$$$     &    -1)*z2*dx - 2996.D0/27.D0*Hr1(0) - 4832.D0/27.D0*Hr1(0)*x + 
c$$$     &    5044.D0/9.D0*Hr1(0)*x**2 + 40.D0*Hr1(0)*dx + 16.D0*Hr1(0)*z3
c$$$     &     + 960.D0*Hr1(0)*z3*x + 96.D0*Hr1(0)*z2 - 484.D0/3.D0*Hr1(0)*
c$$$     &    z2*x + 824.D0/3.D0*Hr1(0)*z2*x**2 + 64.D0/3.D0*Hr1(0)*z2*dx
c$$$     &     - 3940.D0/27.D0*Hr1(1) - 8416.D0/27.D0*Hr1(1)*x + 7448.D0/27.
c$$$     &    D0*Hr1(1)*x**2 + 752.D0/27.D0*Hr1(1)*dx + 248.D0*Hr1(1)*z3 - 
c$$$     &    496.D0*Hr1(1)*z3*x + 496.D0*Hr1(1)*z3*x**2 + 440.D0/3.D0*Hr1(
c$$$     &    1)*z2 - 664.D0/3.D0*Hr1(1)*z2*x + 464.D0/3.D0*Hr1(1)*z2*x**2
c$$$     &     - 64.D0/3.D0*Hr1(1)*z2*dx )
c$$$      PTqg2 = PTqg2 + nf*ca**2 * (  - 208.D0*Hr2(-1,-1)*z2 - 416.D0*
c$$$     &    Hr2(-1,-1)*z2*x - 416.D0*Hr2(-1,-1)*z2*x**2 - 152.D0/3.D0*
c$$$     &    Hr2(-1,0) - 256.D0/3.D0*Hr2(-1,0)*x - 248.D0*Hr2(-1,0)*x**2
c$$$     &     - 32.D0/3.D0*Hr2(-1,0)*dx + 112.D0*Hr2(-1,0)*z2 + 224.D0*
c$$$     &    Hr2(-1,0)*z2*x + 224.D0*Hr2(-1,0)*z2*x**2 + 72.D0*Hr2(0,-1)*
c$$$     &    z2 - 176.D0*Hr2(0,-1)*z2*x + 96.D0*Hr2(0,-1)*z2*x**2 - 598.D0/
c$$$     &    9.D0*Hr2(0,0) - 3508.D0/9.D0*Hr2(0,0)*x - 7208.D0/9.D0*Hr2(0,
c$$$     &    0)*x**2 - 128.D0/9.D0*Hr2(0,0)*dx - 28.D0*Hr2(0,0)*z2 + 632.D0
c$$$     &    *Hr2(0,0)*z2*x - 96.D0*Hr2(0,0)*z2*x**2 - 1262.D0/9.D0*Hr2(0,
c$$$     &    1) + 4828.D0/9.D0*Hr2(0,1)*x - 2900.D0/3.D0*Hr2(0,1)*x**2 + 
c$$$     &    112.D0/9.D0*Hr2(0,1)*dx + 168.D0*Hr2(0,1)*z2 - 432.D0*Hr2(0,1
c$$$     &    )*z2*x + 384.D0*Hr2(0,1)*z2*x**2 - 562.D0/9.D0*Hr2(1,0) + 92.D
c$$$     &    0/9.D0*Hr2(1,0)*x + 704.D0/3.D0*Hr2(1,0)*x**2 - 8.D0/9.D0*
c$$$     &    Hr2(1,0)*dx + 48.D0*Hr2(1,0)*z2 - 96.D0*Hr2(1,0)*z2*x + 96.D0
c$$$     &    *Hr2(1,0)*z2*x**2 + 134.D0/9.D0*Hr2(1,1) - 700.D0/9.D0*Hr2(1,
c$$$     &    1)*x )
c$$$      PTqg2 = PTqg2 + nf*ca**2 * ( 652.D0/3.D0*Hr2(1,1)*x**2 - 128.D0/9.
c$$$     &    D0*Hr2(1,1)*dx - 16.D0*Hr2(1,1)*z2 + 32.D0*Hr2(1,1)*z2*x - 32.
c$$$     &    D0*Hr2(1,1)*z2*x**2 + 320.D0/3.D0*Hr3(-1,-1,0) + 112.D0/3.D0*
c$$$     &    Hr3(-1,-1,0)*x - 208.D0*Hr3(-1,-1,0)*x**2 - 64.D0/3.D0*Hr3(-1
c$$$     &    ,-1,0)*dx - 872.D0/3.D0*Hr3(-1,0,0) - 520.D0/3.D0*Hr3(-1,0,0)
c$$$     &    *x + 304.D0/3.D0*Hr3(-1,0,0)*x**2 + 128.D0/3.D0*Hr3(-1,0,0)*
c$$$     &    dx - 352.D0/3.D0*Hr3(-1,0,1) - 704.D0/3.D0*Hr3(-1,0,1)*x - 80.
c$$$     &    D0*Hr3(-1,0,1)*x**2 - 64.D0/3.D0*Hr3(-1,0,1)*dx + 184.D0/3.D0
c$$$     &    *Hr3(0,-1,0) - 344.D0/3.D0*Hr3(0,-1,0)*x + 832.D0/3.D0*Hr3(0,
c$$$     &    -1,0)*x**2 + 64.D0/3.D0*Hr3(0,-1,0)*dx + 64.D0/3.D0*Hr3(0,0,0
c$$$     &    ) + 2680.D0/3.D0*Hr3(0,0,0)*x - 192.D0*Hr3(0,0,0)*x**2 - 128.D
c$$$     &    0/3.D0*Hr3(0,0,0)*dx + 64.D0/3.D0*Hr3(0,0,1) + 292.D0*Hr3(0,0
c$$$     &    ,1)*x + 120.D0*Hr3(0,0,1)*x**2 - 64.D0/3.D0*Hr3(0,0,1)*dx - 
c$$$     &    100.D0*Hr3(0,1,0) - 312.D0*Hr3(0,1,0)*x + 72.D0*Hr3(0,1,0)*
c$$$     &    x**2 + 64.D0/3.D0*Hr3(0,1,0)*dx + 128.D0*Hr3(0,1,1) - 400.D0*
c$$$     &    Hr3(0,1,1)*x )
c$$$      PTqg2 = PTqg2 + nf*ca**2 * ( 640.D0/3.D0*Hr3(0,1,1)*x**2 + 32.D0/
c$$$     &    3.D0*Hr3(0,1,1)*dx + 728.D0/3.D0*Hr3(1,0,0) + 188.D0/3.D0*
c$$$     &    Hr3(1,0,0)*x - 304.D0/3.D0*Hr3(1,0,0)*x**2 + 224.D0/3.D0*Hr3(
c$$$     &    1,0,0)*dx + 160.D0/3.D0*Hr3(1,0,1) + 304.D0/3.D0*Hr3(1,0,1)*x
c$$$     &     - 472.D0/3.D0*Hr3(1,0,1)*x**2 + 32.D0*Hr3(1,0,1)*dx + 160.D0/
c$$$     &    3.D0*Hr3(1,1,0) + 304.D0/3.D0*Hr3(1,1,0)*x - 472.D0/3.D0*Hr3(
c$$$     &    1,1,0)*x**2 + 32.D0*Hr3(1,1,0)*dx - 196.D0/3.D0*Hr3(1,1,1) + 
c$$$     &    344.D0/3.D0*Hr3(1,1,1)*x - 400.D0/3.D0*Hr3(1,1,1)*x**2 + 32.D0
c$$$     &    /3.D0*Hr3(1,1,1)*dx - 96.D0*Hr4(-1,-1,-1,0) - 192.D0*Hr4(-1,
c$$$     &    -1,-1,0)*x - 192.D0*Hr4(-1,-1,-1,0)*x**2 + 64.D0*Hr4(-1,-1,0,
c$$$     &    0) + 128.D0*Hr4(-1,-1,0,0)*x + 128.D0*Hr4(-1,-1,0,0)*x**2 + 
c$$$     &    160.D0*Hr4(-1,-1,0,1) + 320.D0*Hr4(-1,-1,0,1)*x + 320.D0*Hr4(
c$$$     &    -1,-1,0,1)*x**2 + 112.D0*Hr4(-1,0,-1,0) + 224.D0*Hr4(-1,0,-1,
c$$$     &    0)*x + 224.D0*Hr4(-1,0,-1,0)*x**2 + 8.D0*Hr4(-1,0,0,0) + 16.D0
c$$$     &    *Hr4(-1,0,0,0)*x + 16.D0*Hr4(-1,0,0,0)*x**2 - 64.D0*Hr4(-1,0,
c$$$     &    1,0) )
c$$$      PTqg2 = PTqg2 + nf*ca**2 * (  - 128.D0*Hr4(-1,0,1,0)*x - 128.D0*
c$$$     &    Hr4(-1,0,1,0)*x**2 - 32.D0*Hr4(-1,0,1,1) - 64.D0*Hr4(-1,0,1,1
c$$$     &    )*x - 64.D0*Hr4(-1,0,1,1)*x**2 + 16.D0*Hr4(0,-1,-1,0) - 352.D0
c$$$     &    *Hr4(0,-1,-1,0)*x - 64.D0*Hr4(0,-1,-1,0)*x**2 - 104.D0*Hr4(0,
c$$$     &    -1,0,0) + 432.D0*Hr4(0,-1,0,0)*x - 32.D0*Hr4(0,-1,0,0)*x**2
c$$$     &     - 64.D0*Hr4(0,-1,0,1) - 128.D0*Hr4(0,-1,0,1)*x**2 - 56.D0*
c$$$     &    Hr4(0,0,-1,0) + 496.D0*Hr4(0,0,-1,0)*x + 64.D0*Hr4(0,0,0,0)
c$$$     &     - 1264.D0*Hr4(0,0,0,0)*x - 100.D0*Hr4(0,0,0,1) - 776.D0*Hr4(
c$$$     &    0,0,0,1)*x + 96.D0*Hr4(0,0,0,1)*x**2 + 64.D0*Hr4(0,0,1,0) - 
c$$$     &    96.D0*Hr4(0,0,1,0)*x + 160.D0*Hr4(0,0,1,0)*x**2 + 176.D0*Hr4(
c$$$     &    0,0,1,1) - 192.D0*Hr4(0,0,1,1)*x + 320.D0*Hr4(0,0,1,1)*x**2
c$$$     &     + 92.D0*Hr4(0,1,0,0) + 408.D0*Hr4(0,1,0,0)*x - 32.D0*Hr4(0,1
c$$$     &    ,0,0)*x**2 + 288.D0*Hr4(0,1,0,1)*x - 96.D0*Hr4(0,1,0,1)*x**2
c$$$     &     + 288.D0*Hr4(0,1,1,0)*x - 96.D0*Hr4(0,1,1,0)*x**2 - 80.D0*
c$$$     &    Hr4(0,1,1,1) + 256.D0*Hr4(0,1,1,1)*x - 192.D0*Hr4(0,1,1,1)*
c$$$     &    x**2 )
c$$$      PTqg2 = PTqg2 + nf*ca**2 * (  - 16.D0*Hr4(1,0,-1,0) + 32.D0*Hr4(1
c$$$     &    ,0,-1,0)*x - 32.D0*Hr4(1,0,-1,0)*x**2 - 120.D0*Hr4(1,0,0,0)
c$$$     &     + 240.D0*Hr4(1,0,0,0)*x - 240.D0*Hr4(1,0,0,0)*x**2 - 128.D0*
c$$$     &    Hr4(1,0,0,1) + 256.D0*Hr4(1,0,0,1)*x - 256.D0*Hr4(1,0,0,1)*
c$$$     &    x**2 - 32.D0*Hr4(1,0,1,1) + 64.D0*Hr4(1,0,1,1)*x - 64.D0*Hr4(
c$$$     &    1,0,1,1)*x**2 - 160.D0*Hr4(1,1,0,0) + 320.D0*Hr4(1,1,0,0)*x
c$$$     &     - 320.D0*Hr4(1,1,0,0)*x**2 - 32.D0*Hr4(1,1,0,1) + 64.D0*Hr4(
c$$$     &    1,1,0,1)*x - 64.D0*Hr4(1,1,0,1)*x**2 - 32.D0*Hr4(1,1,1,0) + 
c$$$     &    64.D0*Hr4(1,1,1,0)*x - 64.D0*Hr4(1,1,1,0)*x**2 + 32.D0*Hr4(1,
c$$$     &    1,1,1) - 64.D0*Hr4(1,1,1,1)*x + 64.D0*Hr4(1,1,1,1)*x**2 )
c$$$      PTqg2 = PTqg2 + nf*cf*ca * ( 17597.D0/36.D0 + 2659.D0/6.D0*x - 
c$$$     &    3092.D0/3.D0*x**2 - 220.D0/3.D0*dx + 1340.D0/3.D0*z3 + 4904.D0
c$$$     &    /3.D0*z3*x + 488.D0*z3*x**2 + 64.D0*z3*dx + 28.D0/9.D0*z2 + 
c$$$     &    4540.D0/9.D0*z2*x - 2620.D0/3.D0*z2*x**2 - 44.D0/5.D0*z2**2
c$$$     &     + 2632.D0/5.D0*z2**2*x - 128.D0*z2**2*x**2 - 40.D0*Hr1(-1)*
c$$$     &    z3 - 80.D0*Hr1(-1)*z3*x - 80.D0*Hr1(-1)*z3*x**2 - 96.D0*Hr1(
c$$$     &    -1)*z2 - 160.D0*Hr1(-1)*z2*x + 32.D0*Hr1(-1)*z2*x**2 - 2422.D0
c$$$     &    /27.D0*Hr1(0) + 30203.D0/27.D0*Hr1(0)*x - 17918.D0/27.D0*Hr1(
c$$$     &    0)*x**2 - 56.D0*Hr1(0)*z3 + 16.D0*Hr1(0)*z3*x - 160.D0*Hr1(0)
c$$$     &    *z3*x**2 - 260.D0/3.D0*Hr1(0)*z2 - 176.D0/3.D0*Hr1(0)*z2*x - 
c$$$     &    8.D0*Hr1(0)*z2*x**2 + 2788.D0/27.D0*Hr1(1) + 19660.D0/27.D0*
c$$$     &    Hr1(1)*x - 19970.D0/27.D0*Hr1(1)*x**2 - 80.D0/27.D0*Hr1(1)*dx
c$$$     &     - 792.D0*Hr1(1)*z3 + 1584.D0*Hr1(1)*z3*x - 1584.D0*Hr1(1)*z3
c$$$     &    *x**2 - 188.D0*Hr1(1)*z2 + 272.D0*Hr1(1)*z2*x - 120.D0*Hr1(1)
c$$$     &    *z2*x**2 + 32.D0*Hr1(1)*z2*dx + 80.D0*Hr2(-1,-1)*z2 + 160.D0*
c$$$     &    Hr2(-1,-1)*z2*x )
c$$$      PTqg2 = PTqg2 + nf*cf*ca * ( 160.D0*Hr2(-1,-1)*z2*x**2 + 672.D0*
c$$$     &    Hr2(-1,0) + 544.D0*Hr2(-1,0)*x - 16.D0*Hr2(-1,0)*x**2 - 64.D0
c$$$     &    *Hr2(0,-1)*z2 - 32.D0*Hr2(0,-1)*z2*x - 128.D0*Hr2(0,-1)*z2*
c$$$     &    x**2 - 1007.D0/9.D0*Hr2(0,0) - 6254.D0/9.D0*Hr2(0,0)*x + 2444.
c$$$     &    D0/9.D0*Hr2(0,0)*x**2 + 32.D0*Hr2(0,0)*z2 - 256.D0*Hr2(0,0)*
c$$$     &    z2*x + 256.D0*Hr2(0,0)*z2*x**2 + 2204.D0/9.D0*Hr2(0,1) - 3772.
c$$$     &    D0/9.D0*Hr2(0,1)*x + 3176.D0/3.D0*Hr2(0,1)*x**2 - 200.D0*Hr2(
c$$$     &    0,1)*z2 + 592.D0*Hr2(0,1)*z2*x - 512.D0*Hr2(0,1)*z2*x**2 - 
c$$$     &    212.D0*Hr2(1,0) + 2036.D0/3.D0*Hr2(1,0)*x - 3284.D0/3.D0*Hr2(
c$$$     &    1,0)*x**2 + 304.D0/3.D0*Hr2(1,0)*dx + 80.D0*Hr2(1,0)*z2 - 160.
c$$$     &    D0*Hr2(1,0)*z2*x + 160.D0*Hr2(1,0)*z2*x**2 - 488.D0/9.D0*Hr2(
c$$$     &    1,1) + 868.D0/9.D0*Hr2(1,1)*x - 1904.D0/9.D0*Hr2(1,1)*x**2 + 
c$$$     &    304.D0/9.D0*Hr2(1,1)*dx + 32.D0*Hr2(1,1)*z2 - 64.D0*Hr2(1,1)*
c$$$     &    z2*x + 64.D0*Hr2(1,1)*z2*x**2 - 288.D0*Hr3(-1,-1,0) - 448.D0*
c$$$     &    Hr3(-1,-1,0)*x - 64.D0*Hr3(-1,-1,0)*x**2 + 16.D0*Hr3(-1,0,0)
c$$$     &     - 80.D0*Hr3(-1,0,0)*x )
c$$$      PTqg2 = PTqg2 + nf*cf*ca * (  - 144.D0*Hr3(-1,0,0)*x**2 - 48.D0*
c$$$     &    Hr3(-1,0,1) - 64.D0*Hr3(-1,0,1)*x - 64.D0*Hr3(-1,0,1)*x**2 + 
c$$$     &    336.D0*Hr3(0,-1,0) + 16.D0*Hr3(0,-1,0)*x + 160.D0*Hr3(0,-1,0)
c$$$     &    *x**2 - 148.D0/3.D0*Hr3(0,0,0) + 1832.D0/3.D0*Hr3(0,0,0)*x + 
c$$$     &    640.D0/3.D0*Hr3(0,0,0)*x**2 + 388.D0/3.D0*Hr3(0,0,1) + 256.D0/
c$$$     &    3.D0*Hr3(0,0,1)*x - 40.D0/3.D0*Hr3(0,0,1)*x**2 + 28.D0/3.D0*
c$$$     &    Hr3(0,1,0) + 2032.D0/3.D0*Hr3(0,1,0)*x + 280.D0/3.D0*Hr3(0,1,
c$$$     &    0)*x**2 - 592.D0/3.D0*Hr3(0,1,1) + 1352.D0/3.D0*Hr3(0,1,1)*x
c$$$     &     - 192.D0*Hr3(0,1,1)*x**2 - 16.D0*Hr3(1,0,0) - 504.D0*Hr3(1,0
c$$$     &    ,0)*x + 368.D0*Hr3(1,0,0)*x**2 - 96.D0*Hr3(1,0,0)*dx + 268.D0/
c$$$     &    3.D0*Hr3(1,0,1) - 608.D0/3.D0*Hr3(1,0,1)*x + 728.D0/3.D0*Hr3(
c$$$     &    1,0,1)*x**2 - 32.D0*Hr3(1,0,1)*dx + 452.D0/3.D0*Hr3(1,1,0) - 
c$$$     &    1024.D0/3.D0*Hr3(1,1,0)*x + 1048.D0/3.D0*Hr3(1,1,0)*x**2 - 32.
c$$$     &    D0*Hr3(1,1,0)*dx + 424.D0/3.D0*Hr3(1,1,1) - 536.D0/3.D0*Hr3(1
c$$$     &    ,1,1)*x + 544.D0/3.D0*Hr3(1,1,1)*x**2 - 32.D0/3.D0*Hr3(1,1,1)
c$$$     &    *dx )
c$$$      PTqg2 = PTqg2 + nf*cf*ca * ( 96.D0*Hr4(-1,-1,-1,0) + 192.D0*Hr4(
c$$$     &    -1,-1,-1,0)*x + 192.D0*Hr4(-1,-1,-1,0)*x**2 + 64.D0*Hr4(-1,-1
c$$$     &    ,0,0) + 128.D0*Hr4(-1,-1,0,0)*x + 128.D0*Hr4(-1,-1,0,0)*x**2
c$$$     &     - 32.D0*Hr4(-1,-1,0,1) - 64.D0*Hr4(-1,-1,0,1)*x - 64.D0*Hr4(
c$$$     &    -1,-1,0,1)*x**2 - 80.D0*Hr4(-1,0,-1,0) - 160.D0*Hr4(-1,0,-1,0
c$$$     &    )*x - 160.D0*Hr4(-1,0,-1,0)*x**2 - 24.D0*Hr4(-1,0,0,0) - 48.D0
c$$$     &    *Hr4(-1,0,0,0)*x - 48.D0*Hr4(-1,0,0,0)*x**2 + 96.D0*Hr4(-1,0,
c$$$     &    1,0) + 192.D0*Hr4(-1,0,1,0)*x + 192.D0*Hr4(-1,0,1,0)*x**2 + 
c$$$     &    32.D0*Hr4(-1,0,1,1) + 64.D0*Hr4(-1,0,1,1)*x + 64.D0*Hr4(-1,0,
c$$$     &    1,1)*x**2 - 128.D0*Hr4(0,-1,-1,0) - 64.D0*Hr4(0,-1,-1,0)*x - 
c$$$     &    256.D0*Hr4(0,-1,-1,0)*x**2 - 96.D0*Hr4(0,-1,0,0)*x - 128.D0*
c$$$     &    Hr4(0,-1,0,0)*x**2 + 144.D0*Hr4(0,0,-1,0) + 32.D0*Hr4(0,0,-1,
c$$$     &    0)*x + 384.D0*Hr4(0,0,-1,0)*x**2 - 32.D0*Hr4(0,0,0,0) - 64.D0
c$$$     &    *Hr4(0,0,0,0)*x + 64.D0*Hr4(0,0,0,1) + 96.D0*Hr4(0,0,0,1)*x
c$$$     &     + 128.D0*Hr4(0,0,0,1)*x**2 - 160.D0*Hr4(0,0,1,0) - 160.D0*
c$$$     &    Hr4(0,0,1,0)*x )
c$$$      PTqg2 = PTqg2 + nf*cf*ca * (  - 320.D0*Hr4(0,0,1,0)*x**2 - 176.D0
c$$$     &    *Hr4(0,0,1,1) + 192.D0*Hr4(0,0,1,1)*x - 384.D0*Hr4(0,0,1,1)*
c$$$     &    x**2 + 136.D0*Hr4(0,1,0,0) - 1040.D0*Hr4(0,1,0,0)*x + 672.D0*
c$$$     &    Hr4(0,1,0,0)*x**2 + 168.D0*Hr4(0,1,0,1) - 624.D0*Hr4(0,1,0,1)
c$$$     &    *x + 448.D0*Hr4(0,1,0,1)*x**2 + 184.D0*Hr4(0,1,1,0) - 656.D0*
c$$$     &    Hr4(0,1,1,0)*x + 512.D0*Hr4(0,1,1,0)*x**2 + 136.D0*Hr4(0,1,1,
c$$$     &    1) - 368.D0*Hr4(0,1,1,1)*x + 320.D0*Hr4(0,1,1,1)*x**2 + 80.D0
c$$$     &    *Hr4(1,0,-1,0) - 160.D0*Hr4(1,0,-1,0)*x + 160.D0*Hr4(1,0,-1,0
c$$$     &    )*x**2 - 56.D0*Hr4(1,0,0,0) + 112.D0*Hr4(1,0,0,0)*x - 112.D0*
c$$$     &    Hr4(1,0,0,0)*x**2 + 112.D0*Hr4(1,0,0,1) - 224.D0*Hr4(1,0,0,1)
c$$$     &    *x + 224.D0*Hr4(1,0,0,1)*x**2 - 240.D0*Hr4(1,0,1,0) + 480.D0*
c$$$     &    Hr4(1,0,1,0)*x - 480.D0*Hr4(1,0,1,0)*x**2 - 32.D0*Hr4(1,0,1,1
c$$$     &    ) + 64.D0*Hr4(1,0,1,1)*x - 64.D0*Hr4(1,0,1,1)*x**2 + 208.D0*
c$$$     &    Hr4(1,1,0,0) - 416.D0*Hr4(1,1,0,0)*x + 416.D0*Hr4(1,1,0,0)*
c$$$     &    x**2 - 48.D0*Hr4(1,1,0,1) + 96.D0*Hr4(1,1,0,1)*x - 96.D0*Hr4(
c$$$     &    1,1,0,1)*x**2 )
c$$$      PTqg2 = PTqg2 + nf*cf*ca * (  - 16.D0*Hr4(1,1,1,0) + 32.D0*Hr4(1,
c$$$     &    1,1,0)*x - 32.D0*Hr4(1,1,1,0)*x**2 - 64.D0*Hr4(1,1,1,1) + 128.
c$$$     &    D0*Hr4(1,1,1,1)*x - 128.D0*Hr4(1,1,1,1)*x**2 )
c$$$      PTqg2 = PTqg2 + nf*cf**2 * (  - 227.D0/4.D0 - 1967.D0/2.D0*x + 
c$$$     &    1069.D0*x**2 + 156.D0*z3 - 816.D0*z3*x - 640.D0*z3*x**2 + 106.
c$$$     &    D0*z2 - 564.D0*z2*x + 564.D0*z2*x**2 + 56.D0*z2**2 - 1904.D0/
c$$$     &    5.D0*z2**2*x + 352.D0/5.D0*z2**2*x**2 + 128.D0*Hr1(-1)*z2 + 
c$$$     &    160.D0*Hr1(-1)*z2*x + 32.D0*Hr1(-1)*z2*x**2 + 92.D0*Hr1(0) - 
c$$$     &    1143.D0*Hr1(0)*x + 110.D0*Hr1(0)*x**2 + 40.D0*Hr1(0)*z3 - 304.
c$$$     &    D0*Hr1(0)*z3*x + 160.D0*Hr1(0)*z3*x**2 + 12.D0*Hr1(0)*z2 + 
c$$$     &    272.D0*Hr1(0)*z2*x - 176.D0*Hr1(0)*z2*x**2 + 180.D0*Hr1(1) - 
c$$$     &    480.D0*Hr1(1)*x + 318.D0*Hr1(1)*x**2 + 496.D0*Hr1(1)*z3 - 992.
c$$$     &    D0*Hr1(1)*z3*x + 992.D0*Hr1(1)*z3*x**2 + 108.D0*Hr1(1)*z2 - 
c$$$     &    64.D0*Hr1(1)*z2*x - 80.D0*Hr1(1)*z2*x**2 - 496.D0*Hr2(-1,0)
c$$$     &     - 240.D0*Hr2(-1,0)*x + 224.D0*Hr2(-1,0)*x**2 - 96.D0*Hr2(-1,
c$$$     &    0)*z2 - 192.D0*Hr2(-1,0)*z2*x - 192.D0*Hr2(-1,0)*z2*x**2 + 64.
c$$$     &    D0*Hr2(0,-1)*z2 + 64.D0*Hr2(0,-1)*z2*x + 128.D0*Hr2(0,-1)*z2*
c$$$     &    x**2 + 107.D0*Hr2(0,0) - 42.D0*Hr2(0,0)*x - 76.D0*Hr2(0,0)*
c$$$     &    x**2 )
c$$$      PTqg2 = PTqg2 + nf*cf**2 * (  - 40.D0*Hr2(0,0)*z2 + 144.D0*Hr2(0,
c$$$     &    0)*z2*x - 160.D0*Hr2(0,0)*z2*x**2 - 2.D0*Hr2(0,1) + 116.D0*
c$$$     &    Hr2(0,1)*x - 324.D0*Hr2(0,1)*x**2 + 104.D0*Hr2(0,1)*z2 - 144.D
c$$$     &    0*Hr2(0,1)*z2*x + 224.D0*Hr2(0,1)*z2*x**2 + 294.D0*Hr2(1,0)
c$$$     &     - 560.D0*Hr2(1,0)*x + 468.D0*Hr2(1,0)*x**2 - 112.D0*Hr2(1,0)
c$$$     &    *z2 + 224.D0*Hr2(1,0)*z2*x - 224.D0*Hr2(1,0)*z2*x**2 + 34.D0*
c$$$     &    Hr2(1,1) - 32.D0*Hr2(1,1)*x - 4.D0*Hr2(1,1)*x**2 - 16.D0*Hr2(
c$$$     &    1,1)*z2 + 32.D0*Hr2(1,1)*z2*x - 32.D0*Hr2(1,1)*z2*x**2 + 256.D
c$$$     &    0*Hr3(-1,-1,0) + 320.D0*Hr3(-1,-1,0)*x + 64.D0*Hr3(-1,-1,0)*
c$$$     &    x**2 + 64.D0*Hr3(-1,0,0) + 320.D0*Hr3(-1,0,0)*x + 256.D0*Hr3(
c$$$     &    -1,0,0)*x**2 - 240.D0*Hr3(0,-1,0) + 64.D0*Hr3(0,-1,0)*x - 256.
c$$$     &    D0*Hr3(0,-1,0)*x**2 + 36.D0*Hr3(0,0,0) - 584.D0*Hr3(0,0,0)*x
c$$$     &     - 192.D0*Hr3(0,0,0)*x**2 - 28.D0*Hr3(0,0,1) - 144.D0*Hr3(0,0
c$$$     &    ,1)*x + 16.D0*Hr3(0,0,1)*x**2 + 184.D0*Hr3(0,1,0) - 128.D0*
c$$$     &    Hr3(0,1,0)*x - 112.D0*Hr3(0,1,0)*x**2 + 104.D0*Hr3(0,1,1) - 
c$$$     &    80.D0*Hr3(0,1,1)*x )
c$$$      PTqg2 = PTqg2 + nf*cf**2 * ( 16.D0*Hr3(0,1,1)*x**2 - 288.D0*Hr3(1
c$$$     &    ,0,0) + 352.D0*Hr3(1,0,0)*x - 160.D0*Hr3(1,0,0)*x**2 - 124.D0
c$$$     &    *Hr3(1,0,1) + 64.D0*Hr3(1,0,1)*x - 48.D0*Hr3(1,0,1)*x**2 - 
c$$$     &    180.D0*Hr3(1,1,0) + 192.D0*Hr3(1,1,0)*x - 144.D0*Hr3(1,1,0)*
c$$$     &    x**2 - 76.D0*Hr3(1,1,1) + 64.D0*Hr3(1,1,1)*x - 48.D0*Hr3(1,1,
c$$$     &    1)*x**2 - 192.D0*Hr4(-1,-1,0,0) - 384.D0*Hr4(-1,-1,0,0)*x - 
c$$$     &    384.D0*Hr4(-1,-1,0,0)*x**2 - 96.D0*Hr4(-1,0,-1,0) - 192.D0*
c$$$     &    Hr4(-1,0,-1,0)*x - 192.D0*Hr4(-1,0,-1,0)*x**2 + 176.D0*Hr4(-1
c$$$     &    ,0,0,0) + 352.D0*Hr4(-1,0,0,0)*x + 352.D0*Hr4(-1,0,0,0)*x**2
c$$$     &     + 128.D0*Hr4(0,-1,-1,0) + 128.D0*Hr4(0,-1,-1,0)*x + 256.D0*
c$$$     &    Hr4(0,-1,-1,0)*x**2 + 32.D0*Hr4(0,-1,0,0) + 128.D0*Hr4(0,-1,0
c$$$     &    ,0)*x + 256.D0*Hr4(0,-1,0,0)*x**2 - 96.D0*Hr4(0,0,-1,0) + 64.D
c$$$     &    0*Hr4(0,0,-1,0)*x - 256.D0*Hr4(0,0,-1,0)*x**2 + 16.D0*Hr4(0,0
c$$$     &    ,0,0) - 256.D0*Hr4(0,0,0,0)*x + 128.D0*Hr4(0,0,0,0)*x**2 + 8.D
c$$$     &    0*Hr4(0,0,0,1) - 16.D0*Hr4(0,0,0,1)*x + 32.D0*Hr4(0,0,0,1)*
c$$$     &    x**2 )
c$$$      PTqg2 = PTqg2 + nf*cf**2 * ( 96.D0*Hr4(0,0,1,0) - 192.D0*Hr4(0,0,
c$$$     &    1,0)*x + 352.D0*Hr4(0,0,1,0)*x**2 + 64.D0*Hr4(0,0,1,1) - 128.D
c$$$     &    0*Hr4(0,0,1,1)*x + 192.D0*Hr4(0,0,1,1)*x**2 - 176.D0*Hr4(0,1,
c$$$     &    0,0) + 288.D0*Hr4(0,1,0,0)*x - 416.D0*Hr4(0,1,0,0)*x**2 - 104.
c$$$     &    D0*Hr4(0,1,0,1) + 208.D0*Hr4(0,1,0,1)*x - 224.D0*Hr4(0,1,0,1)
c$$$     &    *x**2 - 120.D0*Hr4(0,1,1,0) + 240.D0*Hr4(0,1,1,0)*x - 288.D0*
c$$$     &    Hr4(0,1,1,0)*x**2 - 56.D0*Hr4(0,1,1,1) + 112.D0*Hr4(0,1,1,1)*
c$$$     &    x - 128.D0*Hr4(0,1,1,1)*x**2 - 96.D0*Hr4(1,0,-1,0) + 192.D0*
c$$$     &    Hr4(1,0,-1,0)*x - 192.D0*Hr4(1,0,-1,0)*x**2 + 336.D0*Hr4(1,0,
c$$$     &    0,0) - 672.D0*Hr4(1,0,0,0)*x + 672.D0*Hr4(1,0,0,0)*x**2 + 80.D
c$$$     &    0*Hr4(1,0,0,1) - 160.D0*Hr4(1,0,0,1)*x + 160.D0*Hr4(1,0,0,1)*
c$$$     &    x**2 + 304.D0*Hr4(1,0,1,0) - 608.D0*Hr4(1,0,1,0)*x + 608.D0*
c$$$     &    Hr4(1,0,1,0)*x**2 + 64.D0*Hr4(1,0,1,1) - 128.D0*Hr4(1,0,1,1)*
c$$$     &    x + 128.D0*Hr4(1,0,1,1)*x**2 + 16.D0*Hr4(1,1,0,0) - 32.D0*
c$$$     &    Hr4(1,1,0,0)*x + 32.D0*Hr4(1,1,0,0)*x**2 + 80.D0*Hr4(1,1,0,1)
c$$$     &     - 160.D0*Hr4(1,1,0,1)*x )
c$$$      PTqg2 = PTqg2 + nf*cf**2 * ( 160.D0*Hr4(1,1,0,1)*x**2 + 48.D0*
c$$$     &    Hr4(1,1,1,0) - 96.D0*Hr4(1,1,1,0)*x + 96.D0*Hr4(1,1,1,0)*x**2
c$$$     &     + 32.D0*Hr4(1,1,1,1) - 64.D0*Hr4(1,1,1,1)*x + 64.D0*Hr4(1,1,
c$$$     &    1,1)*x**2 )
c$$$      PTqg2 = PTqg2 + nf2*ca * (  - 14.D0/9.D0 - 1216.D0/9.D0*x + 916.
c$$$     &    D0/9.D0*x**2 - 44.D0/9.D0*dx - 40.D0*z3 + 352.D0/3.D0*z3*x - 
c$$$     &    272.D0/3.D0*z3*x**2 + 100.D0/3.D0*z2 + 208.D0/3.D0*z2*x - 16.D
c$$$     &    0/9.D0*z2*x**2 + 32.D0/9.D0*z2*dx + 64.D0/3.D0*Hr1(-1)*z2 + 
c$$$     &    128.D0/3.D0*Hr1(-1)*z2*x + 128.D0/3.D0*Hr1(-1)*z2*x**2 + 44.D0
c$$$     &    /27.D0*Hr1(0) + 752.D0/27.D0*Hr1(0)*x - 2212.D0/27.D0*Hr1(0)*
c$$$     &    x**2 - 256.D0/27.D0*Hr1(0)*dx - 24.D0*Hr1(0)*z2 + 112.D0/3.D0
c$$$     &    *Hr1(0)*z2*x - 112.D0/3.D0*Hr1(0)*z2*x**2 + 592.D0/27.D0*Hr1(
c$$$     &    1) - 1232.D0/27.D0*Hr1(1)*x + 2632.D0/27.D0*Hr1(1)*x**2 - 320.
c$$$     &    D0/27.D0*Hr1(1)*dx - 32.D0/3.D0*Hr1(1)*z2 + 64.D0/3.D0*Hr1(1)
c$$$     &    *z2*x - 64.D0/3.D0*Hr1(1)*z2*x**2 - 16.D0*Hr2(-1,0) + 32.D0/3.
c$$$     &    D0*Hr2(-1,0)*x + 512.D0/9.D0*Hr2(-1,0)*x**2 + 32.D0/9.D0*Hr2(
c$$$     &    -1,0)*dx - 296.D0/9.D0*Hr2(0,0) - 416.D0/9.D0*Hr2(0,0)*x - 
c$$$     &    536.D0/9.D0*Hr2(0,0)*x**2 - 64.D0/9.D0*Hr2(0,0)*dx + 356.D0/9.
c$$$     &    D0*Hr2(0,1) - 1360.D0/9.D0*Hr2(0,1)*x + 632.D0/9.D0*Hr2(0,1)*
c$$$     &    x**2 )
c$$$      PTqg2 = PTqg2 + nf2*ca * ( 16.D0/9.D0*Hr2(1,0) + 496.D0/9.D0*
c$$$     &    Hr2(1,0)*x - 664.D0/9.D0*Hr2(1,0)*x**2 + 32.D0/3.D0*Hr2(1,0)*
c$$$     &    dx - 320.D0/9.D0*Hr2(1,1) + 544.D0/9.D0*Hr2(1,1)*x - 200.D0/3.
c$$$     &    D0*Hr2(1,1)*x**2 + 32.D0/9.D0*Hr2(1,1)*dx + 64.D0/3.D0*Hr3(-1
c$$$     &    ,-1,0) + 128.D0/3.D0*Hr3(-1,-1,0)*x + 128.D0/3.D0*Hr3(-1,-1,0
c$$$     &    )*x**2 + 32.D0/3.D0*Hr3(-1,0,0) + 64.D0/3.D0*Hr3(-1,0,0)*x + 
c$$$     &    64.D0/3.D0*Hr3(-1,0,0)*x**2 - 32.D0/3.D0*Hr3(-1,0,1) - 64.D0/
c$$$     &    3.D0*Hr3(-1,0,1)*x - 64.D0/3.D0*Hr3(-1,0,1)*x**2 - 64.D0/3.D0
c$$$     &    *Hr3(0,-1,0) + 128.D0/3.D0*Hr3(0,-1,0)*x - 64.D0/3.D0*Hr3(0,
c$$$     &    -1,0)*x**2 - 64.D0/3.D0*Hr3(0,0,0) - 640.D0/3.D0*Hr3(0,0,0)*x
c$$$     &     + 104.D0/3.D0*Hr3(0,0,1) - 16.D0*Hr3(0,0,1)*x + 176.D0/3.D0*
c$$$     &    Hr3(0,0,1)*x**2 + 8.D0*Hr3(0,1,0) + 80.D0*Hr3(0,1,0)*x - 16.D0
c$$$     &    *Hr3(0,1,0)*x**2 - 32.D0*Hr3(0,1,1) + 96.D0*Hr3(0,1,1)*x - 
c$$$     &    224.D0/3.D0*Hr3(0,1,1)*x**2 - 152.D0/3.D0*Hr3(1,0,0) + 304.D0/
c$$$     &    3.D0*Hr3(1,0,0)*x - 304.D0/3.D0*Hr3(1,0,0)*x**2 - 16.D0/3.D0*
c$$$     &    Hr3(1,0,1) )
c$$$      PTqg2 = PTqg2 + nf2*ca * ( 32.D0/3.D0*Hr3(1,0,1)*x - 32.D0/3.D0
c$$$     &    *Hr3(1,0,1)*x**2 - 16.D0/3.D0*Hr3(1,1,0) + 32.D0/3.D0*Hr3(1,1
c$$$     &    ,0)*x - 32.D0/3.D0*Hr3(1,1,0)*x**2 + 40.D0/3.D0*Hr3(1,1,1) - 
c$$$     &    80.D0/3.D0*Hr3(1,1,1)*x + 80.D0/3.D0*Hr3(1,1,1)*x**2 )
c$$$      PTqg2 = PTqg2 + nf2*cf * (  - 4847.D0/54.D0 - 2375.D0/27.D0*x
c$$$     &     + 4066.D0/27.D0*x**2 + 200.D0/27.D0*dx + 304.D0/3.D0*z3 - 
c$$$     &    608.D0/3.D0*z3*x + 160.D0*z3*x**2 + 272.D0/9.D0*z2 + 8.D0/9.D0
c$$$     &    *z2*x + 32.D0*z2*x**2 - 64.D0/9.D0*z2*dx + 1684.D0/27.D0*Hr1(
c$$$     &    0) - 4022.D0/27.D0*Hr1(0)*x - 4480.D0/27.D0*Hr1(0)*x**2 + 512.
c$$$     &    D0/27.D0*Hr1(0)*dx + 80.D0/3.D0*Hr1(0)*z2 - 352.D0/3.D0*Hr1(0
c$$$     &    )*z2*x + 48.D0*Hr1(0)*z2*x**2 - 808.D0/27.D0*Hr1(1) + 560.D0/
c$$$     &    27.D0*Hr1(1)*x - 248.D0/27.D0*Hr1(1)*x**2 - 8.D0*Hr1(1)*z2 + 
c$$$     &    16.D0*Hr1(1)*z2*x - 16.D0*Hr1(1)*z2*x**2 + 48.D0*Hr2(-1,0) - 
c$$$     &    496.D0/9.D0*Hr2(-1,0)*x**2 - 64.D0/9.D0*Hr2(-1,0)*dx + 902.D0/
c$$$     &    9.D0*Hr2(0,0) + 956.D0/9.D0*Hr2(0,0)*x + 1400.D0/9.D0*Hr2(0,0
c$$$     &    )*x**2 + 128.D0/9.D0*Hr2(0,0)*dx - 416.D0/9.D0*Hr2(0,1) + 376.
c$$$     &    D0/9.D0*Hr2(0,1)*x - 88.D0/3.D0*Hr2(0,1)*x**2 + 64.D0*Hr2(1,0
c$$$     &    ) - 248.D0/3.D0*Hr2(1,0)*x + 232.D0/3.D0*Hr2(1,0)*x**2 + 344.D
c$$$     &    0/9.D0*Hr2(1,1) - 376.D0/9.D0*Hr2(1,1)*x + 328.D0/9.D0*Hr2(1,
c$$$     &    1)*x**2 )
c$$$      PTqg2 = PTqg2 + nf2*cf * ( 32.D0*Hr3(0,-1,0) - 64.D0*Hr3(0,-1,0
c$$$     &    )*x + 64.D0/3.D0*Hr3(0,-1,0)*x**2 - 32.D0/3.D0*Hr3(0,0,0) + 
c$$$     &    112.D0/3.D0*Hr3(0,0,0)*x + 32.D0*Hr3(0,0,0)*x**2 - 64.D0/3.D0
c$$$     &    *Hr3(0,0,1) + 128.D0/3.D0*Hr3(0,0,1)*x - 176.D0/3.D0*Hr3(0,0,
c$$$     &    1)*x**2 + 80.D0/3.D0*Hr3(0,1,0) - 160.D0/3.D0*Hr3(0,1,0)*x + 
c$$$     &    176.D0/3.D0*Hr3(0,1,0)*x**2 + 64.D0/3.D0*Hr3(0,1,1) - 128.D0/
c$$$     &    3.D0*Hr3(0,1,1)*x + 48.D0*Hr3(0,1,1)*x**2 + 32.D0*Hr3(1,0,0)
c$$$     &     - 64.D0*Hr3(1,0,0)*x + 64.D0*Hr3(1,0,0)*x**2 - 40.D0/3.D0*
c$$$     &    Hr3(1,0,1) + 80.D0/3.D0*Hr3(1,0,1)*x - 80.D0/3.D0*Hr3(1,0,1)*
c$$$     &    x**2 - 56.D0/3.D0*Hr3(1,1,0) + 112.D0/3.D0*Hr3(1,1,0)*x - 112.
c$$$     &    D0/3.D0*Hr3(1,1,0)*x**2 - 40.D0/3.D0*Hr3(1,1,1) + 80.D0/3.D0*
c$$$     &    Hr3(1,1,1)*x - 80.D0/3.D0*Hr3(1,1,1)*x**2 - 32.D0*Hr4(0,0,0,0
c$$$     &    ) + 64.D0*Hr4(0,0,0,0)*x )
c$$$      PTqg2 = PTqg2 + nf**3 * (  - 8.D0/9.D0 + 16.D0/3.D0*x - 16.D0/3.D0
c$$$     &    *x**2 + 8.D0/3.D0*z2 - 16.D0/3.D0*z2*x + 16.D0/3.D0*z2*x**2
c$$$     &     + 40.D0/9.D0*Hr1(0) - 32.D0/9.D0*Hr1(0)*x + 32.D0/9.D0*Hr1(0
c$$$     &    )*x**2 - 40.D0/9.D0*Hr1(1) + 32.D0/9.D0*Hr1(1)*x - 32.D0/9.D0
c$$$     &    *Hr1(1)*x**2 + 8.D0/3.D0*Hr2(0,0) - 16.D0/3.D0*Hr2(0,0)*x + 
c$$$     &    16.D0/3.D0*Hr2(0,0)*x**2 - 8.D0/3.D0*Hr2(0,1) + 16.D0/3.D0*
c$$$     &    Hr2(0,1)*x - 16.D0/3.D0*Hr2(0,1)*x**2 - 8.D0/3.D0*Hr2(1,0) + 
c$$$     &    16.D0/3.D0*Hr2(1,0)*x - 16.D0/3.D0*Hr2(1,0)*x**2 + 8.D0/3.D0*
c$$$     &    Hr2(1,1) - 16.D0/3.D0*Hr2(1,1)*x + 16.D0/3.D0*Hr2(1,1)*x**2 )
c$$$      PTqg2 = PTqg2
c$$$     , + 2.*z2*(ca-cf)*b0* 2.*nf* (1.- 2.*x + 2.*x*x )
c$$$     ,     * ( 11.+24.*Hr1(0) ) * IMOD
c$$$*
c$$$       X2QGTA = PTqg2
c$$$*
c$$$       RETURN
c$$$       END
c$$$*
c$$$* ---------------------------------------------------------------------
c$$$*
c$$$*
c$$$* ..The gluon-quark splitting functions P_gq^(2)T
c$$$*
c$$$       FUNCTION X2GQTA (X, NF)
c$$$*
c$$$       IMPLICIT REAL*8 (A - Z)
c$$$       COMPLEX*16 HC1, HC2, HC3, HC4 
c$$$       INTEGER NF, NF2, N1, N2, NW, I1, I2, I3, N
c$$$       PARAMETER ( N1 = -1, N2 = 1, NW = 4 ) 
c$$$       DIMENSION HC1(N1:N2),HC2(N1:N2,N1:N2),HC3(N1:N2,N1:N2,N1:N2), 
c$$$     ,           HC4(N1:N2,N1:N2,N1:N2,N1:N2) 
c$$$       DIMENSION HR1(N1:N2),HR2(N1:N2,N1:N2),HR3(N1:N2,N1:N2,N1:N2), 
c$$$     ,           HR4(N1:N2,N1:N2,N1:N2,N1:N2) 
c$$$       DIMENSION HI1(N1:N2),HI2(N1:N2,N1:N2),HI3(N1:N2,N1:N2,N1:N2), 
c$$$     ,           HI4(N1:N2,N1:N2,N1:N2,N1:N2) 
c$$$       PARAMETER ( Z2 = 1.6449 34066 84822 64365 D0,
c$$$     ,             Z3 = 1.2020 56903 15959 42854 D0 )
c$$$*
c$$$* ...Colour factors and an abbreviation
c$$$*
c$$$       CF  = 4./3.D0
c$$$       CA  = 3.D0
c$$$       NF2 = NF*NF
c$$$*
c$$$       DX = 1.D0/X
c$$$*
c$$$* ...The harmonic polylogs up to weight 4 by Gehrmann and Remiddi
c$$$*
c$$$       CALL HPLOG (X, NW, HC1,HC2,HC3,HC4, HR1,HR2,HR3,HR4,
c$$$     ,            HI1,HI2,HI3,HI4, N1, N2) 
c$$$*
c$$$* ...The splitting function in terms of the harmonic polylogs
c$$$*
c$$$      PTgq2 =
c$$$     &  + cf*ca**2 * (  - 17798.D0/27.D0 - 3436.D0/27.D0*x + 4136.D0/9.D
c$$$     &    0*x**2 + 3514.D0/9.D0*dx - 400.D0*z3 + 148.D0*z3*x - 640.D0/3.
c$$$     &    D0*z3*x**2 + 160.D0*z3*dx - 848.D0*z2 - 2800.D0/9.D0*z2*x - 
c$$$     &    352.D0/9.D0*z2*x**2 - 5176.D0/9.D0*z2*dx + 572.D0/5.D0*z2**2
c$$$     &     + 314.D0*z2**2*x + 1984.D0/5.D0*z2**2*dx - 16.D0*Hr1(-1)*z3
c$$$     &     - 8.D0*Hr1(-1)*z3*x - 16.D0*Hr1(-1)*z3*dx - 168.D0*Hr1(-1)*
c$$$     &    z2 - 160.D0*Hr1(-1)*z2*x - 32.D0/3.D0*Hr1(-1)*z2*x**2 - 176.D0
c$$$     &    /3.D0*Hr1(-1)*z2*dx - 32612.D0/27.D0*Hr1(0) - 23630.D0/27.D0*
c$$$     &    Hr1(0)*x - 12964.D0/27.D0*Hr1(0)*x**2 + 13384.D0/27.D0*Hr1(0)
c$$$     &    *dx + 352.D0*Hr1(0)*z3 + 240.D0*Hr1(0)*z3*x + 416.D0*Hr1(0)*
c$$$     &    z3*dx - 12.D0*Hr1(0)*z2 + 532.D0/3.D0*Hr1(0)*z2*x - 128.D0/3.D
c$$$     &    0*Hr1(0)*z2*x**2 - 752.D0/3.D0*Hr1(0)*z2*dx - 22270.D0/27.D0*
c$$$     &    Hr1(1) - 688.D0/27.D0*Hr1(1)*x - 1924.D0/9.D0*Hr1(1)*x**2 + 
c$$$     &    31504.D0/27.D0*Hr1(1)*dx - 656.D0*Hr1(1)*z3 + 328.D0*Hr1(1)*
c$$$     &    z3*x + 656.D0*Hr1(1)*z3*dx - 280.D0*Hr1(1)*z2 - 32.D0*Hr1(1)*
c$$$     &    z2*x )
c$$$      PTgq2 = PTgq2 + cf*ca**2 * (  - 64.D0*Hr1(1)*z2*x**2 + 344.D0*
c$$$     &    Hr1(1)*z2*dx - 32.D0*Hr2(-1,-1)*z2 - 16.D0*Hr2(-1,-1)*z2*x - 
c$$$     &    32.D0*Hr2(-1,-1)*z2*dx - 1280.D0/9.D0*Hr2(-1,0) + 1052.D0/9.D0
c$$$     &    *Hr2(-1,0)*x + 704.D0/9.D0*Hr2(-1,0)*x**2 - 248.D0/3.D0*Hr2(
c$$$     &    -1,0)*dx + 160.D0*Hr2(-1,0)*z2 + 80.D0*Hr2(-1,0)*z2*x + 160.D0
c$$$     &    *Hr2(-1,0)*z2*dx - 144.D0*Hr2(0,-1)*z2 - 8.D0*Hr2(0,-1)*z2*x
c$$$     &     - 64.D0*Hr2(0,-1)*z2*dx + 11848.D0/9.D0*Hr2(0,0) + 7280.D0/9.
c$$$     &    D0*Hr2(0,0)*x + 3688.D0/9.D0*Hr2(0,0)*x**2 + 3920.D0/9.D0*
c$$$     &    Hr2(0,0)*dx + 152.D0*Hr2(0,0)*z2 - 156.D0*Hr2(0,0)*z2*x - 128.
c$$$     &    D0*Hr2(0,0)*z2*dx + 6896.D0/9.D0*Hr2(0,1) + 4220.D0/9.D0*Hr2(
c$$$     &    0,1)*x + 1552.D0/9.D0*Hr2(0,1)*x**2 + 4664.D0/9.D0*Hr2(0,1)*
c$$$     &    dx + 48.D0*Hr2(0,1)*z2 + 248.D0*Hr2(0,1)*z2*x + 320.D0*Hr2(0,
c$$$     &    1)*z2*dx + 1160.D0/9.D0*Hr2(1,0) - 604.D0/9.D0*Hr2(1,0)*x + 
c$$$     &    760.D0/3.D0*Hr2(1,0)*x**2 - 4136.D0/9.D0*Hr2(1,0)*dx + 128.D0
c$$$     &    *Hr2(1,0)*z2 - 64.D0*Hr2(1,0)*z2*x - 128.D0*Hr2(1,0)*z2*dx + 
c$$$     &    808.D0/9.D0*Hr2(1,1) )
c$$$      PTgq2 = PTgq2 + cf*ca**2 * (  - 440.D0/9.D0*Hr2(1,1)*x + 560.D0/9.
c$$$     &    D0*Hr2(1,1)*x**2 - 532.D0/3.D0*Hr2(1,1)*dx - 96.D0*Hr2(1,1)*
c$$$     &    z2 + 48.D0*Hr2(1,1)*z2*x + 96.D0*Hr2(1,1)*z2*dx - 16.D0*Hr3(
c$$$     &    -1,-1,0) - 160.D0*Hr3(-1,-1,0)*x + 64.D0/3.D0*Hr3(-1,-1,0)*
c$$$     &    x**2 + 352.D0/3.D0*Hr3(-1,-1,0)*dx - 328.D0/3.D0*Hr3(-1,0,0)
c$$$     &     + 544.D0/3.D0*Hr3(-1,0,0)*x - 128.D0/3.D0*Hr3(-1,0,0)*x**2
c$$$     &     - 272.D0*Hr3(-1,0,0)*dx + 160.D0*Hr3(-1,0,1) + 80.D0*Hr3(-1,
c$$$     &    0,1)*x + 64.D0/3.D0*Hr3(-1,0,1)*x**2 + 352.D0/3.D0*Hr3(-1,0,1
c$$$     &    )*dx - 280.D0/3.D0*Hr3(0,-1,0) + 928.D0/3.D0*Hr3(0,-1,0)*x - 
c$$$     &    256.D0/3.D0*Hr3(0,-1,0)*x**2 - 752.D0/3.D0*Hr3(0,-1,0)*dx - 
c$$$     &    848.D0/3.D0*Hr3(0,0,0) - 3344.D0/3.D0*Hr3(0,0,0)*x - 224.D0/3.
c$$$     &    D0*Hr3(0,0,0)*x**2 + 1856.D0/3.D0*Hr3(0,0,0)*dx - 1292.D0/3.D0
c$$$     &    *Hr3(0,0,1) - 284.D0/3.D0*Hr3(0,0,1)*x - 128.D0/3.D0*Hr3(0,0,
c$$$     &    1)*x**2 - 32.D0/3.D0*Hr3(0,0,1)*dx - 576.D0*Hr3(0,1,0) - 424.D
c$$$     &    0*Hr3(0,1,0)*x - 512.D0/3.D0*Hr3(0,1,0)*x**2 + 1504.D0/3.D0*
c$$$     &    Hr3(0,1,0)*dx )
c$$$      PTgq2 = PTgq2 + cf*ca**2 * (  - 296.D0/3.D0*Hr3(0,1,1) - 260.D0/3.
c$$$     &    D0*Hr3(0,1,1)*x - 128.D0/3.D0*Hr3(0,1,1)*x**2 + 448.D0/3.D0*
c$$$     &    Hr3(0,1,1)*dx + 1300.D0/3.D0*Hr3(1,0,0) - 1004.D0/3.D0*Hr3(1,
c$$$     &    0,0)*x - 1024.D0/3.D0*Hr3(1,0,0)*dx + 192.D0*Hr3(1,0,1) - 72.D
c$$$     &    0*Hr3(1,0,1)*x + 32.D0*Hr3(1,0,1)*x**2 - 248.D0*Hr3(1,0,1)*dx
c$$$     &     + 368.D0/3.D0*Hr3(1,1,0) - 352.D0/3.D0*Hr3(1,1,0)*x + 32.D0/
c$$$     &    3.D0*Hr3(1,1,0)*x**2 - 424.D0/3.D0*Hr3(1,1,0)*dx + 104.D0/3.D0
c$$$     &    *Hr3(1,1,1) - 76.D0/3.D0*Hr3(1,1,1)*x + 32.D0/3.D0*Hr3(1,1,1)
c$$$     &    *x**2 - 160.D0/3.D0*Hr3(1,1,1)*dx - 192.D0*Hr4(-1,-1,-1,0) - 
c$$$     &    96.D0*Hr4(-1,-1,-1,0)*x - 192.D0*Hr4(-1,-1,-1,0)*dx + 192.D0*
c$$$     &    Hr4(-1,-1,0,0) + 96.D0*Hr4(-1,-1,0,0)*x + 192.D0*Hr4(-1,-1,0,
c$$$     &    0)*dx - 64.D0*Hr4(-1,-1,0,1) - 32.D0*Hr4(-1,-1,0,1)*x - 64.D0
c$$$     &    *Hr4(-1,-1,0,1)*dx + 352.D0*Hr4(-1,0,-1,0) + 176.D0*Hr4(-1,0,
c$$$     &    -1,0)*x + 352.D0*Hr4(-1,0,-1,0)*dx - 336.D0*Hr4(-1,0,0,0) - 
c$$$     &    168.D0*Hr4(-1,0,0,0)*x - 336.D0*Hr4(-1,0,0,0)*dx + 64.D0*Hr4(
c$$$     &    -1,0,0,1) )
c$$$      PTgq2 = PTgq2 + cf*ca**2 * ( 32.D0*Hr4(-1,0,0,1)*x + 64.D0*Hr4(-1
c$$$     &    ,0,0,1)*dx - 192.D0*Hr4(-1,0,1,0) - 96.D0*Hr4(-1,0,1,0)*x - 
c$$$     &    192.D0*Hr4(-1,0,1,0)*dx - 64.D0*Hr4(-1,0,1,1) - 32.D0*Hr4(-1,
c$$$     &    0,1,1)*x - 64.D0*Hr4(-1,0,1,1)*dx - 160.D0*Hr4(0,-1,-1,0) + 
c$$$     &    48.D0*Hr4(0,-1,-1,0)*x + 240.D0*Hr4(0,-1,0,0) - 136.D0*Hr4(0,
c$$$     &    -1,0,0)*x - 96.D0*Hr4(0,-1,0,0)*dx + 64.D0*Hr4(0,-1,0,1) + 32.
c$$$     &    D0*Hr4(0,-1,0,1)*x + 64.D0*Hr4(0,-1,0,1)*dx + 112.D0*Hr4(0,0,
c$$$     &    -1,0) - 248.D0*Hr4(0,0,-1,0)*x - 256.D0*Hr4(0,0,-1,0)*dx - 
c$$$     &    128.D0*Hr4(0,0,0,0) + 632.D0*Hr4(0,0,0,0)*x + 512.D0*Hr4(0,0,
c$$$     &    0,0)*dx - 152.D0*Hr4(0,0,0,1) - 92.D0*Hr4(0,0,0,1)*x - 128.D0
c$$$     &    *Hr4(0,0,0,1)*dx + 384.D0*Hr4(0,0,1,0) + 480.D0*Hr4(0,0,1,0)*
c$$$     &    x + 576.D0*Hr4(0,0,1,0)*dx + 128.D0*Hr4(0,0,1,1) + 160.D0*
c$$$     &    Hr4(0,0,1,1)*x + 192.D0*Hr4(0,0,1,1)*dx + 232.D0*Hr4(0,1,0,0)
c$$$     &     - 76.D0*Hr4(0,1,0,0)*x - 160.D0*Hr4(0,1,0,0)*dx + 64.D0*Hr4(
c$$$     &    0,1,0,1) - 176.D0*Hr4(0,1,0,1)*x - 256.D0*Hr4(0,1,0,1)*dx + 
c$$$     &    128.D0*Hr4(0,1,1,0) )
c$$$      PTgq2 = PTgq2 + cf*ca**2 * (  - 112.D0*Hr4(0,1,1,0)*x - 192.D0*
c$$$     &    Hr4(0,1,1,0)*dx + 64.D0*Hr4(0,1,1,1) - 80.D0*Hr4(0,1,1,1)*x
c$$$     &     - 128.D0*Hr4(0,1,1,1)*dx + 96.D0*Hr4(1,0,-1,0) - 48.D0*Hr4(1
c$$$     &    ,0,-1,0)*x - 96.D0*Hr4(1,0,-1,0)*dx - 560.D0*Hr4(1,0,0,0) + 
c$$$     &    280.D0*Hr4(1,0,0,0)*x + 560.D0*Hr4(1,0,0,0)*dx - 160.D0*Hr4(1
c$$$     &    ,0,0,1) + 80.D0*Hr4(1,0,0,1)*x + 160.D0*Hr4(1,0,0,1)*dx - 640.
c$$$     &    D0*Hr4(1,0,1,0) + 320.D0*Hr4(1,0,1,0)*x + 640.D0*Hr4(1,0,1,0)
c$$$     &    *dx - 64.D0*Hr4(1,0,1,1) + 32.D0*Hr4(1,0,1,1)*x + 64.D0*Hr4(1
c$$$     &    ,0,1,1)*dx - 160.D0*Hr4(1,1,0,0) + 80.D0*Hr4(1,1,0,0)*x + 160.
c$$$     &    D0*Hr4(1,1,0,0)*dx - 128.D0*Hr4(1,1,0,1) + 64.D0*Hr4(1,1,0,1)
c$$$     &    *x + 128.D0*Hr4(1,1,0,1)*dx - 192.D0*Hr4(1,1,1,0) + 96.D0*
c$$$     &    Hr4(1,1,1,0)*x + 192.D0*Hr4(1,1,1,0)*dx - 64.D0*Hr4(1,1,1,1)
c$$$     &     + 32.D0*Hr4(1,1,1,1)*x + 64.D0*Hr4(1,1,1,1)*dx )
c$$$      PTgq2 = PTgq2 + cf**2*ca * ( 1735.D0/6.D0 - 4811.D0/12.D0*x - 140.
c$$$     &    D0/3.D0*x**2 + 200.D0*dx + 552.D0*z3 + 424.D0*z3*x + 128.D0/3.
c$$$     &    D0*z3*x**2 + 264.D0*z3*dx + 3452.D0/3.D0*z2 + 436.D0/3.D0*z2*
c$$$     &    x + 1760.D0/9.D0*z2*x**2 + 528.D0*z2*dx - 608.D0/5.D0*z2**2
c$$$     &     - 2144.D0/5.D0*z2**2*x - 3664.D0/5.D0*z2**2*dx - 176.D0*Hr1(
c$$$     &    -1)*z3 - 88.D0*Hr1(-1)*z3*x - 176.D0*Hr1(-1)*z3*dx - 64.D0*
c$$$     &    Hr1(-1)*z2 - 32.D0*Hr1(-1)*z2*x - 24.D0*Hr1(-1)*z2*dx + 40301.
c$$$     &    D0/27.D0*Hr1(0) + 18911.D0/27.D0*Hr1(0)*x + 8932.D0/27.D0*
c$$$     &    Hr1(0)*x**2 + 244.D0*Hr1(0)*dx - 272.D0*Hr1(0)*z3 - 376.D0*
c$$$     &    Hr1(0)*z3*x - 832.D0*Hr1(0)*z3*dx + 72.D0*Hr1(0)*z2 + 128.D0*
c$$$     &    Hr1(0)*z2*x - 32.D0/3.D0*Hr1(0)*z2*x**2 + 192.D0*Hr1(0)*z2*dx
c$$$     &     + 26008.D0/27.D0*Hr1(1) + 3988.D0/27.D0*Hr1(1)*x + 5908.D0/
c$$$     &    27.D0*Hr1(1)*x**2 - 39092.D0/27.D0*Hr1(1)*dx + 816.D0*Hr1(1)*
c$$$     &    z3 - 408.D0*Hr1(1)*z3*x - 816.D0*Hr1(1)*z3*dx + 672.D0*Hr1(1)
c$$$     &    *z2 - 180.D0*Hr1(1)*z2*x + 224.D0/3.D0*Hr1(1)*z2*x**2 - 1760.D
c$$$     &    0/3.D0*Hr1(1)*z2*dx )
c$$$      PTgq2 = PTgq2 + cf**2*ca * ( 288.D0*Hr2(-1,-1)*z2 + 144.D0*Hr2(-1
c$$$     &    ,-1)*z2*x + 288.D0*Hr2(-1,-1)*z2*dx - 88.D0*Hr2(-1,0) - 268.D0
c$$$     &    *Hr2(-1,0)*x + 192.D0*Hr2(-1,0)*dx - 192.D0*Hr2(-1,0)*z2 - 96.
c$$$     &    D0*Hr2(-1,0)*z2*x - 192.D0*Hr2(-1,0)*z2*dx - 128.D0*Hr2(0,-1)
c$$$     &    *z2 - 16.D0*Hr2(0,-1)*z2*x + 32.D0*Hr2(0,-1)*z2*dx - 2338.D0/
c$$$     &    3.D0*Hr2(0,0) - 441.D0*Hr2(0,0)*x - 1408.D0/9.D0*Hr2(0,0)*
c$$$     &    x**2 + 80.D0*Hr2(0,0)*dx - 16.D0*Hr2(0,0)*z2 - 40.D0*Hr2(0,0)
c$$$     &    *z2*x - 128.D0*Hr2(0,0)*z2*dx - 8980.D0/9.D0*Hr2(0,1) - 4576.D
c$$$     &    0/9.D0*Hr2(0,1)*x - 352.D0/3.D0*Hr2(0,1)*x**2 - 5800.D0/9.D0*
c$$$     &    Hr2(0,1)*dx + 272.D0*Hr2(0,1)*z2 - 456.D0*Hr2(0,1)*z2*x - 608.
c$$$     &    D0*Hr2(0,1)*z2*dx + 488.D0/3.D0*Hr2(1,0) - 200.D0*Hr2(1,0)*x
c$$$     &     - 704.D0/9.D0*Hr2(1,0)*x**2 + 620.D0/9.D0*Hr2(1,0)*dx - 32.D0
c$$$     &    *Hr2(1,0)*z2 + 16.D0*Hr2(1,0)*z2*x + 32.D0*Hr2(1,0)*z2*dx - 
c$$$     &    160.D0/9.D0*Hr2(1,1) - 280.D0/9.D0*Hr2(1,1)*x - 352.D0/9.D0*
c$$$     &    Hr2(1,1)*x**2 + 848.D0/9.D0*Hr2(1,1)*dx + 192.D0*Hr2(1,1)*z2
c$$$     &     - 96.D0*Hr2(1,1)*z2*x )
c$$$      PTgq2 = PTgq2 + cf**2*ca * (  - 192.D0*Hr2(1,1)*z2*dx + 192.D0*
c$$$     &    Hr3(-1,-1,0) + 128.D0*Hr3(-1,-1,0)*x + 48.D0*Hr3(-1,-1,0)*dx
c$$$     &     - 192.D0*Hr3(-1,0,0) - 72.D0*Hr3(-1,0,0)*x - 48.D0*Hr3(-1,0,
c$$$     &    0)*dx + 160.D0*Hr3(-1,0,1) + 96.D0*Hr3(-1,0,1)*x + 48.D0*Hr3(
c$$$     &    -1,0,1)*dx - 144.D0*Hr3(0,-1,0) + 128.D0*Hr3(0,-1,0)*x + 96.D0
c$$$     &    *Hr3(0,-1,0)*dx + 368.D0/3.D0*Hr3(0,0,0) - 268.D0/3.D0*Hr3(0,
c$$$     &    0,0)*x + 64.D0/3.D0*Hr3(0,0,0)*x**2 + 136.D0/3.D0*Hr3(0,0,1)
c$$$     &     - 680.D0/3.D0*Hr3(0,0,1)*x - 32.D0*Hr3(0,0,1)*x**2 + 640.D0/
c$$$     &    3.D0*Hr3(0,0,1)*dx + 600.D0*Hr3(0,1,0) + 96.D0*Hr3(0,1,0)*x
c$$$     &     + 224.D0/3.D0*Hr3(0,1,0)*x**2 - 1280.D0/3.D0*Hr3(0,1,0)*dx
c$$$     &     + 424.D0/3.D0*Hr3(0,1,1) + 28.D0/3.D0*Hr3(0,1,1)*x + 32.D0*
c$$$     &    Hr3(0,1,1)*x**2 - 464.D0/3.D0*Hr3(0,1,1)*dx - 1760.D0/3.D0*
c$$$     &    Hr3(1,0,0) + 832.D0/3.D0*Hr3(1,0,0)*x - 64.D0/3.D0*Hr3(1,0,0)
c$$$     &    *x**2 + 336.D0*Hr3(1,0,0)*dx - 976.D0/3.D0*Hr3(1,0,1) + 404.D0
c$$$     &    /3.D0*Hr3(1,0,1)*x - 32.D0*Hr3(1,0,1)*x**2 + 904.D0/3.D0*Hr3(
c$$$     &    1,0,1)*dx )
c$$$      PTgq2 = PTgq2 + cf**2*ca * (  - 256.D0/3.D0*Hr3(1,1,0) + 332.D0/3.
c$$$     &    D0*Hr3(1,1,0)*x - 32.D0/3.D0*Hr3(1,1,0)*x**2 + 88.D0*Hr3(1,1,
c$$$     &    0)*dx - 152.D0/3.D0*Hr3(1,1,1) + 232.D0/3.D0*Hr3(1,1,1)*x - 
c$$$     &    32.D0/3.D0*Hr3(1,1,1)*x**2 + 160.D0/3.D0*Hr3(1,1,1)*dx + 192.D
c$$$     &    0*Hr4(-1,-1,-1,0) + 96.D0*Hr4(-1,-1,-1,0)*x + 192.D0*Hr4(-1,
c$$$     &    -1,-1,0)*dx + 64.D0*Hr4(-1,-1,0,0) + 32.D0*Hr4(-1,-1,0,0)*x
c$$$     &     + 64.D0*Hr4(-1,-1,0,0)*dx - 192.D0*Hr4(-1,-1,0,1) - 96.D0*
c$$$     &    Hr4(-1,-1,0,1)*x - 192.D0*Hr4(-1,-1,0,1)*dx - 160.D0*Hr4(-1,0
c$$$     &    ,-1,0) - 80.D0*Hr4(-1,0,-1,0)*x - 160.D0*Hr4(-1,0,-1,0)*dx - 
c$$$     &    16.D0*Hr4(-1,0,0,0) - 8.D0*Hr4(-1,0,0,0)*x - 16.D0*Hr4(-1,0,0
c$$$     &    ,0)*dx - 64.D0*Hr4(-1,0,0,1) - 32.D0*Hr4(-1,0,0,1)*x - 64.D0*
c$$$     &    Hr4(-1,0,0,1)*dx + 128.D0*Hr4(-1,0,1,0) + 64.D0*Hr4(-1,0,1,0)
c$$$     &    *x + 128.D0*Hr4(-1,0,1,0)*dx + 64.D0*Hr4(-1,0,1,1) + 32.D0*
c$$$     &    Hr4(-1,0,1,1)*x + 64.D0*Hr4(-1,0,1,1)*dx - 128.D0*Hr4(0,-1,-1
c$$$     &    ,0) - 96.D0*Hr4(0,-1,-1,0)*x - 64.D0*Hr4(0,-1,-1,0)*dx - 32.D0
c$$$     &    *Hr4(0,-1,0,0) )
c$$$      PTgq2 = PTgq2 + cf**2*ca * (  - 32.D0*Hr4(0,-1,0,0)*x + 64.D0*
c$$$     &    Hr4(0,-1,0,1) - 32.D0*Hr4(0,-1,0,1)*x - 64.D0*Hr4(0,-1,0,1)*
c$$$     &    dx - 32.D0*Hr4(0,0,-1,0) + 48.D0*Hr4(0,0,-1,0)*x + 64.D0*Hr4(
c$$$     &    0,0,0,0) + 64.D0*Hr4(0,0,0,0)*x + 336.D0*Hr4(0,0,0,1) + 376.D0
c$$$     &    *Hr4(0,0,0,1)*x + 512.D0*Hr4(0,0,0,1)*dx - 208.D0*Hr4(0,0,1,0
c$$$     &    ) - 344.D0*Hr4(0,0,1,0)*x - 448.D0*Hr4(0,0,1,0)*dx - 16.D0*
c$$$     &    Hr4(0,0,1,1) - 216.D0*Hr4(0,0,1,1)*x - 320.D0*Hr4(0,0,1,1)*dx
c$$$     &     - 320.D0*Hr4(0,1,0,0) + 208.D0*Hr4(0,1,0,0)*x + 256.D0*Hr4(0
c$$$     &    ,1,0,0)*dx - 336.D0*Hr4(0,1,0,1) + 312.D0*Hr4(0,1,0,1)*x + 
c$$$     &    448.D0*Hr4(0,1,0,1)*dx - 336.D0*Hr4(0,1,1,0) + 216.D0*Hr4(0,1
c$$$     &    ,1,0)*x + 384.D0*Hr4(0,1,1,0)*dx - 272.D0*Hr4(0,1,1,1) + 184.D
c$$$     &    0*Hr4(0,1,1,1)*x + 320.D0*Hr4(0,1,1,1)*dx - 96.D0*Hr4(1,0,-1,
c$$$     &    0) + 48.D0*Hr4(1,0,-1,0)*x + 96.D0*Hr4(1,0,-1,0)*dx + 48.D0*
c$$$     &    Hr4(1,0,0,0) - 24.D0*Hr4(1,0,0,0)*x - 48.D0*Hr4(1,0,0,0)*dx
c$$$     &     - 224.D0*Hr4(1,0,0,1) + 112.D0*Hr4(1,0,0,1)*x + 224.D0*Hr4(1
c$$$     &    ,0,0,1)*dx )
c$$$      PTgq2 = PTgq2 + cf**2*ca * ( 544.D0*Hr4(1,0,1,0) - 272.D0*Hr4(1,0
c$$$     &    ,1,0)*x - 544.D0*Hr4(1,0,1,0)*dx - 64.D0*Hr4(1,0,1,1) + 32.D0
c$$$     &    *Hr4(1,0,1,1)*x + 64.D0*Hr4(1,0,1,1)*dx - 160.D0*Hr4(1,1,0,0)
c$$$     &     + 80.D0*Hr4(1,1,0,0)*x + 160.D0*Hr4(1,1,0,0)*dx + 32.D0*Hr4(
c$$$     &    1,1,0,1) - 16.D0*Hr4(1,1,0,1)*x - 32.D0*Hr4(1,1,0,1)*dx + 224.
c$$$     &    D0*Hr4(1,1,1,0) - 112.D0*Hr4(1,1,1,0)*x - 224.D0*Hr4(1,1,1,0)
c$$$     &    *dx + 128.D0*Hr4(1,1,1,1) - 64.D0*Hr4(1,1,1,1)*x - 128.D0*
c$$$     &    Hr4(1,1,1,1)*dx )
c$$$      PTgq2 = PTgq2 + cf**3 * (  - 1915.D0/2.D0 + 731.D0/4.D0*x + 794.D0
c$$$     &    *dx - 96.D0*z3 - 196.D0*z3*x - 720.D0*z3*dx - 224.D0*z2 + 148.
c$$$     &    D0*z2*x - 256.D0*z2*dx - 776.D0/5.D0*z2**2 + 84.D0/5.D0*z2**2
c$$$     &    *x + 1216.D0/5.D0*z2**2*dx - 64.D0*Hr1(-1)*z2 - 32.D0*Hr1(-1)
c$$$     &    *z2*x + 197.D0*Hr1(0) + 29.D0*Hr1(0)*x + 208.D0*Hr1(0)*dx - 
c$$$     &    80.D0*Hr1(0)*z3 - 200.D0*Hr1(0)*z3*x + 128.D0*Hr1(0)*z3*dx - 
c$$$     &    24.D0*Hr1(0)*z2 - 20.D0*Hr1(0)*z2*x - 192.D0*Hr1(0)*z2*dx - 
c$$$     &    230.D0*Hr1(1) + 16.D0*Hr1(1)*x + 228.D0*Hr1(1)*dx - 64.D0*
c$$$     &    Hr1(1)*z3 + 32.D0*Hr1(1)*z3*x + 64.D0*Hr1(1)*z3*dx - 352.D0*
c$$$     &    Hr1(1)*z2 + 252.D0*Hr1(1)*z2*x + 144.D0*Hr1(1)*z2*dx - 80.D0*
c$$$     &    Hr2(-1,0) + 208.D0*Hr2(-1,0)*x - 256.D0*Hr2(-1,0)*dx + 128.D0
c$$$     &    *Hr2(0,-1)*z2 + 96.D0*Hr2(0,-1)*z2*x + 128.D0*Hr2(0,-1)*z2*dx
c$$$     &     - 290.D0*Hr2(0,0) - 251.D0*Hr2(0,0)*x - 64.D0*Hr2(0,0)*z2 - 
c$$$     &    64.D0*Hr2(0,0)*z2*x + 208.D0*Hr2(0,1) - 140.D0*Hr2(0,1)*x - 
c$$$     &    464.D0*Hr2(0,1)*z2 + 200.D0*Hr2(0,1)*z2*x + 320.D0*Hr2(0,1)*
c$$$     &    z2*dx )
c$$$      PTgq2 = PTgq2 + cf**3 * (  - 264.D0*Hr2(1,0) + 140.D0*Hr2(1,0)*x
c$$$     &     + 180.D0*Hr2(1,0)*dx - 128.D0*Hr2(1,0)*z2 + 64.D0*Hr2(1,0)*
c$$$     &    z2*x + 128.D0*Hr2(1,0)*z2*dx - 72.D0*Hr2(1,1) + 80.D0*Hr2(1,1
c$$$     &    )*x + 60.D0*Hr2(1,1)*dx - 96.D0*Hr2(1,1)*z2 + 48.D0*Hr2(1,1)*
c$$$     &    z2*x + 96.D0*Hr2(1,1)*z2*dx - 128.D0*Hr3(-1,-1,0) - 64.D0*
c$$$     &    Hr3(-1,-1,0)*x + 224.D0*Hr3(-1,0,0) + 96.D0*Hr3(-1,0,0)*x + 
c$$$     &    96.D0*Hr3(-1,0,0)*dx + 128.D0*Hr3(0,-1,0) - 112.D0*Hr3(0,-1,0
c$$$     &    )*x - 192.D0*Hr3(0,-1,0)*dx + 160.D0*Hr3(0,0,0) + 40.D0*Hr3(0
c$$$     &    ,0,0)*x + 248.D0*Hr3(0,0,1) - 28.D0*Hr3(0,0,1)*x - 104.D0*
c$$$     &    Hr3(0,1,0) + 40.D0*Hr3(0,1,0)*x - 72.D0*Hr3(0,1,1) + 56.D0*
c$$$     &    Hr3(0,1,1)*x + 128.D0*Hr3(1,0,0) - 36.D0*Hr3(1,0,0)*x + 48.D0
c$$$     &    *Hr3(1,0,0)*dx + 128.D0*Hr3(1,0,1) - 60.D0*Hr3(1,0,1)*x - 48.D
c$$$     &    0*Hr3(1,0,1)*dx - 32.D0*Hr3(1,1,0) + 4.D0*Hr3(1,1,0)*x + 48.D0
c$$$     &    *Hr3(1,1,0)*dx + 16.D0*Hr3(1,1,1) - 52.D0*Hr3(1,1,1)*x - 128.D
c$$$     &    0*Hr4(-1,-1,0,0) - 64.D0*Hr4(-1,-1,0,0)*x - 128.D0*Hr4(-1,-1,
c$$$     &    0,0)*dx )
c$$$      PTgq2 = PTgq2 + cf**3 * (  - 64.D0*Hr4(-1,0,-1,0) - 32.D0*Hr4(-1,
c$$$     &    0,-1,0)*x - 64.D0*Hr4(-1,0,-1,0)*dx + 32.D0*Hr4(-1,0,0,0) + 
c$$$     &    16.D0*Hr4(-1,0,0,0)*x + 32.D0*Hr4(-1,0,0,0)*dx + 256.D0*Hr4(0
c$$$     &    ,-1,-1,0) + 192.D0*Hr4(0,-1,-1,0)*x + 256.D0*Hr4(0,-1,-1,0)*
c$$$     &    dx - 64.D0*Hr4(0,-1,0,0) - 64.D0*Hr4(0,-1,0,0)*x - 128.D0*
c$$$     &    Hr4(0,-1,0,0)*dx - 64.D0*Hr4(0,0,-1,0) - 96.D0*Hr4(0,0,-1,0)*
c$$$     &    x - 32.D0*Hr4(0,0,0,0) + 96.D0*Hr4(0,0,0,0)*x - 128.D0*Hr4(0,
c$$$     &    0,0,1) + 64.D0*Hr4(0,0,0,1)*x - 176.D0*Hr4(0,0,1,0) + 88.D0*
c$$$     &    Hr4(0,0,1,0)*x + 128.D0*Hr4(0,0,1,0)*dx - 240.D0*Hr4(0,0,1,1)
c$$$     &     + 120.D0*Hr4(0,0,1,1)*x + 256.D0*Hr4(0,0,1,1)*dx - 16.D0*
c$$$     &    Hr4(0,1,0,0) + 40.D0*Hr4(0,1,0,0)*x + 128.D0*Hr4(0,1,0,0)*dx
c$$$     &     + 144.D0*Hr4(0,1,0,1) - 72.D0*Hr4(0,1,0,1)*x - 64.D0*Hr4(0,1
c$$$     &    ,0,1)*dx + 80.D0*Hr4(0,1,1,0) - 40.D0*Hr4(0,1,1,0)*x - 64.D0*
c$$$     &    Hr4(0,1,1,0)*dx + 208.D0*Hr4(0,1,1,1) - 104.D0*Hr4(0,1,1,1)*x
c$$$     &     - 192.D0*Hr4(0,1,1,1)*dx + 64.D0*Hr4(1,0,-1,0) - 32.D0*Hr4(1
c$$$     &    ,0,-1,0)*x )
c$$$      PTgq2 = PTgq2 + cf**3 * (  - 64.D0*Hr4(1,0,-1,0)*dx + 192.D0*Hr4(
c$$$     &    1,0,0,0) - 96.D0*Hr4(1,0,0,0)*x - 192.D0*Hr4(1,0,0,0)*dx + 
c$$$     &    256.D0*Hr4(1,0,0,1) - 128.D0*Hr4(1,0,0,1)*x - 256.D0*Hr4(1,0,
c$$$     &    0,1)*dx - 32.D0*Hr4(1,0,1,0) + 16.D0*Hr4(1,0,1,0)*x + 32.D0*
c$$$     &    Hr4(1,0,1,0)*dx + 128.D0*Hr4(1,0,1,1) - 64.D0*Hr4(1,0,1,1)*x
c$$$     &     - 128.D0*Hr4(1,0,1,1)*dx + 192.D0*Hr4(1,1,0,0) - 96.D0*Hr4(1
c$$$     &    ,1,0,0)*x - 192.D0*Hr4(1,1,0,0)*dx + 96.D0*Hr4(1,1,0,1) - 48.D
c$$$     &    0*Hr4(1,1,0,1)*x - 96.D0*Hr4(1,1,0,1)*dx - 32.D0*Hr4(1,1,1,0)
c$$$     &     + 16.D0*Hr4(1,1,1,0)*x + 32.D0*Hr4(1,1,1,0)*dx - 64.D0*Hr4(1
c$$$     &    ,1,1,1) + 32.D0*Hr4(1,1,1,1)*x + 64.D0*Hr4(1,1,1,1)*dx )
c$$$      PTgq2 = PTgq2 + nf*cf*ca * ( 332.D0/27.D0 - 2624.D0/27.D0*x + 
c$$$     &    1892.D0/81.D0*x**2 + 3778.D0/81.D0*dx + 64.D0*z3 - 40.D0*z3*x
c$$$     &     - 80.D0*z3*dx + 8.D0/3.D0*z2 - 176.D0/9.D0*z2*x - 160.D0/9.D0
c$$$     &    *z2*dx + 1280.D0/27.D0*Hr1(0) + 308.D0/27.D0*Hr1(0)*x + 64.D0/
c$$$     &    9.D0*Hr1(0)*x**2 + 1156.D0/27.D0*Hr1(0)*dx - 16.D0/3.D0*Hr1(0
c$$$     &    )*z2*x - 32.D0/3.D0*Hr1(0)*z2*dx + 448.D0/27.D0*Hr1(1) - 392.D
c$$$     &    0/27.D0*Hr1(1)*x - 448.D0/27.D0*Hr1(1)*dx - 32.D0*Hr1(1)*z2
c$$$     &     + 16.D0*Hr1(1)*z2*x + 32.D0*Hr1(1)*z2*dx - 160.D0/9.D0*Hr2(
c$$$     &    -1,0) - 104.D0/9.D0*Hr2(-1,0)*x - 160.D0/9.D0*Hr2(-1,0)*dx - 
c$$$     &    256.D0/9.D0*Hr2(0,0) - 464.D0/9.D0*Hr2(0,0)*x - 128.D0/9.D0*
c$$$     &    Hr2(0,0)*x**2 + 656.D0/9.D0*Hr2(0,0)*dx + 136.D0/9.D0*Hr2(0,1
c$$$     &    ) - 104.D0/9.D0*Hr2(0,1)*x - 160.D0/9.D0*Hr2(0,1)*dx - 320.D0/
c$$$     &    9.D0*Hr2(1,0) + 136.D0/9.D0*Hr2(1,0)*x + 320.D0/9.D0*Hr2(1,0)
c$$$     &    *dx - 160.D0/9.D0*Hr2(1,1) + 104.D0/9.D0*Hr2(1,1)*x + 160.D0/
c$$$     &    9.D0*Hr2(1,1)*dx - 32.D0/3.D0*Hr3(-1,0,0) - 16.D0/3.D0*Hr3(-1
c$$$     &    ,0,0)*x )
c$$$      PTgq2 = PTgq2 + nf*cf*ca * (  - 32.D0/3.D0*Hr3(-1,0,0)*dx - 32.D0/
c$$$     &    3.D0*Hr3(0,-1,0) - 16.D0/3.D0*Hr3(0,-1,0)*x - 32.D0/3.D0*Hr3(
c$$$     &    0,-1,0)*dx + 32.D0/3.D0*Hr3(0,0,0) + 128.D0/3.D0*Hr3(0,0,0)*x
c$$$     &     + 128.D0/3.D0*Hr3(0,0,0)*dx + 32.D0/3.D0*Hr3(0,0,1) - 16.D0/
c$$$     &    3.D0*Hr3(0,0,1)*x - 32.D0/3.D0*Hr3(0,0,1)*dx - 16.D0/3.D0*
c$$$     &    Hr3(0,1,1) + 8.D0/3.D0*Hr3(0,1,1)*x + 16.D0/3.D0*Hr3(0,1,1)*
c$$$     &    dx - 160.D0/3.D0*Hr3(1,0,0) + 80.D0/3.D0*Hr3(1,0,0)*x + 160.D0
c$$$     &    /3.D0*Hr3(1,0,0)*dx - 32.D0/3.D0*Hr3(1,1,0) + 16.D0/3.D0*Hr3(
c$$$     &    1,1,0)*x + 32.D0/3.D0*Hr3(1,1,0)*dx + 16.D0/3.D0*Hr3(1,1,1)
c$$$     &     - 8.D0/3.D0*Hr3(1,1,1)*x - 16.D0/3.D0*Hr3(1,1,1)*dx )
c$$$      PTgq2 = PTgq2 + nf*cf**2 * (  - 12803.D0/27.D0 + 30475.D0/54.D0*x
c$$$     &     - 21784.D0/81.D0*x**2 + 14008.D0/81.D0*dx - 2330.D0/27.D0*
c$$$     &    Hr1(0) + 2494.D0/27.D0*Hr1(0)*x + 3040.D0/27.D0*Hr1(0)*x**2
c$$$     &     + 200.D0/9.D0*Hr1(0)*dx - 484.D0/27.D0*Hr1(1) + 464.D0/27.D0
c$$$     &    *Hr1(1)*x + 448.D0/27.D0*Hr1(1)*dx + 16.D0*Hr1(1)*z2 - 8.D0*
c$$$     &    Hr1(1)*z2*x - 16.D0*Hr1(1)*z2*dx - 204.D0*Hr2(0,0) + 142.D0/3.
c$$$     &    D0*Hr2(0,0)*x - 256.D0/9.D0*Hr2(0,0)*x**2 - 128.D0/3.D0*Hr2(0
c$$$     &    ,0)*dx - 320.D0/9.D0*Hr2(0,1) + 208.D0/9.D0*Hr2(0,1)*x + 320.D
c$$$     &    0/9.D0*Hr2(0,1)*dx + 80.D0/3.D0*Hr2(1,0) - 40.D0/3.D0*Hr2(1,0
c$$$     &    )*x - 80.D0/3.D0*Hr2(1,0)*dx + 160.D0/9.D0*Hr2(1,1) - 104.D0/
c$$$     &    9.D0*Hr2(1,1)*x - 160.D0/9.D0*Hr2(1,1)*dx - 176.D0/3.D0*Hr3(0
c$$$     &    ,0,0) + 160.D0/3.D0*Hr3(0,0,0)*x - 128.D0/3.D0*Hr3(0,0,0)*dx
c$$$     &     - 64.D0/3.D0*Hr3(0,0,1) + 32.D0/3.D0*Hr3(0,0,1)*x + 64.D0/3.D
c$$$     &    0*Hr3(0,0,1)*dx + 32.D0/3.D0*Hr3(0,1,1) - 16.D0/3.D0*Hr3(0,1,
c$$$     &    1)*x - 32.D0/3.D0*Hr3(0,1,1)*dx + 80.D0/3.D0*Hr3(1,0,0) - 40.D
c$$$     &    0/3.D0*Hr3(1,0,0)*x )
c$$$      PTgq2 = PTgq2 + nf*cf**2 * (  - 80.D0/3.D0*Hr3(1,0,0)*dx + 16.D0/
c$$$     &    3.D0*Hr3(1,0,1) - 8.D0/3.D0*Hr3(1,0,1)*x - 16.D0/3.D0*Hr3(1,0
c$$$     &    ,1)*dx + 16.D0/3.D0*Hr3(1,1,0) - 8.D0/3.D0*Hr3(1,1,0)*x - 16.D
c$$$     &    0/3.D0*Hr3(1,1,0)*dx - 16.D0/3.D0*Hr3(1,1,1) + 8.D0/3.D0*Hr3(
c$$$     &    1,1,1)*x + 16.D0/3.D0*Hr3(1,1,1)*dx + 64.D0*Hr4(0,0,0,0) - 32.
c$$$     &    D0*Hr4(0,0,0,0)*x )
c$$$*
c$$$       X2GQTA = PTgq2
c$$$*
c$$$       RETURN
c$$$       END
c$$$* ---------------------------------------------------------------------
c$$$*
c$$$*
c$$$* ..The regular piece of the gluon-gluon splitting functions P_gg^(2)T
c$$$*
c$$$       FUNCTION X2GGTA (X, NF)
c$$$*
c$$$       IMPLICIT REAL*8 (A - Z)
c$$$       COMPLEX*16 HC1, HC2, HC3, HC4
c$$$       INTEGER NF, NF2, N1, N2, NW, I1, I2, I3, N
c$$$       PARAMETER ( N1 = -1, N2 = 1, NW = 4 )
c$$$       DIMENSION HC1(N1:N2),HC2(N1:N2,N1:N2),HC3(N1:N2,N1:N2,N1:N2),
c$$$     ,           HC4(N1:N2,N1:N2,N1:N2,N1:N2)
c$$$       DIMENSION HR1(N1:N2),HR2(N1:N2,N1:N2),HR3(N1:N2,N1:N2,N1:N2),
c$$$     ,           HR4(N1:N2,N1:N2,N1:N2,N1:N2)
c$$$       DIMENSION HI1(N1:N2),HI2(N1:N2,N1:N2),HI3(N1:N2,N1:N2,N1:N2),
c$$$     ,           HI4(N1:N2,N1:N2,N1:N2,N1:N2)
c$$$       PARAMETER ( Z2 = 1.6449 34066 84822 64365 D0,
c$$$     ,             Z3 = 1.2020 56903 15959 42854 D0,
c$$$     ,             Z5 = 1.0369 27755 14336 99263 D0 )
c$$$*
c$$$* ..The soft coefficient for use in X2GGTB and X2GGTC
c$$$*
c$$$       COMMON / P2GSOFT / A3G
c$$$*
c$$$* ...Colour factors
c$$$*
c$$$       CF  = 4./3.D0
c$$$       CA  = 3.D0
c$$$       NF2 = NF*NF
c$$$*
c$$$* ...Some abbreviations
c$$$*
c$$$       DX = 1.D0/X
c$$$       DM = 1.D0/(1.D0-X)
c$$$       DP = 1.D0/(1.D0+X)
c$$$*
c$$$* ...The harmonic polylogs up to weight 4 by Gehrmann and Remiddi
c$$$*
c$$$       CALL HPLOG (X, NW, HC1,HC2,HC3,HC4, HR1,HR2,HR3,HR4,
c$$$     ,            HI1,HI2,HI3,HI4, N1, N2)
c$$$*
c$$$* ...The splitting function in terms of the harmonic polylogs
c$$$*    (without the delta(1-x) part, but with the soft contribution)
c$$$*
c$$$      PTgg2 =
c$$$     &  + ca**3 * (  - 54088.D0/27.D0 + 49678.D0/27.D0*x - 146182.D0/81.
c$$$     &    D0*x**2 + 146182.D0/81.D0*dx + 490.D0/3.D0*dm + 1064.D0/3.D0*
c$$$     &    z3 + 704.D0/3.D0*z3*x + 1408.D0/3.D0*z3*x**2 - 1144.D0/3.D0*
c$$$     &    z3*dx + 176.D0/3.D0*z3*dm + 976.D0/3.D0*z2 - 2632.D0/9.D0*z2*
c$$$     &    x + 1072.D0/9.D0*z2*x**2 - 3112.D0/9.D0*z2*dx - 1072.D0/9.D0*
c$$$     &    z2*dp - 1072.D0/9.D0*z2*dm - 272.D0*z2**2 - 1152.D0/5.D0*
c$$$     &    z2**2*x - 416.D0/5.D0*z2**2*x**2 - 464.D0/5.D0*z2**2*dx + 88.D
c$$$     &    0*z2**2*dp - 24.D0/5.D0*z2**2*dm - 384.D0*Hr1(-1)*z3 - 192.D0
c$$$     &    *Hr1(-1)*z3*x - 192.D0*Hr1(-1)*z3*x**2 - 192.D0*Hr1(-1)*z3*dx
c$$$     &     + 192.D0*Hr1(-1)*z3*dp - 424.D0*Hr1(-1)*z2 - 424.D0*Hr1(-1)*
c$$$     &    z2*x - 88.D0*Hr1(-1)*z2*x**2 - 88.D0*Hr1(-1)*z2*dx + 1274.D0/
c$$$     &    9.D0*Hr1(0) + 9382.D0/9.D0*Hr1(0)*x + 5524.D0/27.D0*Hr1(0)*
c$$$     &    x**2 + 32608.D0/27.D0*Hr1(0)*dx + 8272.D0/27.D0*Hr1(0)*dm - 
c$$$     &    800.D0*Hr1(0)*z3*x - 288.D0*Hr1(0)*z3*dx + 144.D0*Hr1(0)*z3*
c$$$     &    dp - 144.D0*Hr1(0)*z3*dm + 808.D0/3.D0*Hr1(0)*z2 + 512.D0/3.D0
c$$$     &    *Hr1(0)*z2*x )
c$$$      PTgg2 = PTgg2 + ca**3 * ( 176.D0/3.D0*Hr1(0)*z2*x**2 - 880.D0/3.D0
c$$$     &    *Hr1(0)*z2*dx - 88.D0*Hr1(0)*z2*dp - 176.D0/3.D0*Hr1(0)*z2*dm
c$$$     &     + 124.D0/9.D0*Hr1(1) - 124.D0/9.D0*Hr1(1)*x + 1652.D0/27.D0*
c$$$     &    Hr1(1)*x**2 - 1652.D0/27.D0*Hr1(1)*dx + 192.D0*Hr1(1)*z3 - 96.
c$$$     &    D0*Hr1(1)*z3*x + 96.D0*Hr1(1)*z3*x**2 - 96.D0*Hr1(1)*z3*dx - 
c$$$     &    96.D0*Hr1(1)*z3*dm - 24.D0*Hr1(1)*z2 + 24.D0*Hr1(1)*z2*x + 88.
c$$$     &    D0*Hr1(1)*z2*x**2 - 88.D0*Hr1(1)*z2*dx + 512.D0*Hr2(-1,-1)*z2
c$$$     &     + 256.D0*Hr2(-1,-1)*z2*x + 256.D0*Hr2(-1,-1)*z2*x**2 + 256.D0
c$$$     &    *Hr2(-1,-1)*z2*dx - 256.D0*Hr2(-1,-1)*z2*dp - 872.D0/9.D0*
c$$$     &    Hr2(-1,0) - 3016.D0/9.D0*Hr2(-1,0)*x - 680.D0/3.D0*Hr2(-1,0)*
c$$$     &    x**2 - 680.D0/3.D0*Hr2(-1,0)*dx - 2144.D0/9.D0*Hr2(-1,0)*dp
c$$$     &     - 64.D0*Hr2(-1,0)*z2 - 32.D0*Hr2(-1,0)*z2*x - 32.D0*Hr2(-1,0
c$$$     &    )*z2*x**2 - 32.D0*Hr2(-1,0)*z2*dx + 32.D0*Hr2(-1,0)*z2*dp - 
c$$$     &    256.D0*Hr2(0,-1)*z2 + 160.D0*Hr2(0,-1)*z2*x - 96.D0*Hr2(0,-1)
c$$$     &    *z2*x**2 + 96.D0*Hr2(0,-1)*z2*dx + 96.D0*Hr2(0,-1)*z2*dm - 
c$$$     &    116.D0/9.D0*Hr2(0,0) )
c$$$      PTgg2 = PTgg2 + ca**3 * ( 9044.D0/9.D0*Hr2(0,0)*x + 2248.D0/9.D0*
c$$$     &    Hr2(0,0)*x**2 + 6224.D0/9.D0*Hr2(0,0)*dx + 1072.D0/9.D0*Hr2(0
c$$$     &    ,0)*dp - 1280.D0/9.D0*Hr2(0,0)*dm + 96.D0*Hr2(0,0)*z2 - 608.D0
c$$$     &    *Hr2(0,0)*z2*x - 256.D0*Hr2(0,0)*z2*dx + 128.D0*Hr2(0,0)*z2*
c$$$     &    dp - 128.D0*Hr2(0,0)*z2*dm + 1760.D0/9.D0*Hr2(0,1) - 784.D0/9.
c$$$     &    D0*Hr2(0,1)*x + 1072.D0/9.D0*Hr2(0,1)*x**2 - 1072.D0/9.D0*
c$$$     &    Hr2(0,1)*dx - 2144.D0/9.D0*Hr2(0,1)*dm - 256.D0*Hr2(0,1)*z2
c$$$     &     - 32.D0*Hr2(0,1)*z2*x - 96.D0*Hr2(0,1)*z2*x**2 + 32.D0*Hr2(0
c$$$     &    ,1)*z2*dx + 32.D0*Hr2(0,1)*z2*dp + 64.D0*Hr2(0,1)*z2*dm + 
c$$$     &    2344.D0/9.D0*Hr2(1,0) - 200.D0/9.D0*Hr2(1,0)*x + 1072.D0/9.D0
c$$$     &    *Hr2(1,0)*x**2 - 1072.D0/9.D0*Hr2(1,0)*dx - 2144.D0/9.D0*Hr2(
c$$$     &    1,0)*dm - 64.D0*Hr2(1,0)*z2 + 32.D0*Hr2(1,0)*z2*x - 32.D0*
c$$$     &    Hr2(1,0)*z2*x**2 + 32.D0*Hr2(1,0)*z2*dx + 32.D0*Hr2(1,0)*z2*
c$$$     &    dm - 48.D0*Hr3(-1,-1,0) - 48.D0*Hr3(-1,-1,0)*x + 176.D0*Hr3(
c$$$     &    -1,-1,0)*x**2 + 176.D0*Hr3(-1,-1,0)*dx + 928.D0/3.D0*Hr3(-1,0
c$$$     &    ,0) )
c$$$      PTgg2 = PTgg2 + ca**3 * ( 224.D0/3.D0*Hr3(-1,0,0)*x - 704.D0/3.D0
c$$$     &    *Hr3(-1,0,0)*x**2 - 704.D0/3.D0*Hr3(-1,0,0)*dx - 704.D0/3.D0*
c$$$     &    Hr3(-1,0,0)*dp + 400.D0*Hr3(-1,0,1) + 400.D0*Hr3(-1,0,1)*x + 
c$$$     &    176.D0*Hr3(-1,0,1)*x**2 + 176.D0*Hr3(-1,0,1)*dx + 64.D0*Hr3(0
c$$$     &    ,-1,0) - 320.D0/3.D0*Hr3(0,-1,0)*x - 1232.D0/3.D0*Hr3(0,-1,0)
c$$$     &    *dx - 352.D0/3.D0*Hr3(0,-1,0)*dp + 176.D0/3.D0*Hr3(0,-1,0)*dm
c$$$     &     + 176.D0*Hr3(0,0,0) - 1264.D0/3.D0*Hr3(0,0,0)*x + 352.D0*
c$$$     &    Hr3(0,0,0)*x**2 + 704.D0*Hr3(0,0,0)*dx + 176.D0*Hr3(0,0,0)*dp
c$$$     &     - 528.D0*Hr3(0,0,0)*dm + 376.D0/3.D0*Hr3(0,0,1) - 512.D0/3.D0
c$$$     &    *Hr3(0,0,1)*x + 176.D0/3.D0*Hr3(0,0,1)*x**2 + 704.D0/3.D0*
c$$$     &    Hr3(0,0,1)*dx - 616.D0/3.D0*Hr3(0,0,1)*dm + 944.D0/3.D0*Hr3(0
c$$$     &    ,1,0) - 16.D0/3.D0*Hr3(0,1,0)*x + 352.D0/3.D0*Hr3(0,1,0)*x**2
c$$$     &     + 352.D0/3.D0*Hr3(0,1,0)*dx - 704.D0/3.D0*Hr3(0,1,0)*dm + 
c$$$     &    632.D0/3.D0*Hr3(1,0,0) - 16.D0/3.D0*Hr3(1,0,0)*x - 176.D0/3.D0
c$$$     &    *Hr3(1,0,0)*x**2 + 176.D0/3.D0*Hr3(1,0,0)*dx - 616.D0/3.D0*
c$$$     &    Hr3(1,0,0)*dm )
c$$$      PTgg2 = PTgg2 + ca**3 * ( 256.D0*Hr4(-1,-1,0,0) + 128.D0*Hr4(-1,
c$$$     &    -1,0,0)*x + 128.D0*Hr4(-1,-1,0,0)*x**2 + 128.D0*Hr4(-1,-1,0,0
c$$$     &    )*dx - 128.D0*Hr4(-1,-1,0,0)*dp - 512.D0*Hr4(-1,-1,0,1) - 256.
c$$$     &    D0*Hr4(-1,-1,0,1)*x - 256.D0*Hr4(-1,-1,0,1)*x**2 - 256.D0*
c$$$     &    Hr4(-1,-1,0,1)*dx + 256.D0*Hr4(-1,-1,0,1)*dp + 256.D0*Hr4(-1,
c$$$     &    0,-1,0) + 128.D0*Hr4(-1,0,-1,0)*x + 128.D0*Hr4(-1,0,-1,0)*
c$$$     &    x**2 + 128.D0*Hr4(-1,0,-1,0)*dx - 128.D0*Hr4(-1,0,-1,0)*dp - 
c$$$     &    640.D0*Hr4(-1,0,0,0) - 320.D0*Hr4(-1,0,0,0)*x - 320.D0*Hr4(-1
c$$$     &    ,0,0,0)*x**2 - 320.D0*Hr4(-1,0,0,0)*dx + 320.D0*Hr4(-1,0,0,0)
c$$$     &    *dp - 128.D0*Hr4(-1,0,1,0) - 64.D0*Hr4(-1,0,1,0)*x - 64.D0*
c$$$     &    Hr4(-1,0,1,0)*x**2 - 64.D0*Hr4(-1,0,1,0)*dx + 64.D0*Hr4(-1,0,
c$$$     &    1,0)*dp + 320.D0*Hr4(0,-1,-1,0)*x + 64.D0*Hr4(0,-1,-1,0)*x**2
c$$$     &     + 192.D0*Hr4(0,-1,-1,0)*dx - 128.D0*Hr4(0,-1,-1,0)*dp + 64.D0
c$$$     &    *Hr4(0,-1,-1,0)*dm + 128.D0*Hr4(0,-1,0,0) - 544.D0*Hr4(0,-1,0
c$$$     &    ,0)*x - 96.D0*Hr4(0,-1,0,0)*x**2 - 224.D0*Hr4(0,-1,0,0)*dx + 
c$$$     &    160.D0*Hr4(0,-1,0,0)*dp )
c$$$      PTgg2 = PTgg2 + ca**3 * (  - 64.D0*Hr4(0,-1,0,0)*dm + 256.D0*Hr4(
c$$$     &    0,-1,0,1) + 128.D0*Hr4(0,-1,0,1)*x**2 - 64.D0*Hr4(0,-1,0,1)*
c$$$     &    dp - 64.D0*Hr4(0,-1,0,1)*dm - 64.D0*Hr4(0,0,-1,0) - 704.D0*
c$$$     &    Hr4(0,0,-1,0)*x - 128.D0*Hr4(0,0,-1,0)*x**2 - 256.D0*Hr4(0,0,
c$$$     &    -1,0)*dx + 192.D0*Hr4(0,0,-1,0)*dp - 64.D0*Hr4(0,0,-1,0)*dm
c$$$     &     - 256.D0*Hr4(0,0,0,0) + 1984.D0*Hr4(0,0,0,0)*x - 128.D0*Hr4(
c$$$     &    0,0,0,0)*x**2 + 512.D0*Hr4(0,0,0,0)*dx - 192.D0*Hr4(0,0,0,0)*
c$$$     &    dp + 320.D0*Hr4(0,0,0,0)*dm - 96.D0*Hr4(0,0,0,1) + 800.D0*
c$$$     &    Hr4(0,0,0,1)*x - 256.D0*Hr4(0,0,0,1)*x**2 + 384.D0*Hr4(0,0,0,
c$$$     &    1)*dx - 64.D0*Hr4(0,0,0,1)*dp + 320.D0*Hr4(0,0,0,1)*dm - 128.D
c$$$     &    0*Hr4(0,0,1,0) + 512.D0*Hr4(0,0,1,0)*x - 192.D0*Hr4(0,0,1,0)*
c$$$     &    x**2 + 256.D0*Hr4(0,0,1,0)*dx - 32.D0*Hr4(0,0,1,0)*dp + 224.D0
c$$$     &    *Hr4(0,0,1,0)*dm - 256.D0*Hr4(0,0,1,1) + 128.D0*Hr4(0,0,1,1)*
c$$$     &    x - 128.D0*Hr4(0,0,1,1)*x**2 + 128.D0*Hr4(0,0,1,1)*dx + 128.D0
c$$$     &    *Hr4(0,0,1,1)*dm - 288.D0*Hr4(0,1,0,0) + 384.D0*Hr4(0,1,0,0)*
c$$$     &    x )
c$$$      PTgg2 = PTgg2 + ca**3 * (  - 224.D0*Hr4(0,1,0,0)*x**2 + 224.D0*
c$$$     &    Hr4(0,1,0,0)*dx + 224.D0*Hr4(0,1,0,0)*dm - 256.D0*Hr4(0,1,0,1
c$$$     &    ) + 128.D0*Hr4(0,1,0,1)*x - 128.D0*Hr4(0,1,0,1)*x**2 + 128.D0
c$$$     &    *Hr4(0,1,0,1)*dx + 128.D0*Hr4(0,1,0,1)*dm - 256.D0*Hr4(0,1,1,
c$$$     &    0) + 128.D0*Hr4(0,1,1,0)*x - 128.D0*Hr4(0,1,1,0)*x**2 + 128.D0
c$$$     &    *Hr4(0,1,1,0)*dx + 128.D0*Hr4(0,1,1,0)*dm + 128.D0*Hr4(1,0,-1
c$$$     &    ,0) - 64.D0*Hr4(1,0,-1,0)*x + 64.D0*Hr4(1,0,-1,0)*x**2 - 64.D0
c$$$     &    *Hr4(1,0,-1,0)*dx - 64.D0*Hr4(1,0,-1,0)*dm - 640.D0*Hr4(1,0,0
c$$$     &    ,0) + 320.D0*Hr4(1,0,0,0)*x - 320.D0*Hr4(1,0,0,0)*x**2 + 320.D
c$$$     &    0*Hr4(1,0,0,0)*dx + 320.D0*Hr4(1,0,0,0)*dm - 256.D0*Hr4(1,0,0
c$$$     &    ,1) + 128.D0*Hr4(1,0,0,1)*x - 128.D0*Hr4(1,0,0,1)*x**2 + 128.D
c$$$     &    0*Hr4(1,0,0,1)*dx + 128.D0*Hr4(1,0,0,1)*dm - 256.D0*Hr4(1,0,1
c$$$     &    ,0) + 128.D0*Hr4(1,0,1,0)*x - 128.D0*Hr4(1,0,1,0)*x**2 + 128.D
c$$$     &    0*Hr4(1,0,1,0)*dx + 128.D0*Hr4(1,0,1,0)*dm - 256.D0*Hr4(1,1,0
c$$$     &    ,0) + 128.D0*Hr4(1,1,0,0)*x - 128.D0*Hr4(1,1,0,0)*x**2 + 128.D
c$$$     &    0*Hr4(1,1,0,0)*dx )
c$$$      PTgg2 = PTgg2 + ca**3 * ( 128.D0*Hr4(1,1,0,0)*dm )
c$$$      PTgg2 = PTgg2 + nf*ca**2 * (  - 2174.D0/9.D0 + 7358.D0/27.D0*x - 
c$$$     &    19264.D0/81.D0*x**2 + 19264.D0/81.D0*dx - 836.D0/27.D0*dm + 
c$$$     &    496.D0/3.D0*z3 + 112.D0/3.D0*z3*x + 128.D0/3.D0*z3*x**2 - 176.
c$$$     &    D0/3.D0*z3*dx - 128.D0/3.D0*z3*dm - 1700.D0/9.D0*z2 + 1076.D0/
c$$$     &    9.D0*z2*x - 176.D0/3.D0*z2*x**2 - 128.D0/3.D0*z2*dx + 160.D0/
c$$$     &    9.D0*z2*dp + 160.D0/9.D0*z2*dm + 276.D0/5.D0*z2**2 + 156.D0/5.
c$$$     &    D0*z2**2*x - 88.D0*Hr1(-1)*z2 - 88.D0*Hr1(-1)*z2*x + 16.D0/3.D
c$$$     &    0*Hr1(-1)*z2*x**2 + 16.D0/3.D0*Hr1(-1)*z2*dx + 1408.D0/9.D0*
c$$$     &    Hr1(0) + 202.D0/3.D0*Hr1(0)*x + 340.D0/3.D0*Hr1(0)*x**2 + 
c$$$     &    5912.D0/27.D0*Hr1(0)*dx - 2600.D0/27.D0*Hr1(0)*dm - 96.D0*
c$$$     &    Hr1(0)*z3*x - 196.D0/3.D0*Hr1(0)*z2 + 76.D0/3.D0*Hr1(0)*z2*x
c$$$     &     - 32.D0/3.D0*Hr1(0)*z2*x**2 - 32.D0/3.D0*Hr1(0)*z2*dx + 16.D0
c$$$     &    *Hr1(0)*z2*dp + 32.D0/3.D0*Hr1(0)*z2*dm + 1354.D0/9.D0*Hr1(1)
c$$$     &     - 1354.D0/9.D0*Hr1(1)*x - 820.D0/27.D0*Hr1(1)*x**2 + 820.D0/
c$$$     &    27.D0*Hr1(1)*dx - 56.D0*Hr1(1)*z2 + 56.D0*Hr1(1)*z2*x + 16.D0
c$$$     &    *Hr1(1)*z2*x**2 )
c$$$      PTgg2 = PTgg2 + nf*ca**2 * (  - 16.D0*Hr1(1)*z2*dx + 1208.D0/9.D0
c$$$     &    *Hr2(-1,0) + 1528.D0/9.D0*Hr2(-1,0)*x - 304.D0/3.D0*Hr2(-1,0)
c$$$     &    *x**2 - 304.D0/3.D0*Hr2(-1,0)*dx + 320.D0/9.D0*Hr2(-1,0)*dp
c$$$     &     - 48.D0*Hr2(0,-1)*z2 + 48.D0*Hr2(0,-1)*z2*x + 1264.D0/9.D0*
c$$$     &    Hr2(0,0) + 196.D0/3.D0*Hr2(0,0)*x + 560.D0/9.D0*Hr2(0,0)*x**2
c$$$     &     + 2176.D0/9.D0*Hr2(0,0)*dx - 160.D0/9.D0*Hr2(0,0)*dp - 224.D0
c$$$     &    /9.D0*Hr2(0,0)*dm - 24.D0*Hr2(0,0)*z2 - 40.D0*Hr2(0,0)*z2*x
c$$$     &     + 44.D0/3.D0*Hr2(0,1) + 460.D0/3.D0*Hr2(0,1)*x - 176.D0/3.D0
c$$$     &    *Hr2(0,1)*x**2 + 176.D0/3.D0*Hr2(0,1)*dx + 320.D0/9.D0*Hr2(0,
c$$$     &    1)*dm - 48.D0*Hr2(0,1)*z2 - 48.D0*Hr2(0,1)*z2*x - 784.D0/9.D0
c$$$     &    *Hr2(1,0) + 464.D0/9.D0*Hr2(1,0)*x - 176.D0/3.D0*Hr2(1,0)*
c$$$     &    x**2 + 176.D0/3.D0*Hr2(1,0)*dx + 320.D0/9.D0*Hr2(1,0)*dm - 
c$$$     &    112.D0*Hr3(-1,-1,0) - 112.D0*Hr3(-1,-1,0)*x + 32.D0*Hr3(-1,-1
c$$$     &    ,0)*x**2 + 32.D0*Hr3(-1,-1,0)*dx - 40.D0/3.D0*Hr3(-1,0,0) + 
c$$$     &    88.D0/3.D0*Hr3(-1,0,0)*x - 160.D0/3.D0*Hr3(-1,0,0)*x**2 - 160.
c$$$     &    D0/3.D0*Hr3(-1,0,0)*dx )
c$$$      PTgg2 = PTgg2 + nf*ca**2 * ( 128.D0/3.D0*Hr3(-1,0,0)*dp + 32.D0*
c$$$     &    Hr3(-1,0,1) + 32.D0*Hr3(-1,0,1)*x + 32.D0/3.D0*Hr3(-1,0,1)*
c$$$     &    x**2 + 32.D0/3.D0*Hr3(-1,0,1)*dx + 104.D0*Hr3(0,-1,0) - 136.D0
c$$$     &    /3.D0*Hr3(0,-1,0)*x - 128.D0/3.D0*Hr3(0,-1,0)*x**2 - 32.D0*
c$$$     &    Hr3(0,-1,0)*dx + 64.D0/3.D0*Hr3(0,-1,0)*dp - 32.D0/3.D0*Hr3(0
c$$$     &    ,-1,0)*dm - 32.D0*Hr3(0,0,0) + 1000.D0/3.D0*Hr3(0,0,0)*x - 64.
c$$$     &    D0*Hr3(0,0,0)*x**2 + 128.D0*Hr3(0,0,0)*dx - 32.D0*Hr3(0,0,0)*
c$$$     &    dp + 96.D0*Hr3(0,0,0)*dm - 124.D0/3.D0*Hr3(0,0,1) + 44.D0/3.D0
c$$$     &    *Hr3(0,0,1)*x - 160.D0/3.D0*Hr3(0,0,1)*x**2 + 128.D0/3.D0*
c$$$     &    Hr3(0,0,1)*dx + 112.D0/3.D0*Hr3(0,0,1)*dm - 224.D0/3.D0*Hr3(0
c$$$     &    ,1,0) + 160.D0/3.D0*Hr3(0,1,0)*x - 128.D0/3.D0*Hr3(0,1,0)*
c$$$     &    x**2 + 128.D0/3.D0*Hr3(0,1,0)*dx + 128.D0/3.D0*Hr3(0,1,0)*dm
c$$$     &     - 140.D0/3.D0*Hr3(1,0,0) + 28.D0/3.D0*Hr3(1,0,0)*x - 128.D0/
c$$$     &    3.D0*Hr3(1,0,0)*x**2 + 128.D0/3.D0*Hr3(1,0,0)*dx + 112.D0/3.D0
c$$$     &    *Hr3(1,0,0)*dm - 96.D0*Hr4(0,-1,-1,0) + 96.D0*Hr4(0,-1,-1,0)*
c$$$     &    x )
c$$$      PTgg2 = PTgg2 + nf*ca**2 * ( 48.D0*Hr4(0,-1,0,0) - 48.D0*Hr4(0,-1
c$$$     &    ,0,0)*x + 80.D0*Hr4(0,0,-1,0) - 16.D0*Hr4(0,0,-1,0)*x + 16.D0
c$$$     &    *Hr4(0,0,0,0)*x + 24.D0*Hr4(0,0,0,1) + 24.D0*Hr4(0,0,0,1)*x
c$$$     &     + 24.D0*Hr4(0,1,0,0) + 24.D0*Hr4(0,1,0,0)*x )
c$$$      PTgg2 = PTgg2 + nf*cf*ca * ( 16180.D0/27.D0 - 12022.D0/27.D0*x + 
c$$$     &    50210.D0/81.D0*x**2 - 59714.D0/81.D0*dx - 110.D0/3.D0*dm - 
c$$$     &    776.D0/3.D0*z3 + 352.D0/3.D0*z3*x - 640.D0/3.D0*z3*x**2 + 512.
c$$$     &    D0/3.D0*z3*dx + 32.D0*z3*dm - 1432.D0/9.D0*z2 - 7972.D0/9.D0*
c$$$     &    z2*x - 352.D0/9.D0*z2*x**2 + 608.D0/9.D0*z2*dx + 616.D0/5.D0*
c$$$     &    z2**2 + 936.D0/5.D0*z2**2*x + 208.D0*Hr1(-1)*z2 + 208.D0*Hr1(
c$$$     &    -1)*z2*x - 64.D0/3.D0*Hr1(-1)*z2*x**2 - 64.D0/3.D0*Hr1(-1)*z2
c$$$     &    *dx - 5636.D0/27.D0*Hr1(0) - 25544.D0/27.D0*Hr1(0)*x - 500.D0
c$$$     &    *Hr1(0)*x**2 - 14552.D0/27.D0*Hr1(0)*dx - 16.D0*Hr1(0)*dm + 
c$$$     &    160.D0*Hr1(0)*z3 + 512.D0*Hr1(0)*z3*x - 304.D0/3.D0*Hr1(0)*z2
c$$$     &     + 8.D0/3.D0*Hr1(0)*z2*x - 32.D0/3.D0*Hr1(0)*z2*x**2 + 128.D0/
c$$$     &    3.D0*Hr1(0)*z2*dx + 2080.D0/9.D0*Hr1(1) - 76.D0/9.D0*Hr1(1)*x
c$$$     &     - 8068.D0/27.D0*Hr1(1)*x**2 + 2056.D0/27.D0*Hr1(1)*dx + 224.D
c$$$     &    0*Hr1(1)*z2 - 224.D0*Hr1(1)*z2*x - 448.D0/3.D0*Hr1(1)*z2*x**2
c$$$     &     + 448.D0/3.D0*Hr1(1)*z2*dx - 1664.D0/3.D0*Hr2(-1,0) - 1376.D0
c$$$     &    /3.D0*Hr2(-1,0)*x )
c$$$      PTgg2 = PTgg2 + nf*cf*ca * ( 2528.D0/9.D0*Hr2(-1,0)*x**2 + 1664.D0
c$$$     &    /9.D0*Hr2(-1,0)*dx + 128.D0*Hr2(0,-1)*z2 - 128.D0*Hr2(0,-1)*
c$$$     &    z2*x + 904.D0/3.D0*Hr2(0,0) + 2000.D0/3.D0*Hr2(0,0)*x + 1544.D
c$$$     &    0/9.D0*Hr2(0,0)*x**2 - 3712.D0/9.D0*Hr2(0,0)*dx - 16.D0*Hr2(0
c$$$     &    ,0)*z2 + 112.D0*Hr2(0,0)*z2*x + 1672.D0/9.D0*Hr2(0,1) + 2740.D
c$$$     &    0/9.D0*Hr2(0,1)*x + 272.D0*Hr2(0,1)*x**2 - 176.D0/9.D0*Hr2(0,
c$$$     &    1)*dx + 288.D0*Hr2(0,1)*z2 + 288.D0*Hr2(0,1)*z2*x + 148.D0/3.D
c$$$     &    0*Hr2(1,0) - 292.D0/3.D0*Hr2(1,0)*x + 1000.D0/9.D0*Hr2(1,0)*
c$$$     &    x**2 - 568.D0/9.D0*Hr2(1,0)*dx + 12.D0*Hr2(1,1) - 12.D0*Hr2(1
c$$$     &    ,1)*x - 16.D0/9.D0*Hr2(1,1)*x**2 + 16.D0/9.D0*Hr2(1,1)*dx + 
c$$$     &    288.D0*Hr3(-1,-1,0) + 288.D0*Hr3(-1,-1,0)*x - 256.D0/3.D0*
c$$$     &    Hr3(-1,-1,0)*x**2 - 256.D0/3.D0*Hr3(-1,-1,0)*dx - 208.D0*Hr3(
c$$$     &    -1,0,0) - 208.D0*Hr3(-1,0,0)*x + 224.D0/3.D0*Hr3(-1,0,0)*x**2
c$$$     &     + 224.D0/3.D0*Hr3(-1,0,0)*dx - 64.D0*Hr3(-1,0,1) - 64.D0*
c$$$     &    Hr3(-1,0,1)*x - 64.D0/3.D0*Hr3(-1,0,1)*x**2 - 64.D0/3.D0*Hr3(
c$$$     &    -1,0,1)*dx )
c$$$      PTgg2 = PTgg2 + nf*cf*ca * (  - 416.D0*Hr3(0,-1,0) + 160.D0*Hr3(0
c$$$     &    ,-1,0)*x + 256.D0/3.D0*Hr3(0,-1,0)*dx + 640.D0/3.D0*Hr3(0,0,0
c$$$     &    ) - 920.D0/3.D0*Hr3(0,0,0)*x - 640.D0/3.D0*Hr3(0,0,0)*dx - 
c$$$     &    128.D0*Hr3(0,0,1) - 88.D0*Hr3(0,0,1)*x - 32.D0*Hr3(0,0,1)*
c$$$     &    x**2 - 64.D0/3.D0*Hr3(0,0,1)*dx + 64.D0/3.D0*Hr3(0,1,0) - 104.
c$$$     &    D0/3.D0*Hr3(0,1,0)*x - 160.D0/3.D0*Hr3(0,1,0)*x**2 + 136.D0/3.
c$$$     &    D0*Hr3(0,1,1) + 88.D0/3.D0*Hr3(0,1,1)*x - 32.D0/3.D0*Hr3(0,1,
c$$$     &    1)*x**2 + 32.D0/3.D0*Hr3(0,1,1)*dx - 200.D0*Hr3(1,0,0) + 200.D
c$$$     &    0*Hr3(1,0,0)*x + 512.D0/3.D0*Hr3(1,0,0)*x**2 - 512.D0/3.D0*
c$$$     &    Hr3(1,0,0)*dx - 48.D0*Hr3(1,0,1) + 48.D0*Hr3(1,0,1)*x + 64.D0
c$$$     &    *Hr3(1,0,1)*x**2 - 64.D0*Hr3(1,0,1)*dx - 32.D0*Hr3(1,1,0) + 
c$$$     &    32.D0*Hr3(1,1,0)*x + 128.D0/3.D0*Hr3(1,1,0)*x**2 - 128.D0/3.D0
c$$$     &    *Hr3(1,1,0)*dx - 16.D0*Hr3(1,1,1) + 16.D0*Hr3(1,1,1)*x + 64.D0
c$$$     &    /3.D0*Hr3(1,1,1)*x**2 - 64.D0/3.D0*Hr3(1,1,1)*dx + 256.D0*
c$$$     &    Hr4(0,-1,-1,0) - 256.D0*Hr4(0,-1,-1,0)*x - 192.D0*Hr4(0,-1,0,
c$$$     &    0) )
c$$$      PTgg2 = PTgg2 + nf*cf*ca * ( 192.D0*Hr4(0,-1,0,0)*x - 256.D0*Hr4(
c$$$     &    0,0,-1,0) + 128.D0*Hr4(0,0,-1,0)*x + 128.D0*Hr4(0,0,0,0) - 
c$$$     &    448.D0*Hr4(0,0,0,0)*x - 112.D0*Hr4(0,0,0,1) - 112.D0*Hr4(0,0,
c$$$     &    0,1)*x + 80.D0*Hr4(0,0,1,0) + 80.D0*Hr4(0,0,1,0)*x + 48.D0*
c$$$     &    Hr4(0,0,1,1) + 48.D0*Hr4(0,0,1,1)*x - 304.D0*Hr4(0,1,0,0) - 
c$$$     &    304.D0*Hr4(0,1,0,0)*x - 96.D0*Hr4(0,1,0,1) - 96.D0*Hr4(0,1,0,
c$$$     &    1)*x - 64.D0*Hr4(0,1,1,0) - 64.D0*Hr4(0,1,1,0)*x - 32.D0*Hr4(
c$$$     &    0,1,1,1) - 32.D0*Hr4(0,1,1,1)*x )
c$$$      PTgg2 = PTgg2 + nf*cf**2 * (  - 1682.D0/9.D0 - 382.D0/9.D0*x + 
c$$$     &    292.D0/3.D0*x**2 + 132.D0*dx - 288.D0*z3 - 224.D0*z3*x - 256.D
c$$$     &    0/3.D0*z3*x**2 - 320.D0/3.D0*z3*dx + 824.D0/3.D0*z2 + 2516.D0/
c$$$     &    3.D0*z2*x + 320.D0*z2*x**2 - 152.D0*z2**2 - 184.D0*z2**2*x - 
c$$$     &    64.D0*Hr1(-1)*z2 - 64.D0*Hr1(-1)*z2*x + 64.D0/3.D0*Hr1(-1)*z2
c$$$     &    *x**2 + 64.D0/3.D0*Hr1(-1)*z2*dx + 370.D0/3.D0*Hr1(0) + 1786.D
c$$$     &    0/3.D0*Hr1(0)*x + 3604.D0/9.D0*Hr1(0)*x**2 + 32.D0*Hr1(0)*dx
c$$$     &     - 160.D0*Hr1(0)*z3 - 288.D0*Hr1(0)*z3*x - 80.D0*Hr1(0)*z2*x
c$$$     &     - 128.D0/3.D0*Hr1(0)*z2*x**2 - 1432.D0/3.D0*Hr1(1) + 716.D0/
c$$$     &    3.D0*Hr1(1)*x + 3316.D0/9.D0*Hr1(1)*x**2 - 1168.D0/9.D0*Hr1(1
c$$$     &    )*dx - 144.D0*Hr1(1)*z2 + 144.D0*Hr1(1)*z2*x + 128.D0*Hr1(1)*
c$$$     &    z2*x**2 - 128.D0*Hr1(1)*z2*dx + 1072.D0/3.D0*Hr2(-1,0) + 1072.
c$$$     &    D0/3.D0*Hr2(-1,0)*x - 64.D0*Hr2(0,-1)*z2 + 64.D0*Hr2(0,-1)*z2
c$$$     &    *x - 760.D0/3.D0*Hr2(0,0) - 784.D0*Hr2(0,0)*x - 160.D0/3.D0*
c$$$     &    Hr2(0,0)*x**2 + 16.D0*Hr2(0,0)*z2 + 16.D0*Hr2(0,0)*z2*x - 296.
c$$$     &    D0*Hr2(0,1) )
c$$$      PTgg2 = PTgg2 + nf*cf**2 * (  - 556.D0*Hr2(0,1)*x - 512.D0/3.D0*
c$$$     &    Hr2(0,1)*x**2 - 160.D0/3.D0*Hr2(0,1)*dx - 224.D0*Hr2(0,1)*z2
c$$$     &     - 224.D0*Hr2(0,1)*z2*x - 116.D0/3.D0*Hr2(1,0) - 28.D0/3.D0*
c$$$     &    Hr2(1,0)*x + 416.D0/3.D0*Hr2(1,0)*x**2 - 272.D0/3.D0*Hr2(1,0)
c$$$     &    *dx - 28.D0/3.D0*Hr2(1,1) + 28.D0/3.D0*Hr2(1,1)*x + 64.D0/3.D0
c$$$     &    *Hr2(1,1)*x**2 - 64.D0/3.D0*Hr2(1,1)*dx - 128.D0*Hr3(-1,-1,0)
c$$$     &     - 128.D0*Hr3(-1,-1,0)*x + 128.D0/3.D0*Hr3(-1,-1,0)*x**2 + 
c$$$     &    128.D0/3.D0*Hr3(-1,-1,0)*dx + 64.D0*Hr3(-1,0,0) + 64.D0*Hr3(
c$$$     &    -1,0,0)*x - 64.D0/3.D0*Hr3(-1,0,0)*x**2 - 64.D0/3.D0*Hr3(-1,0
c$$$     &    ,0)*dx + 128.D0*Hr3(0,-1,0) - 128.D0/3.D0*Hr3(0,-1,0)*x**2 - 
c$$$     &    72.D0*Hr3(0,0,0) - 40.D0*Hr3(0,0,0)*x - 160.D0/3.D0*Hr3(0,0,0
c$$$     &    )*x**2 + 16.D0*Hr3(0,0,1)*x - 128.D0/3.D0*Hr3(0,0,1)*dx - 192.
c$$$     &    D0*Hr3(0,1,0) - 176.D0*Hr3(0,1,0)*x - 64.D0/3.D0*Hr3(0,1,0)*
c$$$     &    x**2 - 256.D0/3.D0*Hr3(0,1,0)*dx - 64.D0*Hr3(0,1,1) - 48.D0*
c$$$     &    Hr3(0,1,1)*x - 64.D0/3.D0*Hr3(0,1,1)*dx + 152.D0*Hr3(1,0,0)
c$$$     &     - 152.D0*Hr3(1,0,0)*x )
c$$$      PTgg2 = PTgg2 + nf*cf**2 * (  - 416.D0/3.D0*Hr3(1,0,0)*x**2 + 416.
c$$$     &    D0/3.D0*Hr3(1,0,0)*dx + 48.D0*Hr3(1,0,1) - 48.D0*Hr3(1,0,1)*x
c$$$     &     - 64.D0*Hr3(1,0,1)*x**2 + 64.D0*Hr3(1,0,1)*dx + 32.D0*Hr3(1,
c$$$     &    1,0) - 32.D0*Hr3(1,1,0)*x - 128.D0/3.D0*Hr3(1,1,0)*x**2 + 128.
c$$$     &    D0/3.D0*Hr3(1,1,0)*dx + 16.D0*Hr3(1,1,1) - 16.D0*Hr3(1,1,1)*x
c$$$     &     - 64.D0/3.D0*Hr3(1,1,1)*x**2 + 64.D0/3.D0*Hr3(1,1,1)*dx - 
c$$$     &    128.D0*Hr4(0,-1,-1,0) + 128.D0*Hr4(0,-1,-1,0)*x + 64.D0*Hr4(0
c$$$     &    ,-1,0,0) - 64.D0*Hr4(0,-1,0,0)*x + 128.D0*Hr4(0,0,-1,0) - 32.D
c$$$     &    0*Hr4(0,0,0,0) - 32.D0*Hr4(0,0,0,0)*x - 16.D0*Hr4(0,0,0,1) - 
c$$$     &    16.D0*Hr4(0,0,0,1)*x - 144.D0*Hr4(0,0,1,0) - 144.D0*Hr4(0,0,1
c$$$     &    ,0)*x - 48.D0*Hr4(0,0,1,1) - 48.D0*Hr4(0,0,1,1)*x + 240.D0*
c$$$     &    Hr4(0,1,0,0) + 240.D0*Hr4(0,1,0,0)*x + 96.D0*Hr4(0,1,0,1) + 
c$$$     &    96.D0*Hr4(0,1,0,1)*x + 64.D0*Hr4(0,1,1,0) + 64.D0*Hr4(0,1,1,0
c$$$     &    )*x + 32.D0*Hr4(0,1,1,1) + 32.D0*Hr4(0,1,1,1)*x )
c$$$      PTgg2 = PTgg2 + nf2*ca * (  - 110.D0/27.D0 + 14.D0/3.D0*x - 472.
c$$$     &    D0/81.D0*x**2 + 472.D0/81.D0*dx - 16.D0/27.D0*dm + 32.D0/9.D0
c$$$     &    *z2 + 32.D0/9.D0*z2*x - 104.D0/9.D0*Hr1(0) + 130.D0/9.D0*Hr1(
c$$$     &    0)*x - 88.D0/9.D0*Hr1(0)*x**2 + 368.D0/27.D0*Hr1(0)*dx + 160.D
c$$$     &    0/27.D0*Hr1(0)*dm + 22.D0/9.D0*Hr1(1) - 22.D0/9.D0*Hr1(1)*x
c$$$     &     + 104.D0/27.D0*Hr1(1)*x**2 - 104.D0/27.D0*Hr1(1)*dx - 80.D0/
c$$$     &    9.D0*Hr2(0,0) + 112.D0/9.D0*Hr2(0,0)*x - 64.D0/9.D0*Hr2(0,0)*
c$$$     &    x**2 + 64.D0/9.D0*Hr2(0,0)*dx + 64.D0/9.D0*Hr2(0,0)*dm - 32.D0
c$$$     &    /9.D0*Hr2(0,1) - 32.D0/9.D0*Hr2(0,1)*x )
c$$$      PTgg2 = PTgg2 + nf2*cf * ( 544.D0/9.D0 - 304.D0/9.D0*x - 1216.D0
c$$$     &    /81.D0*x**2 - 944.D0/81.D0*dx - 64.D0/3.D0*z3 - 64.D0/3.D0*z3
c$$$     &    *x - 200.D0/9.D0*z2 - 440.D0/9.D0*z2*x - 224.D0/9.D0*z2*x**2
c$$$     &     - 208.D0/27.D0*Hr1(0) - 64.D0/27.D0*Hr1(0)*x + 16.D0*Hr1(0)*
c$$$     &    x**2 - 736.D0/27.D0*Hr1(0)*dx + 16.D0/3.D0*Hr1(0)*z2 + 16.D0/
c$$$     &    3.D0*Hr1(0)*z2*x - 32.D0/3.D0*Hr1(1) + 80.D0/3.D0*Hr1(1)*x - 
c$$$     &    1168.D0/27.D0*Hr1(1)*x**2 + 736.D0/27.D0*Hr1(1)*dx - 160.D0/3.
c$$$     &    D0*Hr2(0,0) - 64.D0*Hr2(0,0)*x - 160.D0/9.D0*Hr2(0,0)*x**2 - 
c$$$     &    128.D0/9.D0*Hr2(0,0)*dx + 296.D0/9.D0*Hr2(0,1) + 344.D0/9.D0*
c$$$     &    Hr2(0,1)*x + 32.D0/3.D0*Hr2(0,1)*x**2 + 128.D0/9.D0*Hr2(0,1)*
c$$$     &    dx - 8.D0/3.D0*Hr2(1,0) + 8.D0/3.D0*Hr2(1,0)*x + 32.D0/9.D0*
c$$$     &    Hr2(1,0)*x**2 - 32.D0/9.D0*Hr2(1,0)*dx - 8.D0/3.D0*Hr2(1,1)
c$$$     &     + 8.D0/3.D0*Hr2(1,1)*x + 32.D0/9.D0*Hr2(1,1)*x**2 - 32.D0/9.D
c$$$     &    0*Hr2(1,1)*dx - 64.D0/3.D0*Hr3(0,0,0) - 64.D0/3.D0*Hr3(0,0,0)
c$$$     &    *x + 16.D0*Hr3(0,0,1) + 16.D0*Hr3(0,0,1)*x - 16.D0/3.D0*Hr3(0
c$$$     &    ,1,0) )
c$$$      PTgg2 = PTgg2 + nf2*cf * (  - 16.D0/3.D0*Hr3(0,1,0)*x - 16.D0/3.
c$$$     &    D0*Hr3(0,1,1) - 16.D0/3.D0*Hr3(0,1,1)*x )
c$$$*
c$$$* ...The soft (`+'-distribution) part of the splitting function
c$$$*
c$$$       A3G =
c$$$     ,      ca**3    * ( + 490.D0/3.D0 + 88.D0/3.D0*z3 - 1072.D0/9.D0*z2
c$$$     ,                   + 176.D0/5.D0*z2**2 )
c$$$     ,    + ca**2*nf * ( - 836./27.D0 + 160./9.D0*z2 - 112./3.D0*z3 )
c$$$     ,    + ca*cf*nf * ( - 110./3.D0 + 32.*z3 ) - ca*nf2 * 16./27.D0
c$$$*
c$$$       GGG2L = DM * A3G
c$$$*
c$$$* ...The regular piece of the splitting function
c$$$*
c$$$       X2GGTA = PTgg2 - GGG2L
c$$$*
c$$$       RETURN
c$$$       END
c$$$*
c$$$* ---------------------------------------------------------------------
c$$$*
c$$$*
c$$$* ..The singular (soft) piece of P_gg^(2)T  (as in the spacelike case)
c$$$*
c$$$       FUNCTION X2GGTB (Y, NF)
c$$$       IMPLICIT REAL*8 (A - Z)
c$$$       INTEGER NF
c$$$*
c$$$       COMMON / P2GSOFT / A3G
c$$$*
c$$$       X2GGTB  = A3G/(1.D0-Y)
c$$$*
c$$$       RETURN
c$$$       END
c$$$*
c$$$* ---------------------------------------------------------------------
c$$$*
c$$$*
c$$$* ..The 'local' piece of P_gg^(2)T  (as in the spacelike case)
c$$$*
c$$$       FUNCTION X2GGTC (Y, NF)
c$$$*
c$$$       IMPLICIT REAL*8 (A - Z)
c$$$       INTEGER NF, NF2
c$$$       PARAMETER ( Z2 = 1.6449 34066 84822 64365 D0,
c$$$     ,             Z3 = 1.2020 56903 15959 42854 D0,
c$$$     ,             Z5 = 1.0369 27755 14336 99263 D0 )
c$$$*
c$$$       COMMON / P2GSOFT / A3G
c$$$*
c$$$* ...Colour factors
c$$$*
c$$$       CF  = 4./3.D0
c$$$       CA  = 3.D0
c$$$       NF2 = NF*NF
c$$$*
c$$$* ...The coefficient of delta(1-x)
c$$$*
c$$$       P2GDELT = 
c$$$     ,    + 79.D0/2.D0*ca**3
c$$$     ,    + 8.D0/3.D0*z2*ca**3
c$$$     ,    + 22.D0/3.D0*z2**2*ca**3
c$$$     ,    + 536.D0/3.D0*z3*ca**3
c$$$     ,    - 16.D0*z2*z3*ca**3
c$$$     ,    - 80.D0*z5*ca**3
c$$$     ,    + cf**2*nf
c$$$     ,    - 233.D0/18*ca**2*nf
c$$$     ,    - 8.D0/3.D0*z2*ca**2*nf
c$$$     ,    - 80.D0/3.D0*z3*ca**2*nf
c$$$     ,    - 4.D0/3.D0*z2**2*ca**2*nf
c$$$     ,    - 241.D0/18.D0*ca*cf*nf
c$$$     ,    + 29.D0/18.D0*ca*nf2
c$$$     ,    + 11.D0/9.D0*cf*nf2
c$$$*
c$$$       X2GGTC = LOG (1.D0-Y) * A3G + P2GDELT
c$$$*
c$$$       RETURN
c$$$       END
c$$$*
c$$$* =================================================================av==
c$$$*
c$$$* ..File: xpns2e.f
c$$$*
c$$$*
c$$$* ..The exact 3-loop MS(bar) non-singlet splitting functions P_NS^(2)T
c$$$*    for the evolution of unpolarized fragmentation densities at
c$$$*    mu_r = mu_f.  The expansion parameter is alpha_s/(4 pi).
c$$$*
c$$$* ..The distributions (in the mathematical sense) are given as in eq.
c$$$*    (B.26) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.
c$$$*    The name-endings A, B, and C of the functions below correspond to
c$$$*    the kernel superscripts [2], [3], and [1] in that equation.
c$$$*
c$$$* ..The code uses the package of Gehrmann and Remiddi for the harmonic
c$$$*    polylogarithms published in hep-ph/0107173 = CPC 141 (2001) 296.
c$$$*
c$$$* ..References: A. Mitov, S. Moch and A. Vogt,
c$$$*               Phys. Lett. B638 (2006) 61, hep-ph/0604053   (+ and -)
c$$$*               S. Moch, J. Vermaseren and A. Vogt,
c$$$*               Nucl. Phys. B688 (2004) 101, hep-ph/0403192        (S)
c$$$*
c$$$* =====================================================================
c$$$*
c$$$*
c$$$* ..This is the regular piece of P_NS+  (the even-moment combination)
c$$$*
c$$$       FUNCTION X2NSPTA (X, NF)
c$$$*
c$$$       IMPLICIT REAL*8 (A - Z)
c$$$       COMPLEX*16 HC1, HC2, HC3, HC4
c$$$       INTEGER NF, NF2, N1, N2, NW, I1, I2, I3, N
c$$$       PARAMETER ( N1 = -1, N2 = 1, NW = 4 )
c$$$       DIMENSION HC1(N1:N2),HC2(N1:N2,N1:N2),HC3(N1:N2,N1:N2,N1:N2),
c$$$     ,           HC4(N1:N2,N1:N2,N1:N2,N1:N2)
c$$$       DIMENSION HR1(N1:N2),HR2(N1:N2,N1:N2),HR3(N1:N2,N1:N2,N1:N2),
c$$$     ,           HR4(N1:N2,N1:N2,N1:N2,N1:N2)
c$$$       DIMENSION HI1(N1:N2),HI2(N1:N2,N1:N2),HI3(N1:N2,N1:N2,N1:N2),
c$$$     ,           HI4(N1:N2,N1:N2,N1:N2,N1:N2)
c$$$       PARAMETER ( Z2 = 1.6449 34066 84822 64365 D0,
c$$$     ,             Z3 = 1.2020 56903 15959 42854 D0,
c$$$     ,             Z5 = 1.0369 27755 14336 99263 D0 )
c$$$*
c$$$* ..The soft coefficient for use in X2NSPB and X2NSPC
c$$$*
c$$$       COMMON / P2SOFT / A3
c$$$*
c$$$* ...Colour factors
c$$$*
c$$$       CF  = 4./3.D0
c$$$       CA  = 3.D0
c$$$       NF2 = NF*NF
c$$$*
c$$$* ...Some abbreviations
c$$$*
c$$$       DX = 1.D0/X
c$$$       DM = 1.D0/(1.D0-X)
c$$$       DP = 1.D0/(1.D0+X)
c$$$*
c$$$* ...The harmonic polylogs up to weight 4 by Gehrmann and Remiddi
c$$$*
c$$$       CALL HPLOG (X, NW, HC1,HC2,HC3,HC4, HR1,HR2,HR3,HR4,
c$$$     ,            HI1,HI2,HI3,HI4, N1, N2)
c$$$*
c$$$* ...The splitting function in terms of the harmonic polylogs
c$$$*    (without the delta(1-x) part, but with the soft contribution)
c$$$*
c$$$      PTqq2 =
c$$$     &  + cf*ca**2 * ( 5327.D0/27.D0 - 9737.D0/27.D0*x + 490.D0/3.D0*dm
c$$$     &     - 224.D0*z3*x - 88.D0*z3*dp + 16.D0*z3*dm - 112.D0*z2 + 448.D
c$$$     &    0/9.D0*z2*x + 1072.D0/9.D0*z2*dp - 1072.D0/9.D0*z2*dm - 62.D0/
c$$$     &    5.D0*z2**2 - 242.D0/5.D0*z2**2*x - 32.D0*z2**2*dp + 384.D0/5.D
c$$$     &    0*z2**2*dm - 192.D0*Hr1(-1)*z3 + 192.D0*Hr1(-1)*z3*x + 384.D0
c$$$     &    *Hr1(-1)*z3*dp + 208.D0/3.D0*Hr1(-1)*z2 + 560.D0/3.D0*Hr1(-1)
c$$$     &    *z2*x + 352.D0/3.D0*Hr1(-1)*z2*dp + 410.D0/27.D0*Hr1(0) - 
c$$$     &    8686.D0/27.D0*Hr1(0)*x - 24.D0*Hr1(0)*dp + 4172.D0/27.D0*Hr1(
c$$$     &    0)*dm - 144.D0*Hr1(0)*z3*x - 128.D0*Hr1(0)*z3*dp + 128.D0*
c$$$     &    Hr1(0)*z3*dm - 4.D0*Hr1(0)*z2 - 148.D0/3.D0*Hr1(0)*z2*x - 16.D
c$$$     &    0/3.D0*Hr1(0)*z2*dp - 248.D0/3.D0*Hr1(0)*z2*dm + 176.D0*Hr1(1
c$$$     &    ) - 176.D0*Hr1(1)*x - 144.D0*Hr1(1)*z3 - 144.D0*Hr1(1)*z3*x
c$$$     &     + 288.D0*Hr1(1)*z3*dm + 32.D0*Hr1(1)*z2 - 32.D0*Hr1(1)*z2*x
c$$$     &     + 256.D0*Hr2(-1,-1)*z2 - 256.D0*Hr2(-1,-1)*z2*x - 512.D0*
c$$$     &    Hr2(-1,-1)*z2*dp - 688.D0/9.D0*Hr2(-1,0) + 1456.D0/9.D0*Hr2(
c$$$     &    -1,0)*x )
c$$$      PTqq2 = PTqq2 + cf*ca**2 * ( 2144.D0/9.D0*Hr2(-1,0)*dp - 176.D0*
c$$$     &    Hr2(-1,0)*z2 + 176.D0*Hr2(-1,0)*z2*x + 352.D0*Hr2(-1,0)*z2*dp
c$$$     &     - 136.D0*Hr2(0,-1)*z2 + 136.D0*Hr2(0,-1)*z2*x + 256.D0*Hr2(0
c$$$     &    ,-1)*z2*dp - 242.D0/9.D0*Hr2(0,0) - 230.D0/3.D0*Hr2(0,0)*x - 
c$$$     &    1072.D0/9.D0*Hr2(0,0)*dp + 1556.D0/9.D0*Hr2(0,0)*dm + 36.D0*
c$$$     &    Hr2(0,0)*z2 - 68.D0*Hr2(0,0)*z2*x - 96.D0*Hr2(0,0)*z2*dp + 
c$$$     &    112.D0*Hr2(0,1) + 112.D0*Hr2(0,1)*x - 40.D0*Hr2(0,1)*z2 + 24.D
c$$$     &    0*Hr2(0,1)*z2*x + 64.D0*Hr2(0,1)*z2*dp + 16.D0*Hr2(1,0)*z2 + 
c$$$     &    16.D0*Hr2(1,0)*z2*x - 32.D0*Hr2(1,0)*z2*dm + 64.D0*Hr3(-1,-1,
c$$$     &    0) + 64.D0*Hr3(-1,-1,0)*x - 328.D0/3.D0*Hr3(-1,0,0) - 152.D0/
c$$$     &    3.D0*Hr3(-1,0,0)*x + 176.D0/3.D0*Hr3(-1,0,0)*dp - 112.D0/3.D0
c$$$     &    *Hr3(-1,0,1) - 464.D0/3.D0*Hr3(-1,0,1)*x - 352.D0/3.D0*Hr3(-1
c$$$     &    ,0,1)*dp - 88.D0/3.D0*Hr3(0,-1,0) + 40.D0/3.D0*Hr3(0,-1,0)*x
c$$$     &     + 176.D0/3.D0*Hr3(0,-1,0)*dp - 48.D0*Hr3(0,-1,0)*dm - 128.D0/
c$$$     &    3.D0*Hr3(0,0,0)*x - 248.D0/3.D0*Hr3(0,0,0)*dp + 248.D0/3.D0*
c$$$     &    Hr3(0,0,0)*dm )
c$$$      PTqq2 = PTqq2 + cf*ca**2 * ( 4.D0*Hr3(0,0,1) + 188.D0/3.D0*Hr3(0,
c$$$     &    0,1)*x + 176.D0/3.D0*Hr3(0,0,1)*dp + 88.D0/3.D0*Hr3(0,0,1)*dm
c$$$     &     + 12.D0*Hr3(1,0,0) + 76.D0*Hr3(1,0,0)*x - 88.D0*Hr3(1,0,0)*
c$$$     &    dm - 128.D0*Hr4(-1,-1,0,0) + 128.D0*Hr4(-1,-1,0,0)*x + 256.D0
c$$$     &    *Hr4(-1,-1,0,0)*dp - 256.D0*Hr4(-1,-1,0,1) + 256.D0*Hr4(-1,-1
c$$$     &    ,0,1)*x + 512.D0*Hr4(-1,-1,0,1)*dp + 48.D0*Hr4(-1,0,0,0) - 48.
c$$$     &    D0*Hr4(-1,0,0,0)*x - 96.D0*Hr4(-1,0,0,0)*dp + 128.D0*Hr4(-1,0
c$$$     &    ,0,1) - 128.D0*Hr4(-1,0,0,1)*x - 256.D0*Hr4(-1,0,0,1)*dp - 80.
c$$$     &    D0*Hr4(0,-1,-1,0) - 48.D0*Hr4(0,-1,-1,0)*x + 128.D0*Hr4(0,-1,
c$$$     &    -1,0)*dm + 88.D0*Hr4(0,-1,0,0) - 56.D0*Hr4(0,-1,0,0)*x - 128.D
c$$$     &    0*Hr4(0,-1,0,0)*dp - 32.D0*Hr4(0,-1,0,0)*dm + 96.D0*Hr4(0,-1,
c$$$     &    0,1) - 160.D0*Hr4(0,-1,0,1)*x - 256.D0*Hr4(0,-1,0,1)*dp + 64.D
c$$$     &    0*Hr4(0,-1,0,1)*dm + 8.D0*Hr4(0,0,-1,0) - 40.D0*Hr4(0,0,-1,0)
c$$$     &    *x - 32.D0*Hr4(0,0,-1,0)*dp + 32.D0*Hr4(0,0,-1,0)*dm + 40.D0*
c$$$     &    Hr4(0,0,0,0)*x + 32.D0*Hr4(0,0,0,0)*dp - 32.D0*Hr4(0,0,0,0)*
c$$$     &    dm )
c$$$      PTqq2 = PTqq2 + cf*ca**2 * (  - 36.D0*Hr4(0,0,0,1) + 28.D0*Hr4(0,
c$$$     &    0,0,1)*x + 64.D0*Hr4(0,0,0,1)*dp + 32.D0*Hr4(0,0,0,1)*dm + 28.
c$$$     &    D0*Hr4(0,1,0,0) + 28.D0*Hr4(0,1,0,0)*x - 64.D0*Hr4(0,1,0,0)*
c$$$     &    dm - 96.D0*Hr4(1,0,-1,0) - 96.D0*Hr4(1,0,-1,0)*x + 192.D0*
c$$$     &    Hr4(1,0,-1,0)*dm + 48.D0*Hr4(1,0,0,0) + 48.D0*Hr4(1,0,0,0)*x
c$$$     &     - 96.D0*Hr4(1,0,0,0)*dm - 64.D0*Hr4(1,0,0,1) - 64.D0*Hr4(1,0
c$$$     &    ,0,1)*x + 128.D0*Hr4(1,0,0,1)*dm + 64.D0*Hr4(1,1,0,0) + 64.D0
c$$$     &    *Hr4(1,1,0,0)*x - 128.D0*Hr4(1,1,0,0)*dm )
c$$$      PTqq2 = PTqq2 + cf**2*ca * ( 532.D0/9.D0 - 532.D0/9.D0*x - 16.D0/
c$$$     &    3.D0*z3 + 2336.D0/3.D0*z3*x + 248.D0*z3*dp + 80.D0/3.D0*z3*dm
c$$$     &     + 3448.D0/9.D0*z2 + 2024.D0/9.D0*z2*x - 2144.D0/9.D0*z2*dp
c$$$     &     - 24.D0/5.D0*z2**2 + 56.D0/5.D0*z2**2*x + 8.D0*z2**2*dp - 
c$$$     &    552.D0/5.D0*z2**2*dm + 672.D0*Hr1(-1)*z3 - 672.D0*Hr1(-1)*z3*
c$$$     &    x - 1344.D0*Hr1(-1)*z3*dp - 992.D0/3.D0*Hr1(-1)*z2 - 1984.D0/
c$$$     &    3.D0*Hr1(-1)*z2*x - 992.D0/3.D0*Hr1(-1)*z2*dp - 562.D0/3.D0*
c$$$     &    Hr1(0) + 5650.D0/9.D0*Hr1(0)*x + 72.D0*Hr1(0)*dp + 302.D0/3.D0
c$$$     &    *Hr1(0)*dm + 224.D0*Hr1(0)*z3*x + 240.D0*Hr1(0)*z3*dp - 240.D0
c$$$     &    *Hr1(0)*z3*dm + 40.D0/3.D0*Hr1(0)*z2 + 328.D0*Hr1(0)*z2*x + 
c$$$     &    248.D0/3.D0*Hr1(0)*z2*dp + 536.D0/3.D0*Hr1(0)*z2*dm - 1672.D0/
c$$$     &    3.D0*Hr1(1) + 1672.D0/3.D0*Hr1(1)*x + 384.D0*Hr1(1)*z3 + 384.D
c$$$     &    0*Hr1(1)*z3*x - 768.D0*Hr1(1)*z3*dm - 112.D0*Hr1(1)*z2 + 112.D
c$$$     &    0*Hr1(1)*z2*x - 896.D0*Hr2(-1,-1)*z2 + 896.D0*Hr2(-1,-1)*z2*x
c$$$     &     + 1792.D0*Hr2(-1,-1)*z2*dp + 1520.D0/9.D0*Hr2(-1,0) - 2768.D0
c$$$     &    /9.D0*Hr2(-1,0)*x )
c$$$      PTqq2 = PTqq2 + cf**2*ca * (  - 4288.D0/9.D0*Hr2(-1,0)*dp + 544.D0
c$$$     &    *Hr2(-1,0)*z2 - 544.D0*Hr2(-1,0)*z2*x - 1088.D0*Hr2(-1,0)*z2*
c$$$     &    dp + 448.D0*Hr2(0,-1)*z2 - 352.D0*Hr2(0,-1)*z2*x - 768.D0*
c$$$     &    Hr2(0,-1)*z2*dp - 96.D0*Hr2(0,-1)*z2*dm + 32.D0*Hr2(0,0) + 
c$$$     &    1736.D0/9.D0*Hr2(0,0)*x + 2144.D0/9.D0*Hr2(0,0)*dp - 3640.D0/
c$$$     &    9.D0*Hr2(0,0)*dm - 96.D0*Hr2(0,0)*z2 + 32.D0*Hr2(0,0)*z2*x + 
c$$$     &    160.D0*Hr2(0,0)*z2*dp + 128.D0*Hr2(0,0)*z2*dm - 2648.D0/9.D0*
c$$$     &    Hr2(0,1) - 1304.D0/9.D0*Hr2(0,1)*x - 2144.D0/9.D0*Hr2(0,1)*dm
c$$$     &     + 96.D0*Hr2(0,1)*z2 - 128.D0*Hr2(0,1)*z2*x - 224.D0*Hr2(0,1)
c$$$     &    *z2*dp + 64.D0*Hr2(0,1)*z2*dm + 400.D0/9.D0*Hr2(1,0) + 1744.D0
c$$$     &    /9.D0*Hr2(1,0)*x - 2144.D0/9.D0*Hr2(1,0)*dm - 32.D0*Hr2(1,0)*
c$$$     &    z2 - 32.D0*Hr2(1,0)*z2*x + 64.D0*Hr2(1,0)*z2*dm - 224.D0*Hr3(
c$$$     &    -1,-1,0) - 224.D0*Hr3(-1,-1,0)*x + 680.D0/3.D0*Hr3(-1,0,0) + 
c$$$     &    760.D0/3.D0*Hr3(-1,0,0)*x + 80.D0/3.D0*Hr3(-1,0,0)*dp + 656.D0
c$$$     &    /3.D0*Hr3(-1,0,1) + 1648.D0/3.D0*Hr3(-1,0,1)*x + 992.D0/3.D0*
c$$$     &    Hr3(-1,0,1)*dp )
c$$$      PTqq2 = PTqq2 + cf**2*ca * ( 32.D0/3.D0*Hr3(0,-1,0) + 160.D0/3.D0
c$$$     &    *Hr3(0,-1,0)*x - 208.D0/3.D0*Hr3(0,-1,0)*dp + 96.D0*Hr3(0,-1,
c$$$     &    0)*dm + 132.D0*Hr3(0,0,0) + 460.D0/3.D0*Hr3(0,0,0)*x + 280.D0/
c$$$     &    3.D0*Hr3(0,0,0)*dp - 808.D0/3.D0*Hr3(0,0,0)*dm + 136.D0/3.D0*
c$$$     &    Hr3(0,0,1) - 216.D0*Hr3(0,0,1)*x - 496.D0/3.D0*Hr3(0,0,1)*dp
c$$$     &     - 640.D0/3.D0*Hr3(0,0,1)*dm + 176.D0/3.D0*Hr3(0,1,0) + 176.D0
c$$$     &    /3.D0*Hr3(0,1,0)*x - 352.D0/3.D0*Hr3(0,1,0)*dm + 16.D0*Hr3(1,
c$$$     &    0,0) - 112.D0*Hr3(1,0,0)*x + 96.D0*Hr3(1,0,0)*dm + 320.D0*
c$$$     &    Hr4(-1,-1,0,0) - 320.D0*Hr4(-1,-1,0,0)*x - 640.D0*Hr4(-1,-1,0
c$$$     &    ,0)*dp + 896.D0*Hr4(-1,-1,0,1) - 896.D0*Hr4(-1,-1,0,1)*x - 
c$$$     &    1792.D0*Hr4(-1,-1,0,1)*dp - 64.D0*Hr4(-1,0,-1,0) + 64.D0*Hr4(
c$$$     &    -1,0,-1,0)*x + 128.D0*Hr4(-1,0,-1,0)*dp + 16.D0*Hr4(-1,0,0,0)
c$$$     &     - 16.D0*Hr4(-1,0,0,0)*x - 32.D0*Hr4(-1,0,0,0)*dp - 384.D0*
c$$$     &    Hr4(-1,0,0,1) + 384.D0*Hr4(-1,0,0,1)*x + 768.D0*Hr4(-1,0,0,1)
c$$$     &    *dp + 32.D0*Hr4(-1,0,1,0) - 32.D0*Hr4(-1,0,1,0)*x - 64.D0*
c$$$     &    Hr4(-1,0,1,0)*dp )
c$$$      PTqq2 = PTqq2 + cf**2*ca * ( 192.D0*Hr4(0,-1,-1,0) + 256.D0*Hr4(0
c$$$     &    ,-1,-1,0)*x + 128.D0*Hr4(0,-1,-1,0)*dp - 448.D0*Hr4(0,-1,-1,0
c$$$     &    )*dm - 176.D0*Hr4(0,-1,0,0) + 16.D0*Hr4(0,-1,0,0)*x + 224.D0*
c$$$     &    Hr4(0,-1,0,0)*dp + 160.D0*Hr4(0,-1,0,0)*dm - 352.D0*Hr4(0,-1,
c$$$     &    0,1) + 480.D0*Hr4(0,-1,0,1)*x + 832.D0*Hr4(0,-1,0,1)*dp - 128.
c$$$     &    D0*Hr4(0,-1,0,1)*dm + 64.D0*Hr4(0,0,-1,0) - 64.D0*Hr4(0,0,-1,
c$$$     &    0)*x - 96.D0*Hr4(0,0,-1,0)*dp - 32.D0*Hr4(0,0,-1,0)*dm + 160.D
c$$$     &    0*Hr4(0,0,0,0)*x + 96.D0*Hr4(0,0,0,0)*dp - 96.D0*Hr4(0,0,0,0)
c$$$     &    *dm + 96.D0*Hr4(0,0,0,1) - 32.D0*Hr4(0,0,0,1)*x - 128.D0*Hr4(
c$$$     &    0,0,0,1)*dp - 160.D0*Hr4(0,0,0,1)*dm + 32.D0*Hr4(0,0,1,0)*x
c$$$     &     + 32.D0*Hr4(0,0,1,0)*dp - 32.D0*Hr4(0,0,1,0)*dm - 16.D0*Hr4(
c$$$     &    0,1,0,0) - 16.D0*Hr4(0,1,0,0)*x + 96.D0*Hr4(0,1,0,0)*dm + 256.
c$$$     &    D0*Hr4(1,0,-1,0) + 256.D0*Hr4(1,0,-1,0)*x - 512.D0*Hr4(1,0,-1
c$$$     &    ,0)*dm - 80.D0*Hr4(1,0,0,0) - 80.D0*Hr4(1,0,0,0)*x + 160.D0*
c$$$     &    Hr4(1,0,0,0)*dm + 128.D0*Hr4(1,0,0,1) + 128.D0*Hr4(1,0,0,1)*x
c$$$     &     - 256.D0*Hr4(1,0,0,1)*dm )
c$$$      PTqq2 = PTqq2 + cf**2*ca * (  - 128.D0*Hr4(1,1,0,0) - 128.D0*Hr4(
c$$$     &    1,1,0,0)*x + 256.D0*Hr4(1,1,0,0)*dm )
c$$$      PTqq2 = PTqq2 + cf**3 * (  - 62.D0 + 62.D0*x - 48.D0*z3 - 720.D0*
c$$$     &    z3*x - 144.D0*z3*dp - 308.D0*z2 - 372.D0*z2*x - 56.D0/5.D0*
c$$$     &    z2**2 + 504.D0/5.D0*z2**2*x + 112.D0*z2**2*dp + 144.D0/5.D0*
c$$$     &    z2**2*dm - 576.D0*Hr1(-1)*z3 + 576.D0*Hr1(-1)*z3*x + 1152.D0*
c$$$     &    Hr1(-1)*z3*dp + 384.D0*Hr1(-1)*z2 + 576.D0*Hr1(-1)*z2*x + 192.
c$$$     &    D0*Hr1(-1)*z2*dp + 42.D0*Hr1(0) - 590.D0*Hr1(0)*x - 48.D0*
c$$$     &    Hr1(0)*dp + 6.D0*Hr1(0)*dm + 128.D0*Hr1(0)*z3*x + 32.D0*Hr1(0
c$$$     &    )*z3*dp - 32.D0*Hr1(0)*z3*dm + 48.D0*Hr1(0)*z2 - 400.D0*Hr1(0
c$$$     &    )*z2*x - 144.D0*Hr1(0)*z2*dp - 144.D0*Hr1(0)*z2*dm + 560.D0*
c$$$     &    Hr1(1) - 560.D0*Hr1(1)*x - 192.D0*Hr1(1)*z3 - 192.D0*Hr1(1)*
c$$$     &    z3*x + 384.D0*Hr1(1)*z3*dm + 96.D0*Hr1(1)*z2 - 96.D0*Hr1(1)*
c$$$     &    z2*x + 768.D0*Hr2(-1,-1)*z2 - 768.D0*Hr2(-1,-1)*z2*x - 1536.D0
c$$$     &    *Hr2(-1,-1)*z2*dp - 32.D0*Hr2(-1,0) - 32.D0*Hr2(-1,0)*x - 384.
c$$$     &    D0*Hr2(-1,0)*z2 + 384.D0*Hr2(-1,0)*z2*x + 768.D0*Hr2(-1,0)*z2
c$$$     &    *dp - 352.D0*Hr2(0,-1)*z2 + 160.D0*Hr2(0,-1)*z2*x + 512.D0*
c$$$     &    Hr2(0,-1)*z2*dp )
c$$$      PTqq2 = PTqq2 + cf**3 * ( 192.D0*Hr2(0,-1)*z2*dm - 20.D0*Hr2(0,0)
c$$$     &     + 316.D0*Hr2(0,0)*x + 52.D0*Hr2(0,0)*dm + 48.D0*Hr2(0,0)*z2
c$$$     &     + 208.D0*Hr2(0,0)*z2*x + 64.D0*Hr2(0,0)*z2*dp - 256.D0*Hr2(0
c$$$     &    ,0)*z2*dm + 340.D0*Hr2(0,1) + 308.D0*Hr2(0,1)*x - 96.D0*Hr2(0
c$$$     &    ,1)*z2 + 96.D0*Hr2(0,1)*z2*x + 192.D0*Hr2(0,1)*z2*dp + 16.D0*
c$$$     &    Hr2(1,0) - 16.D0*Hr2(1,0)*x + 192.D0*Hr3(-1,-1,0) + 192.D0*
c$$$     &    Hr3(-1,-1,0)*x - 16.D0*Hr3(-1,0,0) - 304.D0*Hr3(-1,0,0)*x - 
c$$$     &    288.D0*Hr3(-1,0,0)*dp - 288.D0*Hr3(-1,0,1) - 480.D0*Hr3(-1,0,
c$$$     &    1)*x - 192.D0*Hr3(-1,0,1)*dp + 96.D0*Hr3(0,-1,0) - 160.D0*
c$$$     &    Hr3(0,-1,0)*x - 96.D0*Hr3(0,-1,0)*dp + 144.D0*Hr3(0,0,0) + 
c$$$     &    176.D0*Hr3(0,0,0)*x + 144.D0*Hr3(0,0,0)*dp - 288.D0*Hr3(0,0,0
c$$$     &    )*dm + 128.D0*Hr3(0,0,1) + 288.D0*Hr3(0,0,1)*x + 96.D0*Hr3(0,
c$$$     &    0,1)*dp + 88.D0*Hr3(0,1,0) + 24.D0*Hr3(0,1,0)*x - 96.D0*Hr3(0
c$$$     &    ,1,0)*dm + 96.D0*Hr3(1,0,0) + 96.D0*Hr3(1,0,0)*x - 192.D0*
c$$$     &    Hr3(1,0,0)*dm - 128.D0*Hr4(-1,-1,0,0) + 128.D0*Hr4(-1,-1,0,0)
c$$$     &    *x )
c$$$      PTqq2 = PTqq2 + cf**3 * ( 256.D0*Hr4(-1,-1,0,0)*dp - 768.D0*Hr4(
c$$$     &    -1,-1,0,1) + 768.D0*Hr4(-1,-1,0,1)*x + 1536.D0*Hr4(-1,-1,0,1)
c$$$     &    *dp + 128.D0*Hr4(-1,0,-1,0) - 128.D0*Hr4(-1,0,-1,0)*x - 256.D0
c$$$     &    *Hr4(-1,0,-1,0)*dp - 224.D0*Hr4(-1,0,0,0) + 224.D0*Hr4(-1,0,0
c$$$     &    ,0)*x + 448.D0*Hr4(-1,0,0,0)*dp + 256.D0*Hr4(-1,0,0,1) - 256.D
c$$$     &    0*Hr4(-1,0,0,1)*x - 512.D0*Hr4(-1,0,0,1)*dp - 64.D0*Hr4(-1,0,
c$$$     &    1,0) + 64.D0*Hr4(-1,0,1,0)*x + 128.D0*Hr4(-1,0,1,0)*dp - 64.D0
c$$$     &    *Hr4(0,-1,-1,0) - 320.D0*Hr4(0,-1,-1,0)*x - 256.D0*Hr4(0,-1,
c$$$     &    -1,0)*dp + 384.D0*Hr4(0,-1,-1,0)*dm + 192.D0*Hr4(0,-1,0,0)*x
c$$$     &     + 64.D0*Hr4(0,-1,0,0)*dp - 192.D0*Hr4(0,-1,0,0)*dm + 320.D0*
c$$$     &    Hr4(0,-1,0,1) - 320.D0*Hr4(0,-1,0,1)*x - 640.D0*Hr4(0,-1,0,1)
c$$$     &    *dp - 160.D0*Hr4(0,0,-1,0) + 288.D0*Hr4(0,0,-1,0)*x + 320.D0*
c$$$     &    Hr4(0,0,-1,0)*dp - 64.D0*Hr4(0,0,-1,0)*dm - 112.D0*Hr4(0,0,0,
c$$$     &    0) - 592.D0*Hr4(0,0,0,0)*x - 320.D0*Hr4(0,0,0,0)*dp + 448.D0*
c$$$     &    Hr4(0,0,0,0)*dm - 240.D0*Hr4(0,0,0,1) - 240.D0*Hr4(0,0,0,1)*x
c$$$     &     + 448.D0*Hr4(0,0,0,1)*dm )
c$$$      PTqq2 = PTqq2 + cf**3 * (  - 128.D0*Hr4(0,0,1,0) - 192.D0*Hr4(0,0
c$$$     &    ,1,0)*x - 64.D0*Hr4(0,0,1,0)*dp + 256.D0*Hr4(0,0,1,0)*dm - 64.
c$$$     &    D0*Hr4(0,0,1,1) - 64.D0*Hr4(0,0,1,1)*x + 128.D0*Hr4(0,0,1,1)*
c$$$     &    dm - 144.D0*Hr4(0,1,0,0) - 144.D0*Hr4(0,1,0,0)*x + 192.D0*
c$$$     &    Hr4(0,1,0,0)*dm - 64.D0*Hr4(0,1,0,1) - 64.D0*Hr4(0,1,0,1)*x
c$$$     &     + 128.D0*Hr4(0,1,0,1)*dm - 64.D0*Hr4(0,1,1,0) - 64.D0*Hr4(0,
c$$$     &    1,1,0)*x + 128.D0*Hr4(0,1,1,0)*dm - 128.D0*Hr4(1,0,-1,0) - 
c$$$     &    128.D0*Hr4(1,0,-1,0)*x + 256.D0*Hr4(1,0,-1,0)*dm - 128.D0*
c$$$     &    Hr4(1,0,0,0) - 128.D0*Hr4(1,0,0,0)*x + 256.D0*Hr4(1,0,0,0)*dm
c$$$     &     - 128.D0*Hr4(1,0,0,1) - 128.D0*Hr4(1,0,0,1)*x + 256.D0*Hr4(1
c$$$     &    ,0,0,1)*dm - 64.D0*Hr4(1,0,1,0) - 64.D0*Hr4(1,0,1,0)*x + 128.D
c$$$     &    0*Hr4(1,0,1,0)*dm )
c$$$      PTqq2 = PTqq2 + nf*cf*ca * (  - 182.D0/3.D0 + 2474.D0/27.D0*x - 
c$$$     &    836.D0/27.D0*dm + 16.D0*z3 + 32.D0*z3*x + 16.D0*z3*dp - 48.D0
c$$$     &    *z3*dm + 8.D0*z2 - 184.D0/9.D0*z2*x - 160.D0/9.D0*z2*dp + 160.
c$$$     &    D0/9.D0*z2*dm + 32.D0/3.D0*Hr1(-1)*z2 - 32.D0/3.D0*Hr1(-1)*z2
c$$$     &    *x - 64.D0/3.D0*Hr1(-1)*z2*dp + 68.D0/27.D0*Hr1(0) + 1700.D0/
c$$$     &    27.D0*Hr1(0)*x - 1336.D0/27.D0*Hr1(0)*dm - 8.D0*Hr1(0)*z2 - 8.
c$$$     &    D0/3.D0*Hr1(0)*z2*x + 16.D0/3.D0*Hr1(0)*z2*dp + 32.D0/3.D0*
c$$$     &    Hr1(0)*z2*dm - 16.D0*Hr1(1) + 16.D0*Hr1(1)*x + 64.D0/9.D0*
c$$$     &    Hr2(-1,0) - 256.D0/9.D0*Hr2(-1,0)*x - 320.D0/9.D0*Hr2(-1,0)*
c$$$     &    dp + 88.D0/9.D0*Hr2(0,0) + 272.D0/9.D0*Hr2(0,0)*x + 160.D0/9.D
c$$$     &    0*Hr2(0,0)*dp - 112.D0/3.D0*Hr2(0,0)*dm - 8.D0*Hr2(0,1) - 8.D0
c$$$     &    *Hr2(0,1)*x + 16.D0/3.D0*Hr3(-1,0,0) - 16.D0/3.D0*Hr3(-1,0,0)
c$$$     &    *x - 32.D0/3.D0*Hr3(-1,0,0)*dp - 32.D0/3.D0*Hr3(-1,0,1) + 32.D
c$$$     &    0/3.D0*Hr3(-1,0,1)*x + 64.D0/3.D0*Hr3(-1,0,1)*dp + 16.D0/3.D0
c$$$     &    *Hr3(0,-1,0) - 16.D0/3.D0*Hr3(0,-1,0)*x - 32.D0/3.D0*Hr3(0,-1
c$$$     &    ,0)*dp )
c$$$      PTqq2 = PTqq2 + nf*cf*ca * ( 32.D0/3.D0*Hr3(0,0,0)*x + 32.D0/3.D0
c$$$     &    *Hr3(0,0,0)*dp - 32.D0/3.D0*Hr3(0,0,0)*dm + 8.D0*Hr3(0,0,1)
c$$$     &     - 8.D0/3.D0*Hr3(0,0,1)*x - 32.D0/3.D0*Hr3(0,0,1)*dp - 16.D0/
c$$$     &    3.D0*Hr3(0,0,1)*dm - 8.D0*Hr3(1,0,0) - 8.D0*Hr3(1,0,0)*x + 16.
c$$$     &    D0*Hr3(1,0,0)*dm )
c$$$      PTqq2 = PTqq2 + nf*cf**2 * ( 5.D0/9.D0 + 325.D0/9.D0*x - 110.D0/3.
c$$$     &    D0*dm - 32.D0/3.D0*z3 - 128.D0/3.D0*z3*x - 32.D0*z3*dp + 160.D
c$$$     &    0/3.D0*z3*dm - 160.D0/9.D0*z2 + 160.D0/9.D0*z2*x + 320.D0/9.D0
c$$$     &    *z2*dp - 64.D0/3.D0*Hr1(-1)*z2 + 64.D0/3.D0*Hr1(-1)*z2*x + 
c$$$     &    128.D0/3.D0*Hr1(-1)*z2*dp + 64.D0/3.D0*Hr1(0) - 256.D0/9.D0*
c$$$     &    Hr1(0)*x - 68.D0/3.D0*Hr1(0)*dm + 32.D0/3.D0*Hr1(0)*z2 - 32.D0
c$$$     &    /3.D0*Hr1(0)*z2*dp - 32.D0/3.D0*Hr1(0)*z2*dm + 64.D0/3.D0*
c$$$     &    Hr1(1) - 64.D0/3.D0*Hr1(1)*x - 128.D0/9.D0*Hr2(-1,0) + 512.D0/
c$$$     &    9.D0*Hr2(-1,0)*x + 640.D0/9.D0*Hr2(-1,0)*dp - 704.D0/9.D0*
c$$$     &    Hr2(0,0)*x - 320.D0/9.D0*Hr2(0,0)*dp + 496.D0/9.D0*Hr2(0,0)*
c$$$     &    dm + 32.D0/9.D0*Hr2(0,1) - 160.D0/9.D0*Hr2(0,1)*x + 320.D0/9.D
c$$$     &    0*Hr2(0,1)*dm - 64.D0/9.D0*Hr2(1,0) - 256.D0/9.D0*Hr2(1,0)*x
c$$$     &     + 320.D0/9.D0*Hr2(1,0)*dm - 32.D0/3.D0*Hr3(-1,0,0) + 32.D0/3.
c$$$     &    D0*Hr3(-1,0,0)*x + 64.D0/3.D0*Hr3(-1,0,0)*dp + 64.D0/3.D0*
c$$$     &    Hr3(-1,0,1) - 64.D0/3.D0*Hr3(-1,0,1)*x - 128.D0/3.D0*Hr3(-1,0
c$$$     &    ,1)*dp )
c$$$      PTqq2 = PTqq2 + nf*cf**2 * (  - 32.D0/3.D0*Hr3(0,-1,0) + 32.D0/3.D
c$$$     &    0*Hr3(0,-1,0)*x + 64.D0/3.D0*Hr3(0,-1,0)*dp - 24.D0*Hr3(0,0,0
c$$$     &    ) - 136.D0/3.D0*Hr3(0,0,0)*x - 64.D0/3.D0*Hr3(0,0,0)*dp + 160.
c$$$     &    D0/3.D0*Hr3(0,0,0)*dm - 64.D0/3.D0*Hr3(0,0,1) + 64.D0/3.D0*
c$$$     &    Hr3(0,0,1)*dp + 64.D0/3.D0*Hr3(0,0,1)*dm - 32.D0/3.D0*Hr3(0,1
c$$$     &    ,0) - 32.D0/3.D0*Hr3(0,1,0)*x + 64.D0/3.D0*Hr3(0,1,0)*dm )
c$$$      PTqq2 = PTqq2 + nf2*cf * ( 112.D0/27.D0 - 32.D0/9.D0*x - 16.D0/
c$$$     &    27.D0*dm + 8.D0/27.D0*Hr1(0) - 88.D0/27.D0*Hr1(0)*x + 80.D0/
c$$$     &    27.D0*Hr1(0)*dm - 8.D0/9.D0*Hr2(0,0) - 8.D0/9.D0*Hr2(0,0)*x
c$$$     &     + 16.D0/9.D0*Hr2(0,0)*dm )
c$$$*
c$$$* ...The soft (`+'-distribution) part of the splitting function
c$$$*
c$$$       A3 =
c$$$     ,      ca**2*cf * ( + 490.D0/3.D0 + 88.D0/3.D0*z3 - 1072.D0/9.D0*z2
c$$$     ,                   + 176.D0/5.D0*z2**2 )
c$$$     ,    + ca*cf*nf * ( - 836./27.D0 + 160./9.D0*z2 - 112./3.D0*z3 )
c$$$     ,    + cf**2*nf * ( - 110./3.D0 + 32.*z3 ) - cf*nf2 * 16./27.D0
c$$$*
c$$$       GQQ2L = DM * A3
c$$$*
c$$$* ...The regular piece of the splitting function
c$$$*
c$$$       X2NSPTA = PTqq2 - GQQ2L
c$$$*
c$$$       RETURN
c$$$       END
c$$$*
c$$$* ---------------------------------------------------------------------
c$$$*
c$$$*
c$$$       FUNCTION X2NSMTA (X, NF)
c$$$*
c$$$* ..This is the regular piece of P_NS-  (the odd-moment combination)
c$$$*
c$$$       IMPLICIT REAL*8 (A - Z)
c$$$       COMPLEX*16 HC1, HC2, HC3, HC4
c$$$       INTEGER NF, NF2, N1, N2, NW, I1, I2, I3, N
c$$$       PARAMETER ( N1 = -1, N2 = 1, NW = 4 )
c$$$       DIMENSION HC1(N1:N2),HC2(N1:N2,N1:N2),HC3(N1:N2,N1:N2,N1:N2),
c$$$     ,           HC4(N1:N2,N1:N2,N1:N2,N1:N2)
c$$$       DIMENSION HR1(N1:N2),HR2(N1:N2,N1:N2),HR3(N1:N2,N1:N2,N1:N2),
c$$$     ,           HR4(N1:N2,N1:N2,N1:N2,N1:N2)
c$$$       DIMENSION HI1(N1:N2),HI2(N1:N2,N1:N2),HI3(N1:N2,N1:N2,N1:N2),
c$$$     ,           HI4(N1:N2,N1:N2,N1:N2,N1:N2)
c$$$       PARAMETER ( Z2 = 1.6449 34066 84822 64365 D0,
c$$$     ,             Z3 = 1.2020 56903 15959 42854 D0,
c$$$     ,             Z5 = 1.0369 27755 14336 99263 D0 )
c$$$*
c$$$* ..The soft coefficient for use in X2NSPB and X2NSMTC
c$$$*
c$$$       COMMON / P2SOFT / A3
c$$$*
c$$$* ...Colour factors
c$$$*
c$$$       CF  = 4./3.D0
c$$$       CA  = 3.D0
c$$$       NF2 = NF*NF
c$$$*
c$$$* ...Some abbreviations
c$$$*
c$$$       DX = 1.D0/X
c$$$       DM = 1.D0/(1.D0-X)
c$$$       DP = 1.D0/(1.D0+X)
c$$$*
c$$$* ...The harmonic polylogs up to weight 4 by Gehrmann and Remiddi
c$$$*
c$$$       CALL HPLOG (X, NW, HC1,HC2,HC3,HC4, HR1,HR2,HR3,HR4,
c$$$     ,            HI1,HI2,HI3,HI4, N1, N2)
c$$$*
c$$$* ...The splitting function in terms of the harmonic polylogs
c$$$*    (without the delta(1-x) part, but with the soft contribution)
c$$$*
c$$$      PTqq2m = cf*ca**2 * ( 923.D0/27.D0 - 5333.D0/27.D0*x + 
c$$$     &    490.D0/3.D0*dm + 112.D0*z3 + 48.D0*z3*x + 88.D0*z3*dp + 16.D0
c$$$     &    *z3*dm + 1504.D0/9.D0*z2 + 224.D0/3.D0*z2*x - 1072.D0/9.D0*z2
c$$$     &    *dp - 1072.D0/9.D0*z2*dm - 242.D0/5.D0*z2**2 - 62.D0/5.D0*
c$$$     &    z2**2*x + 32.D0*z2**2*dp + 384.D0/5.D0*z2**2*dm + 192.D0*Hr1(
c$$$     &    -1)*z3 - 192.D0*Hr1(-1)*z3*x - 384.D0*Hr1(-1)*z3*dp - 208.D0/
c$$$     &    3.D0*Hr1(-1)*z2 - 560.D0/3.D0*Hr1(-1)*z2*x - 352.D0/3.D0*Hr1(
c$$$     &    -1)*z2*dp - 1378.D0/27.D0*Hr1(0) - 2266.D0/27.D0*Hr1(0)*x + 
c$$$     &    24.D0*Hr1(0)*dp + 4172.D0/27.D0*Hr1(0)*dm - 144.D0*Hr1(0)*z3
c$$$     &     + 128.D0*Hr1(0)*z3*dp + 128.D0*Hr1(0)*z3*dm + 332.D0/3.D0*
c$$$     &    Hr1(0)*z2 + 188.D0*Hr1(0)*z2*x + 16.D0/3.D0*Hr1(0)*z2*dp - 
c$$$     &    248.D0/3.D0*Hr1(0)*z2*dm - 592.D0/3.D0*Hr1(1) + 592.D0/3.D0*
c$$$     &    Hr1(1)*x - 144.D0*Hr1(1)*z3 - 144.D0*Hr1(1)*z3*x + 288.D0*
c$$$     &    Hr1(1)*z3*dm - 32.D0*Hr1(1)*z2 + 32.D0*Hr1(1)*z2*x - 256.D0*
c$$$     &    Hr2(-1,-1)*z2 + 256.D0*Hr2(-1,-1)*z2*x + 512.D0*Hr2(-1,-1)*z2
c$$$     &    *dp )
c$$$      PTqq2m = PTqq2m + cf*ca**2 * ( 832.D0/9.D0*Hr2(-1,0) - 1312.D0/9.D
c$$$     &    0*Hr2(-1,0)*x - 2144.D0/9.D0*Hr2(-1,0)*dp + 176.D0*Hr2(-1,0)*
c$$$     &    z2 - 176.D0*Hr2(-1,0)*z2*x - 352.D0*Hr2(-1,0)*z2*dp + 136.D0*
c$$$     &    Hr2(0,-1)*z2 - 136.D0*Hr2(0,-1)*z2*x - 256.D0*Hr2(0,-1)*z2*dp
c$$$     &     - 130.D0*Hr2(0,0) + 238.D0/9.D0*Hr2(0,0)*x + 1072.D0/9.D0*
c$$$     &    Hr2(0,0)*dp + 1556.D0/9.D0*Hr2(0,0)*dm - 68.D0*Hr2(0,0)*z2 + 
c$$$     &    36.D0*Hr2(0,0)*z2*x + 96.D0*Hr2(0,0)*z2*dp - 224.D0/3.D0*Hr2(
c$$$     &    0,1) - 224.D0/3.D0*Hr2(0,1)*x + 24.D0*Hr2(0,1)*z2 - 40.D0*
c$$$     &    Hr2(0,1)*z2*x - 64.D0*Hr2(0,1)*z2*dp + 16.D0*Hr2(1,0)*z2 + 16.
c$$$     &    D0*Hr2(1,0)*z2*x - 32.D0*Hr2(1,0)*z2*dm + 64.D0*Hr3(-1,-1,0)
c$$$     &     + 64.D0*Hr3(-1,-1,0)*x + 232.D0/3.D0*Hr3(-1,0,0) + 56.D0/3.D0
c$$$     &    *Hr3(-1,0,0)*x - 176.D0/3.D0*Hr3(-1,0,0)*dp + 304.D0/3.D0*
c$$$     &    Hr3(-1,0,1) + 656.D0/3.D0*Hr3(-1,0,1)*x + 352.D0/3.D0*Hr3(-1,
c$$$     &    0,1)*dp + 328.D0/3.D0*Hr3(0,-1,0) - 376.D0/3.D0*Hr3(0,-1,0)*x
c$$$     &     - 176.D0/3.D0*Hr3(0,-1,0)*dp - 48.D0*Hr3(0,-1,0)*dm - 416.D0/
c$$$     &    3.D0*Hr3(0,0,0) )
c$$$      PTqq2m = PTqq2m + cf*ca**2 * ( 248.D0/3.D0*Hr3(0,0,0)*dp + 248.D0/
c$$$     &    3.D0*Hr3(0,0,0)*dm - 4.D0/3.D0*Hr3(0,0,1) - 188.D0*Hr3(0,0,1)
c$$$     &    *x - 176.D0/3.D0*Hr3(0,0,1)*dp + 88.D0/3.D0*Hr3(0,0,1)*dm + 
c$$$     &    12.D0*Hr3(1,0,0) + 76.D0*Hr3(1,0,0)*x - 88.D0*Hr3(1,0,0)*dm
c$$$     &     + 128.D0*Hr4(-1,-1,0,0) - 128.D0*Hr4(-1,-1,0,0)*x - 256.D0*
c$$$     &    Hr4(-1,-1,0,0)*dp + 256.D0*Hr4(-1,-1,0,1) - 256.D0*Hr4(-1,-1,
c$$$     &    0,1)*x - 512.D0*Hr4(-1,-1,0,1)*dp - 48.D0*Hr4(-1,0,0,0) + 48.D
c$$$     &    0*Hr4(-1,0,0,0)*x + 96.D0*Hr4(-1,0,0,0)*dp - 128.D0*Hr4(-1,0,
c$$$     &    0,1) + 128.D0*Hr4(-1,0,0,1)*x + 256.D0*Hr4(-1,0,0,1)*dp - 48.D
c$$$     &    0*Hr4(0,-1,-1,0) - 80.D0*Hr4(0,-1,-1,0)*x + 128.D0*Hr4(0,-1,
c$$$     &    -1,0)*dm - 56.D0*Hr4(0,-1,0,0) + 88.D0*Hr4(0,-1,0,0)*x + 128.D
c$$$     &    0*Hr4(0,-1,0,0)*dp - 32.D0*Hr4(0,-1,0,0)*dm - 160.D0*Hr4(0,-1
c$$$     &    ,0,1) + 96.D0*Hr4(0,-1,0,1)*x + 256.D0*Hr4(0,-1,0,1)*dp + 64.D
c$$$     &    0*Hr4(0,-1,0,1)*dm - 40.D0*Hr4(0,0,-1,0) + 8.D0*Hr4(0,0,-1,0)
c$$$     &    *x + 32.D0*Hr4(0,0,-1,0)*dp + 32.D0*Hr4(0,0,-1,0)*dm + 40.D0*
c$$$     &    Hr4(0,0,0,0) )
c$$$      PTqq2m = PTqq2m + cf*ca**2 * (  - 32.D0*Hr4(0,0,0,0)*dp - 32.D0*
c$$$     &    Hr4(0,0,0,0)*dm + 28.D0*Hr4(0,0,0,1) - 36.D0*Hr4(0,0,0,1)*x
c$$$     &     - 64.D0*Hr4(0,0,0,1)*dp + 32.D0*Hr4(0,0,0,1)*dm + 28.D0*Hr4(
c$$$     &    0,1,0,0) + 28.D0*Hr4(0,1,0,0)*x - 64.D0*Hr4(0,1,0,0)*dm - 96.D
c$$$     &    0*Hr4(1,0,-1,0) - 96.D0*Hr4(1,0,-1,0)*x + 192.D0*Hr4(1,0,-1,0
c$$$     &    )*dm + 48.D0*Hr4(1,0,0,0) + 48.D0*Hr4(1,0,0,0)*x - 96.D0*Hr4(
c$$$     &    1,0,0,0)*dm - 64.D0*Hr4(1,0,0,1) - 64.D0*Hr4(1,0,0,1)*x + 128.
c$$$     &    D0*Hr4(1,0,0,1)*dm + 64.D0*Hr4(1,1,0,0) + 64.D0*Hr4(1,1,0,0)*
c$$$     &    x - 128.D0*Hr4(1,1,0,0)*dm )
c$$$      PTqq2m = PTqq2m + cf**2*ca * ( 1516.D0/3.D0 - 1516.D0/3.D0*x - 
c$$$     &    832.D0/3.D0*z3 - 880.D0/3.D0*z3*x - 248.D0*z3*dp + 80.D0/3.D0
c$$$     &    *z3*dm - 3880.D0/9.D0*z2 - 152.D0/9.D0*z2*x + 2144.D0/9.D0*z2
c$$$     &    *dp + 56.D0/5.D0*z2**2 - 24.D0/5.D0*z2**2*x - 8.D0*z2**2*dp
c$$$     &     - 552.D0/5.D0*z2**2*dm - 672.D0*Hr1(-1)*z3 + 672.D0*Hr1(-1)*
c$$$     &    z3*x + 1344.D0*Hr1(-1)*z3*dp + 704.D0/3.D0*Hr1(-1)*z2 + 1696.D
c$$$     &    0/3.D0*Hr1(-1)*z2*x + 992.D0/3.D0*Hr1(-1)*z2*dp - 350.D0/9.D0
c$$$     &    *Hr1(0) - 1798.D0/9.D0*Hr1(0)*x - 72.D0*Hr1(0)*dp + 302.D0/3.D
c$$$     &    0*Hr1(0)*dm + 224.D0*Hr1(0)*z3 - 240.D0*Hr1(0)*z3*dp - 240.D0
c$$$     &    *Hr1(0)*z3*dm - 184.D0*Hr1(0)*z2 - 1688.D0/3.D0*Hr1(0)*z2*x
c$$$     &     - 248.D0/3.D0*Hr1(0)*z2*dp + 536.D0/3.D0*Hr1(0)*z2*dm + 2008.
c$$$     &    D0/3.D0*Hr1(1) - 2008.D0/3.D0*Hr1(1)*x + 384.D0*Hr1(1)*z3 + 
c$$$     &    384.D0*Hr1(1)*z3*x - 768.D0*Hr1(1)*z3*dm + 112.D0*Hr1(1)*z2
c$$$     &     - 112.D0*Hr1(1)*z2*x + 896.D0*Hr2(-1,-1)*z2 - 896.D0*Hr2(-1,
c$$$     &    -1)*z2*x - 1792.D0*Hr2(-1,-1)*z2*dp - 1232.D0/9.D0*Hr2(-1,0)
c$$$     &     + 3056.D0/9.D0*Hr2(-1,0)*x )
c$$$      PTqq2m = PTqq2m + cf**2*ca * ( 4288.D0/9.D0*Hr2(-1,0)*dp - 544.D0
c$$$     &    *Hr2(-1,0)*z2 + 544.D0*Hr2(-1,0)*z2*x + 1088.D0*Hr2(-1,0)*z2*
c$$$     &    dp - 352.D0*Hr2(0,-1)*z2 + 448.D0*Hr2(0,-1)*z2*x + 768.D0*
c$$$     &    Hr2(0,-1)*z2*dp - 96.D0*Hr2(0,-1)*z2*dm + 1280.D0/9.D0*Hr2(0,
c$$$     &    0) + 776.D0/3.D0*Hr2(0,0)*x - 2144.D0/9.D0*Hr2(0,0)*dp - 3640.
c$$$     &    D0/9.D0*Hr2(0,0)*dm + 32.D0*Hr2(0,0)*z2 - 96.D0*Hr2(0,0)*z2*x
c$$$     &     - 160.D0*Hr2(0,0)*z2*dp + 128.D0*Hr2(0,0)*z2*dm + 2296.D0/9.D
c$$$     &    0*Hr2(0,1) + 4792.D0/9.D0*Hr2(0,1)*x - 2144.D0/9.D0*Hr2(0,1)*
c$$$     &    dm - 128.D0*Hr2(0,1)*z2 + 96.D0*Hr2(0,1)*z2*x + 224.D0*Hr2(0,
c$$$     &    1)*z2*dp + 64.D0*Hr2(0,1)*z2*dm - 176.D0/9.D0*Hr2(1,0) + 2320.
c$$$     &    D0/9.D0*Hr2(1,0)*x - 2144.D0/9.D0*Hr2(1,0)*dm - 32.D0*Hr2(1,0
c$$$     &    )*z2 - 32.D0*Hr2(1,0)*z2*x + 64.D0*Hr2(1,0)*z2*dm - 224.D0*
c$$$     &    Hr3(-1,-1,0) - 224.D0*Hr3(-1,-1,0)*x - 200.D0/3.D0*Hr3(-1,0,0
c$$$     &    ) - 280.D0/3.D0*Hr3(-1,0,0)*x - 80.D0/3.D0*Hr3(-1,0,0)*dp - 
c$$$     &    1040.D0/3.D0*Hr3(-1,0,1) - 2032.D0/3.D0*Hr3(-1,0,1)*x - 992.D0
c$$$     &    /3.D0*Hr3(-1,0,1)*dp )
c$$$      PTqq2m = PTqq2m + cf**2*ca * (  - 416.D0/3.D0*Hr3(0,-1,0) + 992.D0
c$$$     &    /3.D0*Hr3(0,-1,0)*x + 208.D0/3.D0*Hr3(0,-1,0)*dp + 96.D0*Hr3(
c$$$     &    0,-1,0)*dm + 652.D0/3.D0*Hr3(0,0,0) + 36.D0*Hr3(0,0,0)*x - 
c$$$     &    280.D0/3.D0*Hr3(0,0,0)*dp - 808.D0/3.D0*Hr3(0,0,0)*dm + 40.D0
c$$$     &    *Hr3(0,0,1) + 1672.D0/3.D0*Hr3(0,0,1)*x + 496.D0/3.D0*Hr3(0,0
c$$$     &    ,1)*dp - 640.D0/3.D0*Hr3(0,0,1)*dm + 80.D0/3.D0*Hr3(0,1,0) + 
c$$$     &    80.D0/3.D0*Hr3(0,1,0)*x - 352.D0/3.D0*Hr3(0,1,0)*dm + 16.D0*
c$$$     &    Hr3(1,0,0) - 112.D0*Hr3(1,0,0)*x + 96.D0*Hr3(1,0,0)*dm - 320.D
c$$$     &    0*Hr4(-1,-1,0,0) + 320.D0*Hr4(-1,-1,0,0)*x + 640.D0*Hr4(-1,-1
c$$$     &    ,0,0)*dp - 896.D0*Hr4(-1,-1,0,1) + 896.D0*Hr4(-1,-1,0,1)*x + 
c$$$     &    1792.D0*Hr4(-1,-1,0,1)*dp + 64.D0*Hr4(-1,0,-1,0) - 64.D0*Hr4(
c$$$     &    -1,0,-1,0)*x - 128.D0*Hr4(-1,0,-1,0)*dp - 16.D0*Hr4(-1,0,0,0)
c$$$     &     + 16.D0*Hr4(-1,0,0,0)*x + 32.D0*Hr4(-1,0,0,0)*dp + 384.D0*
c$$$     &    Hr4(-1,0,0,1) - 384.D0*Hr4(-1,0,0,1)*x - 768.D0*Hr4(-1,0,0,1)
c$$$     &    *dp - 32.D0*Hr4(-1,0,1,0) + 32.D0*Hr4(-1,0,1,0)*x + 64.D0*
c$$$     &    Hr4(-1,0,1,0)*dp )
c$$$      PTqq2m = PTqq2m + cf**2*ca * ( 256.D0*Hr4(0,-1,-1,0) + 192.D0*
c$$$     &    Hr4(0,-1,-1,0)*x - 128.D0*Hr4(0,-1,-1,0)*dp - 448.D0*Hr4(0,-1
c$$$     &    ,-1,0)*dm + 16.D0*Hr4(0,-1,0,0) - 176.D0*Hr4(0,-1,0,0)*x - 
c$$$     &    224.D0*Hr4(0,-1,0,0)*dp + 160.D0*Hr4(0,-1,0,0)*dm + 480.D0*
c$$$     &    Hr4(0,-1,0,1) - 352.D0*Hr4(0,-1,0,1)*x - 832.D0*Hr4(0,-1,0,1)
c$$$     &    *dp - 128.D0*Hr4(0,-1,0,1)*dm - 64.D0*Hr4(0,0,-1,0) + 64.D0*
c$$$     &    Hr4(0,0,-1,0)*x + 96.D0*Hr4(0,0,-1,0)*dp - 32.D0*Hr4(0,0,-1,0
c$$$     &    )*dm + 160.D0*Hr4(0,0,0,0) - 96.D0*Hr4(0,0,0,0)*dp - 96.D0*
c$$$     &    Hr4(0,0,0,0)*dm - 32.D0*Hr4(0,0,0,1) + 96.D0*Hr4(0,0,0,1)*x
c$$$     &     + 128.D0*Hr4(0,0,0,1)*dp - 160.D0*Hr4(0,0,0,1)*dm + 32.D0*
c$$$     &    Hr4(0,0,1,0) - 32.D0*Hr4(0,0,1,0)*dp - 32.D0*Hr4(0,0,1,0)*dm
c$$$     &     - 16.D0*Hr4(0,1,0,0) - 16.D0*Hr4(0,1,0,0)*x + 96.D0*Hr4(0,1,
c$$$     &    0,0)*dm + 256.D0*Hr4(1,0,-1,0) + 256.D0*Hr4(1,0,-1,0)*x - 512.
c$$$     &    D0*Hr4(1,0,-1,0)*dm - 80.D0*Hr4(1,0,0,0) - 80.D0*Hr4(1,0,0,0)
c$$$     &    *x + 160.D0*Hr4(1,0,0,0)*dm + 128.D0*Hr4(1,0,0,1) + 128.D0*
c$$$     &    Hr4(1,0,0,1)*x )
c$$$      PTqq2m = PTqq2m + cf**2*ca * (  - 256.D0*Hr4(1,0,0,1)*dm - 128.D0
c$$$     &    *Hr4(1,1,0,0) - 128.D0*Hr4(1,1,0,0)*x + 256.D0*Hr4(1,1,0,0)*
c$$$     &    dm )
c$$$      PTqq2m = PTqq2m + cf**3 * (  - 302.D0 + 302.D0*x + 48.D0*z3 + 336.
c$$$     &    D0*z3*x + 144.D0*z3*dp + 204.D0*z2 + 12.D0*z2*x + 504.D0/5.D0
c$$$     &    *z2**2 - 56.D0/5.D0*z2**2*x - 112.D0*z2**2*dp + 144.D0/5.D0*
c$$$     &    z2**2*dm + 576.D0*Hr1(-1)*z3 - 576.D0*Hr1(-1)*z3*x - 1152.D0*
c$$$     &    Hr1(-1)*z3*dp - 192.D0*Hr1(-1)*z2 - 384.D0*Hr1(-1)*z2*x - 192.
c$$$     &    D0*Hr1(-1)*z2*dp + 10.D0*Hr1(0) + 114.D0*Hr1(0)*x + 48.D0*
c$$$     &    Hr1(0)*dp + 6.D0*Hr1(0)*dm + 128.D0*Hr1(0)*z3 - 32.D0*Hr1(0)*
c$$$     &    z3*dp - 32.D0*Hr1(0)*z3*dm - 16.D0*Hr1(0)*z2 + 432.D0*Hr1(0)*
c$$$     &    z2*x + 144.D0*Hr1(0)*z2*dp - 144.D0*Hr1(0)*z2*dm - 400.D0*
c$$$     &    Hr1(1) + 400.D0*Hr1(1)*x - 192.D0*Hr1(1)*z3 - 192.D0*Hr1(1)*
c$$$     &    z3*x + 384.D0*Hr1(1)*z3*dm - 96.D0*Hr1(1)*z2 + 96.D0*Hr1(1)*
c$$$     &    z2*x - 768.D0*Hr2(-1,-1)*z2 + 768.D0*Hr2(-1,-1)*z2*x + 1536.D0
c$$$     &    *Hr2(-1,-1)*z2*dp - 96.D0*Hr2(-1,0) - 96.D0*Hr2(-1,0)*x + 384.
c$$$     &    D0*Hr2(-1,0)*z2 - 384.D0*Hr2(-1,0)*z2*x - 768.D0*Hr2(-1,0)*z2
c$$$     &    *dp + 160.D0*Hr2(0,-1)*z2 - 352.D0*Hr2(0,-1)*z2*x - 512.D0*
c$$$     &    Hr2(0,-1)*z2*dp )
c$$$      PTqq2m = PTqq2m + cf**3 * ( 192.D0*Hr2(0,-1)*z2*dm + 172.D0*Hr2(0
c$$$     &    ,0) - 228.D0*Hr2(0,0)*x + 52.D0*Hr2(0,0)*dm + 208.D0*Hr2(0,0)
c$$$     &    *z2 + 48.D0*Hr2(0,0)*z2*x - 64.D0*Hr2(0,0)*z2*dp - 256.D0*
c$$$     &    Hr2(0,0)*z2*dm - 12.D0*Hr2(0,1) - 300.D0*Hr2(0,1)*x + 96.D0*
c$$$     &    Hr2(0,1)*z2 - 96.D0*Hr2(0,1)*z2*x - 192.D0*Hr2(0,1)*z2*dp + 
c$$$     &    144.D0*Hr2(1,0) - 144.D0*Hr2(1,0)*x + 192.D0*Hr3(-1,-1,0) + 
c$$$     &    192.D0*Hr3(-1,-1,0)*x - 176.D0*Hr3(-1,0,0) + 112.D0*Hr3(-1,0,
c$$$     &    0)*x + 288.D0*Hr3(-1,0,0)*dp + 288.D0*Hr3(-1,0,1) + 480.D0*
c$$$     &    Hr3(-1,0,1)*x + 192.D0*Hr3(-1,0,1)*dp - 160.D0*Hr3(0,-1,0) - 
c$$$     &    160.D0*Hr3(0,-1,0)*x + 96.D0*Hr3(0,-1,0)*dp + 528.D0*Hr3(0,0,
c$$$     &    0) + 240.D0*Hr3(0,0,0)*x - 144.D0*Hr3(0,0,0)*dp - 288.D0*Hr3(
c$$$     &    0,0,0)*dm + 160.D0*Hr3(0,0,1) - 256.D0*Hr3(0,0,1)*x - 96.D0*
c$$$     &    Hr3(0,0,1)*dp + 152.D0*Hr3(0,1,0) + 88.D0*Hr3(0,1,0)*x - 96.D0
c$$$     &    *Hr3(0,1,0)*dm + 96.D0*Hr3(1,0,0) + 96.D0*Hr3(1,0,0)*x - 192.D
c$$$     &    0*Hr3(1,0,0)*dm + 128.D0*Hr4(-1,-1,0,0) - 128.D0*Hr4(-1,-1,0,
c$$$     &    0)*x )
c$$$      PTqq2m = PTqq2m + cf**3 * (  - 256.D0*Hr4(-1,-1,0,0)*dp + 768.D0*
c$$$     &    Hr4(-1,-1,0,1) - 768.D0*Hr4(-1,-1,0,1)*x - 1536.D0*Hr4(-1,-1,
c$$$     &    0,1)*dp - 128.D0*Hr4(-1,0,-1,0) + 128.D0*Hr4(-1,0,-1,0)*x + 
c$$$     &    256.D0*Hr4(-1,0,-1,0)*dp + 224.D0*Hr4(-1,0,0,0) - 224.D0*Hr4(
c$$$     &    -1,0,0,0)*x - 448.D0*Hr4(-1,0,0,0)*dp - 256.D0*Hr4(-1,0,0,1)
c$$$     &     + 256.D0*Hr4(-1,0,0,1)*x + 512.D0*Hr4(-1,0,0,1)*dp + 64.D0*
c$$$     &    Hr4(-1,0,1,0) - 64.D0*Hr4(-1,0,1,0)*x - 128.D0*Hr4(-1,0,1,0)*
c$$$     &    dp - 320.D0*Hr4(0,-1,-1,0) - 64.D0*Hr4(0,-1,-1,0)*x + 256.D0*
c$$$     &    Hr4(0,-1,-1,0)*dp + 384.D0*Hr4(0,-1,-1,0)*dm + 192.D0*Hr4(0,
c$$$     &    -1,0,0) - 64.D0*Hr4(0,-1,0,0)*dp - 192.D0*Hr4(0,-1,0,0)*dm - 
c$$$     &    320.D0*Hr4(0,-1,0,1) + 320.D0*Hr4(0,-1,0,1)*x + 640.D0*Hr4(0,
c$$$     &    -1,0,1)*dp + 288.D0*Hr4(0,0,-1,0) - 160.D0*Hr4(0,0,-1,0)*x - 
c$$$     &    320.D0*Hr4(0,0,-1,0)*dp - 64.D0*Hr4(0,0,-1,0)*dm - 592.D0*
c$$$     &    Hr4(0,0,0,0) - 112.D0*Hr4(0,0,0,0)*x + 320.D0*Hr4(0,0,0,0)*dp
c$$$     &     + 448.D0*Hr4(0,0,0,0)*dm - 240.D0*Hr4(0,0,0,1) - 240.D0*Hr4(
c$$$     &    0,0,0,1)*x )
c$$$      PTqq2m = PTqq2m + cf**3 * ( 448.D0*Hr4(0,0,0,1)*dm - 192.D0*Hr4(0
c$$$     &    ,0,1,0) - 128.D0*Hr4(0,0,1,0)*x + 64.D0*Hr4(0,0,1,0)*dp + 256.
c$$$     &    D0*Hr4(0,0,1,0)*dm - 64.D0*Hr4(0,0,1,1) - 64.D0*Hr4(0,0,1,1)*
c$$$     &    x + 128.D0*Hr4(0,0,1,1)*dm - 144.D0*Hr4(0,1,0,0) - 144.D0*
c$$$     &    Hr4(0,1,0,0)*x + 192.D0*Hr4(0,1,0,0)*dm - 64.D0*Hr4(0,1,0,1)
c$$$     &     - 64.D0*Hr4(0,1,0,1)*x + 128.D0*Hr4(0,1,0,1)*dm - 64.D0*Hr4(
c$$$     &    0,1,1,0) - 64.D0*Hr4(0,1,1,0)*x + 128.D0*Hr4(0,1,1,0)*dm - 
c$$$     &    128.D0*Hr4(1,0,-1,0) - 128.D0*Hr4(1,0,-1,0)*x + 256.D0*Hr4(1,
c$$$     &    0,-1,0)*dm - 128.D0*Hr4(1,0,0,0) - 128.D0*Hr4(1,0,0,0)*x + 
c$$$     &    256.D0*Hr4(1,0,0,0)*dm - 128.D0*Hr4(1,0,0,1) - 128.D0*Hr4(1,0
c$$$     &    ,0,1)*x + 256.D0*Hr4(1,0,0,1)*dm - 64.D0*Hr4(1,0,1,0) - 64.D0
c$$$     &    *Hr4(1,0,1,0)*x + 128.D0*Hr4(1,0,1,0)*dm )
c$$$      PTqq2m = PTqq2m + nf*cf*ca * (  - 1034.D0/9.D0 + 3938.D0/27.D0*x
c$$$     &     - 836.D0/27.D0*dm + 32.D0*z3 + 16.D0*z3*x - 16.D0*z3*dp - 48.
c$$$     &    D0*z3*dm - 88.D0/9.D0*z2 - 8.D0/3.D0*z2*x + 160.D0/9.D0*z2*dp
c$$$     &     + 160.D0/9.D0*z2*dm - 32.D0/3.D0*Hr1(-1)*z2 + 32.D0/3.D0*
c$$$     &    Hr1(-1)*z2*x + 64.D0/3.D0*Hr1(-1)*z2*dp - 916.D0/27.D0*Hr1(0)
c$$$     &     + 716.D0/27.D0*Hr1(0)*x - 1336.D0/27.D0*Hr1(0)*dm - 8.D0/3.D0
c$$$     &    *Hr1(0)*z2 - 8.D0*Hr1(0)*z2*x - 16.D0/3.D0*Hr1(0)*z2*dp + 32.D
c$$$     &    0/3.D0*Hr1(0)*z2*dm + 16.D0/3.D0*Hr1(1) - 16.D0/3.D0*Hr1(1)*x
c$$$     &     - 64.D0/9.D0*Hr2(-1,0) + 256.D0/9.D0*Hr2(-1,0)*x + 320.D0/9.D
c$$$     &    0*Hr2(-1,0)*dp + 104.D0/9.D0*Hr2(0,0) - 32.D0/9.D0*Hr2(0,0)*x
c$$$     &     - 160.D0/9.D0*Hr2(0,0)*dp - 112.D0/3.D0*Hr2(0,0)*dm + 8.D0/3.
c$$$     &    D0*Hr2(0,1) + 8.D0/3.D0*Hr2(0,1)*x - 16.D0/3.D0*Hr3(-1,0,0)
c$$$     &     + 16.D0/3.D0*Hr3(-1,0,0)*x + 32.D0/3.D0*Hr3(-1,0,0)*dp + 32.D
c$$$     &    0/3.D0*Hr3(-1,0,1) - 32.D0/3.D0*Hr3(-1,0,1)*x - 64.D0/3.D0*
c$$$     &    Hr3(-1,0,1)*dp - 16.D0/3.D0*Hr3(0,-1,0) + 16.D0/3.D0*Hr3(0,-1
c$$$     &    ,0)*x )
c$$$      PTqq2m = PTqq2m + nf*cf*ca * ( 32.D0/3.D0*Hr3(0,-1,0)*dp + 32.D0/
c$$$     &    3.D0*Hr3(0,0,0) - 32.D0/3.D0*Hr3(0,0,0)*dp - 32.D0/3.D0*Hr3(0
c$$$     &    ,0,0)*dm - 8.D0/3.D0*Hr3(0,0,1) + 8.D0*Hr3(0,0,1)*x + 32.D0/3.
c$$$     &    D0*Hr3(0,0,1)*dp - 16.D0/3.D0*Hr3(0,0,1)*dm - 8.D0*Hr3(1,0,0)
c$$$     &     - 8.D0*Hr3(1,0,0)*x + 16.D0*Hr3(1,0,0)*dm )
c$$$      PTqq2m = PTqq2m + nf*cf**2 * ( 109.D0 - 217.D0/3.D0*x - 110.D0/3.D
c$$$     &    0*dm - 128.D0/3.D0*z3 - 32.D0/3.D0*z3*x + 32.D0*z3*dp + 160.D0
c$$$     &    /3.D0*z3*dm + 160.D0/9.D0*z2 - 160.D0/9.D0*z2*x - 320.D0/9.D0
c$$$     &    *z2*dp + 64.D0/3.D0*Hr1(-1)*z2 - 64.D0/3.D0*Hr1(-1)*z2*x - 
c$$$     &    128.D0/3.D0*Hr1(-1)*z2*dp + 848.D0/9.D0*Hr1(0) + 400.D0/9.D0*
c$$$     &    Hr1(0)*x - 68.D0/3.D0*Hr1(0)*dm + 32.D0/3.D0*Hr1(0)*z2*x + 32.
c$$$     &    D0/3.D0*Hr1(0)*z2*dp - 32.D0/3.D0*Hr1(0)*z2*dm - 64.D0/3.D0*
c$$$     &    Hr1(1) + 64.D0/3.D0*Hr1(1)*x + 128.D0/9.D0*Hr2(-1,0) - 512.D0/
c$$$     &    9.D0*Hr2(-1,0)*x - 640.D0/9.D0*Hr2(-1,0)*dp - 32.D0/9.D0*Hr2(
c$$$     &    0,0) - 32.D0/3.D0*Hr2(0,0)*x + 320.D0/9.D0*Hr2(0,0)*dp + 496.D
c$$$     &    0/9.D0*Hr2(0,0)*dm - 160.D0/9.D0*Hr2(0,1) - 352.D0/9.D0*Hr2(0
c$$$     &    ,1)*x + 320.D0/9.D0*Hr2(0,1)*dm - 64.D0/9.D0*Hr2(1,0) - 256.D0
c$$$     &    /9.D0*Hr2(1,0)*x + 320.D0/9.D0*Hr2(1,0)*dm + 32.D0/3.D0*Hr3(
c$$$     &    -1,0,0) - 32.D0/3.D0*Hr3(-1,0,0)*x - 64.D0/3.D0*Hr3(-1,0,0)*
c$$$     &    dp - 64.D0/3.D0*Hr3(-1,0,1) + 64.D0/3.D0*Hr3(-1,0,1)*x + 128.D
c$$$     &    0/3.D0*Hr3(-1,0,1)*dp )
c$$$      PTqq2m = PTqq2m + nf*cf**2 * ( 32.D0/3.D0*Hr3(0,-1,0) - 32.D0/3.D0
c$$$     &    *Hr3(0,-1,0)*x - 64.D0/3.D0*Hr3(0,-1,0)*dp - 136.D0/3.D0*Hr3(
c$$$     &    0,0,0) - 24.D0*Hr3(0,0,0)*x + 64.D0/3.D0*Hr3(0,0,0)*dp + 160.D
c$$$     &    0/3.D0*Hr3(0,0,0)*dm - 64.D0/3.D0*Hr3(0,0,1)*x - 64.D0/3.D0*
c$$$     &    Hr3(0,0,1)*dp + 64.D0/3.D0*Hr3(0,0,1)*dm - 32.D0/3.D0*Hr3(0,1
c$$$     &    ,0) - 32.D0/3.D0*Hr3(0,1,0)*x + 64.D0/3.D0*Hr3(0,1,0)*dm )
c$$$      PTqq2m = PTqq2m + nf2*cf * ( 112.D0/27.D0 - 32.D0/9.D0*x - 16.D0
c$$$     &    /27.D0*dm + 8.D0/27.D0*Hr1(0) - 88.D0/27.D0*Hr1(0)*x + 80.D0/
c$$$     &    27.D0*Hr1(0)*dm - 8.D0/9.D0*Hr2(0,0) - 8.D0/9.D0*Hr2(0,0)*x
c$$$     &     + 16.D0/9.D0*Hr2(0,0)*dm )
c$$$*
c$$$* ...The soft (`+'-distribution) part of the splitting function
c$$$*
c$$$       A3 =
c$$$     ,      ca**2*cf * ( + 490.D0/3.D0 + 88.D0/3.D0*z3 - 1072.D0/9.D0*z2
c$$$     ,                   + 176.D0/5.D0*z2**2 )
c$$$     ,    + ca*cf*nf * ( - 836./27.D0 + 160./9.D0*z2 - 112./3.D0*z3 )
c$$$     ,    + cf**2*nf * ( - 110./3.D0 + 32.*z3 ) - cf*nf2 * 16./27.D0
c$$$*
c$$$       GQQ2L = DM * A3
c$$$*
c$$$* ...The regular piece of the splitting function
c$$$*
c$$$       X2NSMTA = PTqq2m - GQQ2L
c$$$*
c$$$       RETURN
c$$$       END
c$$$*
c$$$* ---------------------------------------------------------------------
c$$$*
c$$$*
c$$$* ..This is the singular (soft) piece  (as in the spacelike case)
c$$$*
c$$$       FUNCTION X2NSTB (Y, NF)
c$$$       IMPLICIT REAL*8 (A - Z)
c$$$       INTEGER NF
c$$$*
c$$$       COMMON / P2SOFT / A3
c$$$*
c$$$       X2NSTB  = A3/(1.D0-Y)
c$$$*
c$$$       RETURN
c$$$       END
c$$$*
c$$$* ---------------------------------------------------------------------
c$$$*
c$$$*
c$$$* ..This is the 'local' piece  (as in the spacelike case)
c$$$*
c$$$       FUNCTION X2NSTC (Y, NF)
c$$$*
c$$$       IMPLICIT REAL*8 (A - Z)
c$$$       INTEGER NF, NF2
c$$$       PARAMETER ( Z2 = 1.6449 34066 84822 64365 D0,
c$$$     ,             Z3 = 1.2020 56903 15959 42854 D0,
c$$$     ,             Z5 = 1.0369 27755 14336 99263 D0 )
c$$$*
c$$$       COMMON / P2SOFT / A3
c$$$*
c$$$* ...Colour factors
c$$$*
c$$$       CF  = 4./3.D0
c$$$       CA  = 3.D0
c$$$       NF2 = NF*NF
c$$$*
c$$$* ...The coefficient of delta(1-x)
c$$$*
c$$$       P2DELT = 
c$$$     &     + 29.D0/2.D0*cf**3
c$$$     &     + 151.D0/4.D0*ca*cf**2
c$$$     &     - 1657.D0/36.D0*ca**2*cf
c$$$     &     - 240.D0*z5*cf**3
c$$$     &     + 120.D0*z5*ca*cf**2
c$$$     &     + 40.D0*z5*ca**2*cf
c$$$     &     + 68.D0*z3*cf**3
c$$$     &     + 844.D0/3.D0*z3*ca*cf**2
c$$$     &     - 1552.D0/9.D0*z3*ca**2*cf
c$$$     &     + 18.D0*z2*cf**3
c$$$     &     - 410.D0/3.D0*z2*ca*cf**2
c$$$     &     + 4496.D0/27.D0*z2*ca**2*cf
c$$$     &     - 32.D0*z2*z3*cf**3
c$$$     &     + 16.D0*z2*z3*ca*cf**2
c$$$     &     + 288.D0/5.D0*z2**2*cf**3
c$$$     &     - 988.D0/15.D0*z2**2*ca*cf**2
c$$$     &     - 2.D0*z2**2*ca**2*cf
c$$$*
c$$$     &     - 1336.D0/27.D0*z2*ca*cf*nf
c$$$     &     + 4.D0/5.D0*z2**2*ca*cf*nf
c$$$     &     + 200.D0/9.D0*z3*ca*cf*nf
c$$$     &     + 20.D0*ca*cf*nf
c$$$     &     + 20.D0/3.D0*z2*cf**2*nf
c$$$     &     + 232.D0/15.D0*z2**2*cf**2*nf
c$$$     &     - 136.D0/3.D0*z3*cf**2*nf
c$$$     &     - 23.D0*cf**2*nf
c$$$     &     + 80.D0/27.D0*z2*cf*nf2
c$$$     &     - 16.D0/9.D0*z3*cf*nf2
c$$$     &     - 17.D0/9.D0*cf*nf2
c$$$*
c$$$       X2NSTC = LOG (1.D0-Y) * A3 + P2DELT
c$$$*
c$$$       RETURN
c$$$       END
c$$$*
c$$$* ---------------------------------------------------------------------
c$$$*
c$$$*
c$$$* ..This is P_NSS, the difference of P_NSV and P_NS-  (odd moments).
c$$$*    Identical to the spacelike result, given here for completeness.
c$$$*
c$$$       FUNCTION X2NSSTA (X, NF)
c$$$*
c$$$       IMPLICIT REAL*8 (A - Z)
c$$$       COMPLEX*16 HC1, HC2, HC3, HC4 
c$$$       INTEGER NF, NF2, N1, N2, NW, I1, I2, I3, N
c$$$       PARAMETER ( N1 = -1, N2 = 1, NW = 4 ) 
c$$$       DIMENSION HC1(N1:N2),HC2(N1:N2,N1:N2),HC3(N1:N2,N1:N2,N1:N2), 
c$$$     ,           HC4(N1:N2,N1:N2,N1:N2,N1:N2) 
c$$$       DIMENSION HR1(N1:N2),HR2(N1:N2,N1:N2),HR3(N1:N2,N1:N2,N1:N2), 
c$$$     ,           HR4(N1:N2,N1:N2,N1:N2,N1:N2) 
c$$$       DIMENSION HI1(N1:N2),HI2(N1:N2,N1:N2),HI3(N1:N2,N1:N2,N1:N2), 
c$$$     ,           HI4(N1:N2,N1:N2,N1:N2,N1:N2) 
c$$$       PARAMETER ( Z2 = 1.6449 34066 84822 64365 D0,
c$$$     ,             Z3 = 1.2020 56903 15959 42854 D0 )
c$$$*
c$$$* ...An abbreviation
c$$$*
c$$$       DX = 1.D0/X
c$$$*
c$$$* ...The harmonic polylogs up to weight 4 by Gehrmann and Remiddi
c$$$*
c$$$       CALL HPLOG (X, NW, HC1,HC2,HC3,HC4, HR1,HR2,HR3,HR4,
c$$$     ,            HI1,HI2,HI3,HI4, N1, N2) 
c$$$*
c$$$* ...The splitting function in terms of the harmonic polylogs
c$$$*
c$$$      gqq2 =
c$$$     &  + 5./18.D0 * NF * ( 6400.D0/3.D0 - 6400.D0/3.D0*x + 256.D0*z3
c$$$     &     + 1280.D0/3.D0*z3*x**2 - 2144.D0/3.D0*z2 - 1312.D0/3.D0*z2*x
c$$$     &     + 96.D0*z2**2 + 160.D0*z2**2*x - 192.D0*Hr1(-1)*z2 - 192.D0*
c$$$     &    Hr1(-1)*z2*x - 256.D0*Hr1(-1)*z2*x**2 - 256.D0*Hr1(-1)*z2*dx
c$$$     &     + 3200.D0/3.D0*Hr1(0) + 96.D0*Hr1(0)*x - 256.D0*Hr1(0)*z3 +
c$$$     &    32.D0*Hr1(0)*z2 + 288.D0*Hr1(0)*z2*x + 1024.D0/3.D0*Hr1(0)*z2
c$$$     &    *x**2 + 2912.D0/3.D0*Hr1(1) - 2912.D0/3.D0*Hr1(1)*x - 64.D0*
c$$$     &    Hr1(1)*z2 + 64.D0*Hr1(1)*z2*x + 256.D0/3.D0*Hr1(1)*z2*x**2 -
c$$$     &    256.D0/3.D0*Hr1(1)*z2*dx - 832.D0/3.D0*Hr2(-1,0) - 832.D0/3.D0
c$$$     &    *Hr2(-1,0)*x + 128.D0*Hr2(0,-1)*z2 - 128.D0*Hr2(0,-1)*z2*x +
c$$$     &    1216.D0/3.D0*Hr2(0,0) + 928.D0/3.D0*Hr2(0,0)*x - 320.D0*Hr2(0
c$$$     &    ,0)*z2 - 192.D0*Hr2(0,0)*z2*x + 1312.D0/3.D0*Hr2(0,1) + 1312.D
c$$$     &    0/3.D0*Hr2(0,1)*x - 128.D0*Hr2(0,1)*z2 - 128.D0*Hr2(0,1)*z2*x
c$$$     &     + 128.D0*Hr3(-1,-1,0) + 128.D0*Hr3(-1,-1,0)*x - 512.D0/3.D0*
c$$$     &    Hr3(-1,-1,0)*x**2 - 512.D0/3.D0*Hr3(-1,-1,0)*dx + 64.D0*Hr3(
c$$$     &    -1,0,0) )
c$$$      gqq2 = gqq2 + 5./18.D0 * NF* ( 64.D0*Hr3(-1,0,0)*x + 512.D0/3.D
c$$$     &    0*Hr3(-1,0,0)*x**2 + 512.D0/3.D0*Hr3(-1,0,0)*dx + 256.D0*Hr3(
c$$$     &    -1,0,1) + 256.D0*Hr3(-1,0,1)*x + 512.D0/3.D0*Hr3(-1,0,1)*x**2
c$$$     &     + 512.D0/3.D0*Hr3(-1,0,1)*dx + 64.D0*Hr3(0,-1,0) - 192.D0*
c$$$     &    Hr3(0,-1,0)*x + 512.D0/3.D0*Hr3(0,-1,0)*x**2 - 64.D0*Hr3(0,0,
c$$$     &    0) - 512.D0/3.D0*Hr3(0,0,0)*x**2 + 32.D0*Hr3(0,0,1) - 288.D0*
c$$$     &    Hr3(0,0,1)*x - 512.D0/3.D0*Hr3(0,0,1)*x**2 - 96.D0*Hr3(1,0,0)
c$$$     &     + 96.D0*Hr3(1,0,0)*x + 256.D0*Hr4(0,-1,-1,0) - 256.D0*Hr4(0,
c$$$     &    -1,-1,0)*x - 128.D0*Hr4(0,-1,0,0) + 128.D0*Hr4(0,-1,0,0)*x -
c$$$     &    128.D0*Hr4(0,0,-1,0) + 128.D0*Hr4(0,0,-1,0)*x + 128.D0*Hr4(0,
c$$$     &    0,0,0) + 192.D0*Hr4(0,0,0,1) + 192.D0*Hr4(0,0,0,1)*x - 64.D0*
c$$$     &    Hr4(0,1,0,0) - 64.D0*Hr4(0,1,0,0)*x )
c$$$*
c$$$       X2NSSTA = GQQ2 
c$$$*
c$$$       RETURN
c$$$       END
c$$$*
c$$$* =================================================================av==

