************************************************************************
*
*     SplittingFunctions.f:
*
*     Collections of the QCD splitting fuctions up to NNLO
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
* ..The soft coefficient for use in X2NSPB and X2NSPC
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
* ..The soft coefficient for use in X2NSPB and X2NSPC
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
