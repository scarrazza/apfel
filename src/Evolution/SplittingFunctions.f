***********************************************************************
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
      include "../commons/ColorFactors.h"
*
      INTEGER NF
*
* ..The soft coefficient for use in X2NSB and X2NSC
*
      COMMON / P1SOFT / A2
*
* ...some abbreviations
*
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
     3     / 6d0 + lnx**2 / 2d0 - pi**2 / 6d0 ) * pqq
     4     + ( 20d0 * ( 1d0 - x ) ) / 3d0 + lnx * ( 1d0 + x) )
     5     + 4d0 * CF**2 * ( ( ( - 3d0 * lnx ) / 2d0 - 2d0 * ln1mx 
     6     * lnx ) * pqq - 5d0 * ( 1d0 - x ) - ( lnx**2 * ( 1d0 
     7     + x ) ) / 2d0 - lnx * ( 1.5d0 + ( 7d0 * x ) / 2d0 ) )
     8     + 4d0 * CF * ( CF - CA / 2d0 ) * ( 2d0 * pqqmx * S2x 
     9     + 4d0 * ( 1d0 - x ) + 2d0 * lnx * ( 1d0 + x ) )
*
* ...The soft (`+'-distribution) part of the splitting function
*
       A2 = - 40D0/9D0*cf*nf + 268D0/9D0*ca*cf - 8D0*zeta2*ca*cf
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
      include "../commons/ColorFactors.h"
*
      INTEGER NF
*
* ..The soft coefficient for use in X2NSB and X2NSC
*
      COMMON / P1SOFT / A2
*
* ...some abbreviations
*
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
     3     / 6d0 + lnx**2 / 2d0 - pi**2 / 6d0 ) * pqq
     4     + ( 20d0 * ( 1d0 - x ) ) / 3d0 + lnx * ( 1d0 + x) )
     5     + 4d0 * CF**2 * ( ( ( - 3d0 * lnx ) / 2d0 - 2d0 * ln1mx 
     6     * lnx ) * pqq - 5d0 * ( 1d0 - x ) - ( lnx**2 * ( 1d0 
     7     + x ) ) / 2d0 - lnx * ( 1.5d0 + ( 7d0 * x ) / 2d0 ) )
     8     - 4d0 * CF * ( CF - CA / 2d0 ) * ( 2d0 * pqqmx * S2x 
     9     + 4d0 * ( 1d0 - x ) + 2d0 * lnx * ( 1d0 + x ) )
*
* ...The soft (`+'-distribution) part of the splitting function
*
       A2 = - 40D0/9D0*cf*nf + 268D0/9D0*ca*cf - 8D0*zeta2*ca*cf
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
*
       include "../commons/consts.h"
       include "../commons/ColorFactors.h"
*
       INTEGER NF
*
* ...The coefficient of delta(1-x)
*
       P1DELT = 
     &     - 1D0/3D0*cf*nf
     &     + 3D0/2D0*cf**2
     &     + 17D0/6D0*ca*cf
     &     + 24D0*zeta3*cf**2
     &     - 12D0*zeta3*ca*cf
     &     - 8D0/3D0*zeta2*cf*nf
     &     - 12D0*zeta2*cf**2
     &     + 44D0/3D0*zeta2*ca*cf
*
* ...The soft (`+'-distribution) part of the splitting function
*
       A2 = - 40D0/9D0*cf*nf + 268D0/9D0*ca*cf - 8D0*zeta2*ca*cf
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
       include "../commons/ColorFactors.h"
*
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
       include "../commons/ColorFactors.h"
*
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
       include "../commons/ColorFactors.h"
*
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
*
      include "../commons/ColorFactors.h"
*
      INTEGER NF
*
* ...some abbreviations
*
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
      include "../commons/ColorFactors.h"
*
      INTEGER NF
*
      lnx = dlog(x)
      ln1mx = dlog(1d0 - x)
      pqg   = x**2 + ( 1d0 - x )**2
      pqgmx = x**2 + ( 1d0 + x )**2
      S2x   = S2(x)
*
      X1QGA = 2d0 * CF * NF * ( 4d0  + 4d0 * ln1mx + ( 10d0 - 4d0 
     1      * ( ln1mx - lnx ) + 2d0 * ( - ln1mx + lnx )**2 
     2      - 2d0 * pi**2 / 3d0 ) * pqg - lnx * ( 1d0 - 4d0 * x )
     3      - lnx**2 * ( 1d0  - 2d0 * x ) - 9d0 * x )
     4      + 2d0 * CA * NF * ( 20.22222222222222d0 - 4d0 * ln1mx
     5      + ( - 24.22222222222222d0 + 4d0 * ln1mx - 2d0 * ln1mx**2
     6      + ( 44d0 * lnx ) /3d0 - lnx**2 + pi**2 / 3d0 ) * pqg 
     7      + 2d0 * pqgmx * S2x + 40d0 / ( 9d0 * x ) 
     8      + ( 14d0 * x ) / 9d0 - lnx**2 * ( 2d0 + 8d0 * x ) 
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
      include "../commons/ColorFactors.h"
*
      INTEGER NF
*
* ...some abbreviations
*
      lnx = dlog(x)
      ln1mx = dlog(1d0 - x)
      pgq   = ( 1d0 + ( 1d0 - x )**2 ) / x
      pgqmx = - ( 1d0 + ( 1d0 + x )**2 ) / x
      S2x   = S2(x)
*
      X1GQA = 2d0 * CF * NF * ( - ( ( 2.2222222222222223d0 
     1      + ( 4d0 * ln1mx ) / 3d0 ) * pgq ) - ( 4d0 * x ) / 3d0 )
     2      + 4d0 * CF**2 * ( - 2.5d0 - ( 3d0 * ln1mx + ln1mx**2 )
     3      * pgq - lnx**2 * ( 1d0 - x / 2d0 ) - ( 7d0 * x ) / 2d0 
     4      - 2d0 * ln1mx * x + lnx * ( 2d0 + ( 7d0 * x ) / 2d0 ) ) 
     5      + 4d0 * CA * CF * ( 3.111111111111111d0 + pgq * ( 0.5 d0 
     6      + ( 11d0 * ln1mx ) / 3d0 + ln1mx**2 - 2d0 * ln1mx * lnx
     7      + lnx**2 / 2d0 - pi**2 / 6d0 ) + pgqmx * S2x 
     8      + ( 65d0 * x ) / 18d0 + 2d0 * ln1mx * x + ( 44d0 * x**2 ) 
     9      / 9d0 + lnx**2 * ( 4d0 + x ) - lnx * ( 12d0 + 5d0 * x 
     1      + ( 8d0 * x**2 ) / 3d0 ) )
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
      include "../commons/ColorFactors.h"
*
      INTEGER NF
*
* ..The soft coefficient for use in X1GGB and X1GGC
*
       COMMON / P1GSOFT / A2G
*
* ...some abbreviations
*
      lnx = dlog(x)
      ln1mx = dlog(1d0 - x)
      pgg   = ( 1d0 / ( 1d0 - x ) +  1d0 / x - 2d0 + x * ( 1d0 - x ) )
      pggmx = ( 1d0 / ( 1d0 + x ) -  1d0 / x - 2d0 - x * ( 1d0 + x ) )
      S2x   = S2(x)
      DM    = 1D0/(1D0-X)
*
      ggg1  = 2d0 * CF * NF * ( - 16d0 + 4d0 / ( 3d0 * x ) + 8d0 * x 
     1      + ( 20d0 * x**2 ) / 3d0 - lnx**2 * ( 2d0 + 2d0 * x ) 
     2      - lnx * ( 6d0 + 10d0 * x ) )
     3      + 2d0 * CA * NF * ( 2d0 - ( 20d0 * pgg ) / 9d0 - 2d0 * x 
     4      - ( 4d0 * lnx * ( 1d0 + x ) ) / 3d0 + ( 26d0 * ( 
     5      - ( 1d0 / x ) + x**2 ) ) / 9d0 ) 
     6      + 4d0 * CA**2 * ( pgg * ( 7.444444444444445d0 
     7      - 4d0 * ln1mx * lnx + lnx**2 - pi**2 / 3d0 ) 
     8      + 2d0 * pggmx * S2x + ( 27d0 * ( 1d0 - x ) ) / 2d0 
     9      + 4d0 * lnx**2 * ( 1d0 + x ) + ( 67d0 * ( - ( 1d0 / x ) 
     1      + x**2 ) ) / 9d0 - lnx * ( 8.333333333333334d0 
     2      - ( 11d0 * x ) / 3d0 + ( 44d0 * x**2 ) / 3d0 ) )
*
* ...The soft (`+'-distribution) part of the splitting function
*
       A2G = - 40D0/9D0*ca*nf + 268D0/9D0*ca**2 - 8D0*zeta2*ca**2
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
*
       include "../commons/consts.h"
       include "../commons/ColorFactors.h"
*
       INTEGER NF
*
* ...The coefficient of delta(1-x)
*
       P1DELT = 
     ,    - 2D0*cf*nf
     ,    - 8D0/3D0*ca*nf
     ,    + 32D0/3D0*ca**2
     ,    + 12D0*zeta3*ca**2
*
* ...The soft (`+'-distribution) part of the splitting function
*
       A2G = - 40D0/9D0*ca*nf + 268D0/9D0*ca**2 - 8D0*zeta2*ca**2
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
       X0QGA = 2D0 * NF * ( 1. - 2. * X + 2. * X**2 )
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
       include "../commons/ColorFactors.h"
*
       X0GQA = 4D0 * CF * ( - 1. + 0.5 * X + 1D0 / X )
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
       include "../commons/ColorFactors.h"
*
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
       include "../commons/ColorFactors.h"
*
       X0GGB = 4D0 * CA / ( 1D0 - X )
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
*
       include "../commons/ColorFactors.h"
*
       INTEGER NF
*
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
      S2 = - 2d0 * ddilog(-x) + lnx**2 / 2d0 
     1     - 2d0 * lnx * log(1d0+x) - pi**2 / 6d0
c      S2 = - 2d0 * wgplg(1,1,-x) + lnx**2 / 2d0 
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
      include "../commons/ColorFactors.h"
*
      INTEGER NF
*
* ..The soft coefficient for use in X2NSBT and X2NSCT
*
      COMMON / P1SOFTT / A2T
*
* ...some abbreviations
*
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
     3     / 6d0 + lnx**2 / 2d0 - pi**2 / 6d0 ) * pqq
     4     + ( 20d0 * ( 1d0 - x ) ) / 3d0 + lnx * ( 1d0 + x ) )
     5     + 4d0 * CF**2 * ( ( ( 3d0 * lnx ) / 2d0 + 2d0 * ln1mx 
     6     * lnx - 2d0 * lnx**2 ) * pqq - 5d0 * ( 1d0 - x ) 
     7     + ( lnx**2 * ( 1d0 + x ) ) / 2d0 
     8     - lnx * ( 3.5d0 + ( 3d0 * x ) / 2d0 ) )
     9     + 4d0 * CF * ( CF - CA / 2d0 ) * ( 2d0 * pqqmx * S2x 
     1     + 4d0 * ( 1d0 - x ) + 2d0 * lnx * ( 1d0 + x ) )
*
* ...The soft (`+'-distribution) part of the splitting function
*
       A2T = - 40D0/9D0*cf*nf + 268D0/9D0*ca*cf - 8D0*zeta2*ca*cf
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
      include "../commons/ColorFactors.h"
*
      INTEGER NF
*
* ..The soft coefficient for use in X2NSBT and X2NSCT
*
      COMMON / P1SOFTT / A2T
*
* ...some abbreviations
*
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
     3     / 6d0 + lnx**2 / 2d0 - pi**2 / 6d0 ) * pqq
     4     + ( 20d0 * ( 1d0 - x ) ) / 3d0 + lnx * ( 1d0 + x) )
     5     + 4d0 * CF**2 * ( ( ( 3d0 * lnx ) / 2d0 + 2d0 * ln1mx 
     6     * lnx - 2d0 * lnx**2 ) * pqq - 5d0 * ( 1d0 - x ) 
     7     + ( lnx**2 * ( 1d0 + x ) ) / 2d0 
     8     - lnx * ( 3.5d0 + ( 3d0 * x ) / 2d0 ) )
     9     - 4d0 * CF * ( CF - CA / 2d0 ) * ( 2d0 * pqqmx * S2x 
     1     + 4d0 * ( 1d0 - x ) + 2d0 * lnx * ( 1d0 + x ) )
*
* ...The soft (`+'-distribution) part of the splitting function
*
       A2T = - 40D0/9D0*cf*nf + 268D0/9D0*ca*cf - 8D0*zeta2*ca*cf
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
*
       include "../commons/consts.h"
       include "../commons/ColorFactors.h"
*
       INTEGER NF
*
* ...The coefficient of delta(1-x)
*
       P1DELT = 
     &     - 1D0/3D0*cf*nf
     &     + 3D0/2D0*cf**2
     &     + 17D0/6D0*ca*cf
     &     + 24D0*zeta3*cf**2
     &     - 12D0*zeta3*ca*cf
     &     - 8D0/3D0*zeta2*cf*nf
     &     - 12D0*zeta2*cf**2
     &     + 44D0/3D0*zeta2*ca*cf
*
* ...The soft (`+'-distribution) part of the splitting function
*
       A2T = - 40D0/9D0*cf*nf + 268D0/9D0*ca*cf - 8D0*zeta2*ca*cf
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
*
      include "../commons/ColorFactors.h"
*
      INTEGER NF
*
* ...some abbreviations
*
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
      include "../commons/ColorFactors.h"
*
      INTEGER NF
*
      lnx = dlog(x)
      ln1mx = dlog(1d0 - x)
      pqg   = x**2 + ( 1d0 - x )**2
      pqgmx = x**2 + ( 1d0 + x )**2
      S1x   = - ddilog(1d0-x)
      S2x   = S2(x)
*
      X1QGTA = ( NF**2 * ( - 8d0 / 3d0 - ( 16d0 / 9d0 
     1       + 8d0 * lnx / 3d0 + 8d0 * ln1mx / 3d0 ) * pqg )
     2       + 2d0 * CF * NF * ( - 2d0 + 3d0 * x 
     3       + ( - 7d0 + 8d0 * x ) * lnx - 4d0 * ln1mx 
     4       + ( 1d0 - 2d0 * x ) * lnx**2 
     5       + ( - 2d0 * ( lnx + ln1mx )**2 - 2d0 * ( ln1mx - lnx )
     6       + 16d0 * S1x + 2d0 * pi**2 - 10d0 ) * pqg )
     7       + 2d0 * CA * NF * ( - 152d0 / 9d0 + 166d0 * x / 9d0 
     8       - 40d0 / 9d0 / x + ( - 4d0 / 3d0 - 76d0 * x / 3d0 ) * lnx 
     9       + 4d0 * ln1mx + ( 2d0 + 8d0 * x ) * lnx**2 
     1       + ( 8d0 * lnx * ln1mx - lnx**2 - 4d0 * lnx / 3d0 
     2       + 10d0 * ln1mx / 3d0 + 2d0 * ln1mx**2 - 16d0 * S1x 
     3       - 7d0 * pi**2 / 3d0 + 178d0 / 9d0 ) * pqg 
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
      include "../commons/ColorFactors.h"
*
      INTEGER NF
*
* ...some abbreviations
*
      lnx = dlog(x)
      ln1mx = dlog(1d0 - x)
      pgq   = ( 1d0 + ( 1d0 - x )**2 ) / x
      pgqmx = - ( 1d0 + ( 1d0 + x )**2 ) / x
      S1x   = - ddilog(1d0-x)
      S2x   = S2(x)
*
      X1GQTA = 2d0 * NF * ( 4d0 * CF**2 * ( - 1d0 / 2d0 + 9d0 * x 
     1       / 2d0 + ( - 8d0 + x / 2d0 ) * lnx + 2d0 * x * ln1mx 
     2       + ( 1d0 - x / 2d0 ) * lnx**2 + ( ln1mx**2 
     3       + 4d0 * lnx * ln1mx - 8d0 * S1x 
     4       - 4d0 * pi**2 / 3d0 ) * pgq )
     5       + 4d0 * CF * CA * ( 62d0 / 9d0 - 35d0 * x / 18d0 
     6       - 44d0 * x**2 / 9d0 
     7       + ( 2d0 + 12d0 * x + 8d0 * x**2 / 3d0 ) * lnx 
     8       - 2d0 * x * ln1mx - ( 4d0 + x ) * lnx**2 + pgqmx * S2x 
     9       + ( - 2d0 * lnx * ln1mx - 3d0 * lnx - 3d0 * lnx**2 / 2d0 
     1       - ln1mx**2 + 8d0 * S1x + 7d0 * pi**2 / 6d0 
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
      include "../commons/consts.h"
      include "../commons/ColorFactors.h"
*
      INTEGER NF
*
* ..The soft coefficient for use in X1GGB and X1GGC
*
       COMMON / P1GSOFTT / A2GT
*
* ...some abbreviations
*
      lnx = dlog(x)
      ln1mx = dlog(1d0 - x)
      pgg   = ( 1d0 / ( 1d0 - x ) +  1d0 / x - 2d0 + x * ( 1d0 - x ) )
      pggmx = ( 1d0 / ( 1d0 + x ) -  1d0 / x - 2d0 - x * ( 1d0 + x ) )
      S2x   = S2(x)
      DM    = 1D0/(1D0-X)
*
      ggg1  = 2d0 * CF * NF * ( - 4d0 + 12d0 * x - 164d0 * x**2 / 9d0 
     1      + ( 10d0 + 14d0 * x + 16d0 * x**2 / 3d0 + 16d0 / 3d0 / x )
     2      * lnx + 92d0 / 9d0 / x + 2d0 * ( 1d0 + x ) * lnx**2 )
     3      + 2d0 * CA * NF * ( 2d0 - 2d0 * x 
     4      + 26d0 * ( x**2 - 1d0 / x ) / 9d0 
     5      - 4d0 * ( 1d0 + x ) * lnx / 3d0 
     6      - ( 20d0 / 9d0 + 8d0 * lnx / 3d0 ) * pgg )
     7      + 4d0 * CA * CA * ( 27d0 * ( 1d0 - x ) / 2d0 
     8      + 67d0 * ( x**2 - 1d0 / x ) / 9d0
     9      + ( 11d0 / 3d0 - 25d0 * x / 3d0 - 44d0 / 3d0 / x ) * lnx
     1      - 4d0 * ( 1d0 + x ) * lnx**2
     2      + ( 4d0 * lnx * ln1mx - 3d0 * lnx**2 + 22d0 * lnx / 3d0
     3      - 2d0 * zeta2 + 67d0 / 9d0 ) * pgg + 2d0 * pggmx * S2x )
*
* ...The soft (`+'-distribution) part of the splitting function
*
       A2GT = - 40D0/9D0*ca*nf + 268D0/9D0*ca**2 - 8D0*zeta2*ca**2
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
*
       include "../commons/consts.h"
       include "../commons/ColorFactors.h"
*
       INTEGER NF
*
* ...The coefficient of delta(1-x)
*
       P1DELT = 
     ,    - 2D0*cf*nf
     ,    - 8D0/3D0*ca*nf
     ,    + 32D0/3D0*ca**2
     ,    + 12D0*zeta3*ca**2
*
* ...The soft (`+'-distribution) part of the splitting function
*
       A2GT = - 40D0/9D0*ca*nf + 268D0/9D0*ca**2 - 8D0*zeta2*ca**2
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
       P2QGTA = ( P2QG1 + NF * P2QG2 + NF**2 * P2QG3 ) / 2D0
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
       P2GQTA = 2D0 * NF * ( P2GQ0 + NF * P2GQ1 )
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
************************************************************************
*
*     Space-like polarized splitting functions.
*     References:
*     - Vogt website: http://www.liv.ac.uk/~avogt/split.html
*     - hep-ph/9603366
*
***********************************************************************
*
* ..This is the regular 1-loop piece. 
*
       FUNCTION X0NSPA (X)
       IMPLICIT REAL*8 (A - Z)
*
       include "../commons/ColorFactors.h"
*
       X0NSPA = - 2D0 * CF * ( 1D0 + X )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the singular (soft) piece.
*
       FUNCTION X0NSPB (Y)
       IMPLICIT REAL*8 (A - Z)
*
       include "../commons/ColorFactors.h"
*
       X0NSPB = 4D0 * CF / ( 1D0 - Y )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' piece.
*
       FUNCTION X0NSPC (Y)
       IMPLICIT REAL*8 (A - Z)
*
       include "../commons/ColorFactors.h"
*
       X0NSPC = 4D0 * CF * LOG ( 1D0 - Y ) + 3D0 * CF
*
       RETURN
       END
*
* =====================================================================
*
*
* ..The 1-loop gluon->quark splitting functions P_qg^(0)
*
       FUNCTION X0QGPA (X, NF)
*
       IMPLICIT REAL*8 (A - Z)
*
       INTEGER NF
*
       X0QGPA = 2D0 * NF * ( 2D0 * X - 1D0 )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The 1-loop quark->gluon splitting functions P_gq^(0)
*
       FUNCTION X0GQPA (X)
*
       IMPLICIT REAL*8 (A - Z)
*
       include "../commons/ColorFactors.h"
*
       X0GQPA = 2D0 * CF * ( 2D0 - X )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The regular piece of the 1-loop gg splitting function P_gg^(0)
*
       FUNCTION X0GGPA (X)
       IMPLICIT REAL*8 (A - Z)
*
       include "../commons/ColorFactors.h"
*
       X0GGPA = 4D0 * CA * ( - 2D0 * X + 1D0 )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the singular (soft) piece.
*
       FUNCTION X0GGPB (X)
       IMPLICIT REAL*8 (A - Z)
*
       include "../commons/ColorFactors.h"
*
       X0GGPB = 4D0 * CA / ( 1D0 - X )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' piece.
*
       FUNCTION X0GGPC (X,NF)
       IMPLICIT REAL*8 (A - Z)
*
       include "../commons/ColorFactors.h"
*
       INTEGER NF
*
       X0GGPC = 4D0 * CA * LOG ( 1D0 - X )
     1      - 2D0 / 3D0 * NF + 11D0 / 3D0 * CA
*
       RETURN
       END
*
* =====================================================================
*
*
* ..This is the regular 2-loop piece for P_NS^-. 
*
      FUNCTION X1NSMPA (X, NF)
*
      IMPLICIT REAL*8 (A - Z)
*
      include "../commons/consts.h"
      include "../commons/ColorFactors.h"
*
      INTEGER NF
*
* ..The soft coefficient for use in X2NSB and X2NSC
*
      COMMON / P1SOFTP / A2P
*
* ...some abbreviations
*
      lnx   = dlog(x)
      ln1mx = dlog(1d0 - x)
      pqq   = 2d0 / ( 1d0 - x ) - 1d0 - x
      pqqmx = 2d0 / ( 1d0 + x ) - 1d0 + x
      S2x   = S2(x)
      DM    = 1d0 / ( 1d0 - x )
*
      gqq1 = 2d0 * CF * NF * ( ( - 1.1111111111111112d0 - ( 2d0 * lnx ) 
     1     / 3d0 ) * pqq - ( 4d0 * ( 1d0 - x ) ) / 3d0 ) 
     2     + 4d0 * CA * CF * ( ( 3.7222222222222223d0 + ( 11d0 * lnx ) 
     3     / 6d0 + lnx**2 / 2d0 - pi**2 / 6d0 ) * pqq
     4     + ( 20d0 * ( 1d0 - x ) ) / 3d0 + lnx * ( 1d0 + x) )
     5     + 4d0 * CF**2 * ( ( ( - 3d0 * lnx ) / 2d0 - 2d0 * ln1mx 
     6     * lnx ) * pqq - 5d0 * ( 1d0 - x ) - ( lnx**2 * ( 1d0 
     7     + x ) ) / 2d0 - lnx * ( 1.5d0 + ( 7d0 * x ) / 2d0 ) )
     8     + 4d0 * CF * ( CF - CA / 2d0 ) * ( 2d0 * pqqmx * S2x 
     9     + 4d0 * ( 1d0 - x ) + 2d0 * lnx * ( 1d0 + x ) )
*
* ...The soft (`+'-distribution) part of the splitting function
*
       A2P = - 40D0/9D0*cf*nf + 268D0/9D0*ca*cf - 8D0*zeta2*ca*cf
*
       GQQ1L = DM * A2P
*
* ...The regular piece of the coefficient function
*
       X1NSMPA = GQQ1 - GQQ1L
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the regular 2-loop piece for P_NS^+. 
*
      FUNCTION X1NSPPA (X, NF)
*
      IMPLICIT REAL*8 (A - Z)
*
      include "../commons/consts.h"
      include "../commons/ColorFactors.h"
*
      INTEGER NF
*
* ..The soft coefficient for use in X2NSB and X2NSC
*
      COMMON / P1SOFTP / A2P
*
* ...some abbreviations
*
      lnx   = dlog(x)
      ln1mx = dlog(1d0 - x)
      pqq   = 2d0 / ( 1d0 - x ) - 1d0 - x
      pqqmx = 2d0 / ( 1d0 + x ) - 1d0 + x
      S2x   = S2(x)
      DM    = 1d0 / ( 1d0 - x )
*
      gqq1 = 2d0 * CF * NF * ( ( - 1.1111111111111112d0 - ( 2d0 * lnx ) 
     1     / 3d0 ) * pqq - ( 4d0 * ( 1d0 - x ) ) / 3d0 ) 
     2     + 4d0 * CA * CF * ( ( 3.7222222222222223d0 + ( 11d0 * lnx ) 
     3     / 6d0 + lnx**2 / 2d0 - pi**2 / 6d0 ) * pqq
     4     + ( 20d0 * ( 1d0 - x ) ) / 3d0 + lnx * ( 1d0 + x) )
     5     + 4d0 * CF**2 * ( ( ( - 3d0 * lnx ) / 2d0 - 2d0 * ln1mx 
     6     * lnx ) * pqq - 5d0 * ( 1d0 - x ) - ( lnx**2 * ( 1d0 
     7     + x ) ) / 2d0 - lnx * ( 1.5d0 + ( 7d0 * x ) / 2d0 ) )
     8     - 4d0 * CF * ( CF - CA / 2d0 ) * ( 2d0 * pqqmx * S2x 
     9     + 4d0 * ( 1d0 - x ) + 2d0 * lnx * ( 1d0 + x ) )
*
* ...The soft (`+'-distribution) part of the splitting function
*
       A2P = - 40D0/9D0*cf*nf + 268D0/9D0*ca*cf - 8D0*zeta2*ca*cf
*
       GQQ1L = DM * A2P
*
* ...The regular piece of the coefficient function
*
       X1NSPPA = GQQ1 - GQQ1L
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the singular (soft) piece.
*
       FUNCTION X1NSPB (Y)
       IMPLICIT REAL*8 (A - Z)
*
       COMMON / P1SOFTP / A2P
*
       X1NSPB  = A2P / ( 1D0 - Y )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' piece.
*
       FUNCTION X1NSPC (Y, NF)
*
       IMPLICIT REAL*8 (A - Z)
*
       include "../commons/consts.h"
       include "../commons/ColorFactors.h"
*
       INTEGER NF
*
* ...The coefficient of delta(1-x)
*
       P1DELT = 
     &     - 1D0/3D0*cf*nf
     &     + 3D0/2D0*cf**2
     &     + 17D0/6D0*ca*cf
     &     + 24D0*zeta3*cf**2
     &     - 12D0*zeta3*ca*cf
     &     - 8D0/3D0*zeta2*cf*nf
     &     - 12D0*zeta2*cf**2
     &     + 44D0/3D0*zeta2*ca*cf
*
* ...The soft (`+'-distribution) part of the splitting function
*
       A2 = - 40D0/9D0*cf*nf + 268D0/9D0*ca*cf - 8D0*zeta2*ca*cf
*
       X1NSPC = LOG (1D0-Y) * A2 + P1DELT
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
* ..The 2-loop pure-singlet splitting functions P_ps^(1)
*
      function X1PSPA(x,nf)
*
      implicit none
*
      include "../commons/ColorFactors.h"
**
*     Input Variables
*
      integer nf
      double precision x
**
*     Internal Variables
*
      double precision lnx,lnx2
**
*     Output Variables
*
      double precision X1PSPA
*
*     Abbreviation
*
      lnx  = dlog(x)
      lnx2 = lnx * lnx
*
* ...The splitting function in terms of the harmonic polylogs
*
      X1PSPA = 8d0 * CF * TR * nf * ( ( 1d0 - x ) 
     1       - ( 1d0 - 3d0 * x ) * lnx - ( 1d0 + x ) * lnx2 )
*
      return
      end
*
* ---------------------------------------------------------------------
*
*
* ..The 2-loop gluon->quark splitting function DeltaP_qg^(1)
*
      FUNCTION X1QGPA (X, NF)
*
      IMPLICIT REAL*8 (A - Z)
*
      include "../commons/consts.h"
      include "../commons/ColorFactors.h"
*
      INTEGER NF
*
      lnx    = dlog(x)
      ln1mx  = dlog(1d0 - x)
      dpqg   =   2d0 * x - 1d0 
      dpqgmx = - 2d0 * x - 1d0
      S2x    = S2(x)
*
      X1QGPA =  4d0 * TR * nf * ( CF * ( - 22d0 + 27d0 * x - 9d0 * lnx
     1     + 8d0 * ( 1d0 - x ) * ln1mx + dpqg * ( 2d0 * ln1mx**2
     3     - 4d0 * ln1mx * lnx + lnx**2 - 4d0 * zeta2 ) )
     3     + CA * ( ( 24d0 - 22d0 * x )
     4     - 8d0 * ( 1d0 - x ) * ln1mx + ( 2d0 + 16d0 * x ) * lnx 
     5     - 2d0 * ( ln1mx**2 - zeta2 ) * dpqg
     6     - ( 2d0 * S2x - 3d0 * lnx**2 ) * dpqgmx ) )
*
      RETURN
      END
*
* ---------------------------------------------------------------------
*
*
* ..The 2-loop quark->gluon splitting function DeltaP_gq^(1)
*
      FUNCTION X1GQPA (X, NF)
*
      IMPLICIT REAL*8 (A - Z)
*
      include "../commons/consts.h"
      include "../commons/ColorFactors.h"
*
      INTEGER NF
*
* ...some abbreviations
*
      lnx    = dlog(x)
      ln1mx  = dlog(1d0 - x)
      dpgq   = 2d0 - x
      dpgqmx = 2d0 + x
      S2x    = S2(x)
*
      X1GQPA = 4d0 * ( CF * TR * nf * (- 4d0 / 9d0 * ( x + 4d0 )
     1     - 4d0 / 3d0 * dpgq * ln1mx )
     2     + CF * CF * ( - 1d0 / 2d0 - ( 4d0 - x ) * lnx / 2d0
     3     - dpgqmx * ln1mx + ( - 4d0 - ln1mx**2 + lnx**2 / 2d0 )
     4     * dpgq )
     5     + CF * CA * ( ( 4d0 - 13d0 * x ) * lnx + ( 10d0 + x )
     6     * ln1mx / 3d0 + ( 41d0 + 35d0 * x ) / 9d0
     7     + ( - 2d0 * S2x + 3d0 * lnx**2 ) * dpgqmx / 2d0
     8     + ( ln1mx**2 - 2d0 * ln1mx * lnx - zeta2 ) * dpgq ) )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..2-loop gluon-gluon splitting function DeltaP_gg^(1), regular part
*
      FUNCTION X1GGPA (X, NF)
*
      IMPLICIT REAL*8 (A - Z)
*
      include "../commons/consts.h"
      include "../commons/ColorFactors.h"
*
      INTEGER NF
*
* ..The soft coefficient for use in X1GGB and X1GGC
*
       COMMON / P1GSOFTP / A2GP
*
* ...some abbreviations
*
      lnx    = dlog(x)
      ln1mx  = dlog(1d0 - x)
      dpgg   = 1d0 / ( 1d0 - x ) - 2d0 * x + 1d0
      dpggmx = 1d0 / ( 1d0 + x ) + 2d0 * x + 1d0
      S2x    = S2(x)
      DM     = 1d0 / ( 1d0 - x )
*
      ggg1  = 4d0 * ( - CA * TR * nf * ( 4d0 * ( 1d0 - x )
     1     + 4d0 / 3d0 * ( 1d0 + x ) * lnx + 20d0 / 9d0 * dpgg )
     2     - CF * TR * nf * ( 10d0 * ( 1d0 - x )
     3     + 2d0 * ( 5d0 - x ) * lnx + 2d0 * ( 1d0 + x ) * lnx**2 )
     4     + CA * CA * ( ( 29d0 - 67d0 * x ) * lnx / 3d0
     5     - 19d0 * ( 1d0 - x ) / 2d0 + 4d0 * ( 1d0 + x ) * lnx**2
     4     - 2d0 * S2x * dpggmx + ( 67d0 / 9d0 - 4d0 * ln1mx * lnx
     5     + lnx**2 - 2d0 * zeta2 ) * dpgg ) )
*
* ...The soft (`+'-distribution) part of the splitting function
*
      A2GP = - 80d0 / 9d0 * CA * TR * nf + 268d0 / 9d0 * CA**2
     1     - 8d0 * zeta2 * CA**2
*
       GGG1L = DM * A2GP
*
* ...The regular piece of the coefficient function
*
       X1GGPA = GGG1 - GGG1L
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
* ..The singular (soft) piece of DeltaP_gg^(1)
*
       FUNCTION X1GGPB (Y, NF)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
*
       COMMON / P1GSOFTP / A2GP
*
       X1GGPB  = A2GP/(1.D0-Y)
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
* ..The 'local' piece of DeltaP_gg^(1)
*
       FUNCTION X1GGPC (Y, NF)
       IMPLICIT REAL*8 (A - Z)
*
       include "../commons/consts.h"
       include "../commons/ColorFactors.h"
*
       INTEGER NF
*
* ...The coefficient of delta(1-x)
*
       P1GDELT = 
     ,    - 4d0 * CF * TR * nf
     ,    - 16d0 / 3d0 * CA * TR * nf
     ,    + 32d0 / 3d0 * CA**2
     ,    + 12d0 * zeta3 * CA**2
*
* ...The soft (`+'-distribution) part of the splitting function
*
      A2GP = - 80d0 / 9d0 * CA * TR * nf + 268d0 / 9d0 * CA**2
     1     - 8d0 * zeta2 * CA**2
*
       X1GGPC = DLOG (1d0-Y) * A2GP + P1GDELT
*
       RETURN
       END
*
* =====================================================================
*
*
* ..The pure-singlet splitting function DeltaP_ps^(2). 
*    The quark-quark splitting function DeltaP_qq^(2) is obtained by 
*    adding the unpolarized quantity P_ns^(2)- given in hep-ph/0403192.

       FUNCTION P2PSPA (x, nf)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER nf
*
       L  = log (x)
       L1 = log (1.d0-x)
*
       P2ps1 = - 344./27.d0 * L**4 - (90.9198 + 81.50* x)* L**3 
     ,         - (368.6 - 349.9* x)* L*L - (739.0 - 232.57* L1)* L 
     ,         - 1362.6 + 1617.4 * x - 674.8 * x*x + 167.41 * x**3
     ,         - 204.76 * L1 - 12.61 * L1*L1 - 6.541 * L1**3
       P2ps2 = (1.1741 - 0.8253* x)* L**3  + (13.287 + 10.657* x)* L*L 
     ,         + 45.482 * L + 49.13 - 30.77 * x - 4.307 * x*x 
     ,         - 0.5094 *x**3 + 9.517 * L1 + 1.7805 * L1*L1
*
       P2PSPA = (1.-x) * nf * ( P2ps1 + nf * P2ps2 ) 
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The gluon->quark splitting function DeltaP_qg^(2).
*
       FUNCTION P2QGPA (x, nf)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER nf 
*
       L  = log (x)
       L1 = log (1.d0-x)
*
       P2qg1 = - 151./3.d0 * L**4 - (385.64 + 73.30* x)* L**3
     ,         - (894.8 - 1145.3* x)* L*L - (1461.2 - 825.4* L1)* L
     ,         - 2972.4 + 4672.* x - 1221.6 * x*x - 18.0 * x**3
     ,         + 278.32* L1 - 90.26* L1*L1 - 5.30* L1**3 + 3.784*L1**4
       P2qg2 =   16./9.d0 * L**4 + (30.739  + 10.186* x) * L**3 
     ,         + (196.96 + 179.1* x)* L*L + (526.3  - 47.30* L1)* L 
     ,         + 499.65 - 432.18 * x - 141.63 * x*x - 11.34 * x**3 
     ,         - 6.256 * L1 + 7.32 * L1*L1 + 0.7374 * L1**3
* 
       P2QGPA = nf * ( P2qg1 + nf * P2qg2 )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The quark->gluon splitting functions DeltaP_gq^(2). P2gq2 is exact.
*
       FUNCTION P2GQPA (x, nf)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER nf 
*
       L  = log (x)
       L1 = log (1.d0-x)
*
       P2gq0 =   11512./81.d0* L**4 + (888.003  + 175.1* x)* L**3
     ,         + (2140. - 850.7* x)* L*L + (4046.6 - 1424.8* L1)* L
     ,         + 6159. - 3825.9 * x + 1942.* x*x - 742.1 * x**3
     ,         + 1843.7* L1 + 451.55* L1*L1 + 59.3* L1**3 + 5.143* L1**4
       P2gq1 = - 128./27.d0 * L**4 - (39.3872 + 30.023*x)* L**3
     ,         - (202.46 + 126.53* x)* L*L - (308.98 + 16.18* L1)* L
     ,         - 301.07 - 296.0 * x + 406.13 * x*x - 101.62 * x**3
     ,         - 171.78* L1 - 47.86 * L1*L1 - 4.963 * L1**3
       P2gq2 =   16./27.d0 * ( - 12. + 10.* x + ( 8.+ 2.*x)* L1 
     ,         + (6.- 3.*x)* L1*L1 )
*
       P2GQPA = P2gq0 + nf * (P2gq1 + nf * P2gq2) 
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The gluon-gluon splitting function DeltaP_gg^(2), regular piece.
*
       FUNCTION P2GGPA (x, nf)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       L  = log (x)
       L1 = log (1.d0-x)
*
       P2ggA0 =   504.d0 * L**4 + (3777.5  + 1167.* x)* L**3
     ,          + (10902. - 863.* x)* L*L + (23091. - 12292.* L1)* L
     ,          + 30988. - 39925.* x + 13447.* x*x - 4576.* x**3 
     ,          - 13247.* (1.-x)*L1 + 3801.* L1
       P2ggA1 = - 766./27.d0 * L**4 - (357.798 - 131.* x)* L**3
     ,          - (1877.2 - 613.1* x)* L*L - (3524. + 7932.* L1)* L
     ,          - 1173.5 + 2648.6 * x - 2160.8 * x*x + 1251.7 * x**3
     ,          - 6746.* (1.-x)*L1 - 295.7* L1
       P2ggA2 = - 1.1809 * L**3 - (6.679 - 15.764* x)* L*L 
     ,          - (13.29 + 16.944* L1) * L - 16.606 + 32.905 * x 
     ,          - 18.30 * x*x + 2.637 * x**3 - 0.210 * L1
*
       P2GGPA = P2ggA0 + nf * ( P2ggA1 + nf * P2ggA2 )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The gluon-gluon splitting function DeltaP_gg^(2), singular piece.
*
       FUNCTION P2GGPB (x, nf)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER nf
*
       P2GGPB = ( 2643.521 - nf * (412.172 + nf * 16./9.D0 ) ) 
     1           / (1.d0-x)
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The gluon-gluon splitting function DeltaP_gg^(2), `local' piece.  
*
       FUNCTION P2GGPC (x, nf)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER nf
*
       L1 = log (1.d0-x)
*
       P2GGPC =             2643.521 * L1 + 4425.448 + 2.314   
     ,           - nf *    (  412.172 * L1 +  528.720 - 0.184 ) 
     ,           - nf*nf * ( 16./9.D0 * L1 -   6.4630 + 0.0023 )
*
       RETURN
       END
*
* =====================================================================
*
*
       FUNCTION P2NSMPA (Y, NF)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
       D81 = 1./81.D0
*
       P2NSMPA =   1641.1 - 3135.* Y + 243.6 * Y**2 - 522.1 * Y**3
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
       FUNCTION P2NSPPA (Y, NF)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
       D81 = 1./81.D0
*
       P2NSPPA =   1860.2 - 3505.* Y + 297.0 * Y**2 - 433.2 * Y**3
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
       FUNCTION P2NSPB (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       P2NSPB = ( 1174.898 - NF * 183.187 - NF**2 * 64./81.D0 ) / (1.-Y)
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
       FUNCTION P2NSMPC (Y, NF)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
*
       DL1 = LOG (1.-Y)
*
       P2NSMPC =       1174.898 * DL1 + 1295.624 - 0.24
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
       FUNCTION P2NSPPC (Y, NF)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
*
       DL1 = LOG (1.-Y)
*
       P2NSPPC =       1174.898 * DL1 + 1295.624 - 0.154
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
*     (For the polarized case, this contributions has not been calculated
*     yet. We include the function but we set it to zero.)
*
       FUNCTION P2NSSPA (Y, NF)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       P2NSSPA = 0D0
*
       RETURN
       END
*
***********************************************************************
*
*     O(alphas alpha) splitting functions for the NLO QED evolution
*
***********************************************************************
*
*
* ..This is the regular 2-loop piece for P_NS^+. 
*
      FUNCTION X1NSPA_ASA (X)
*
      IMPLICIT REAL*8 (A - Z)
*
      include "../commons/ColorFactors.h"
*
* ...some abbreviations
*
      lnx   = dlog(x)
      ln1mx = dlog(1d0 - x)
      pqq   = 2d0 / ( 1d0 - x ) - 1d0 - x
      pqqmx = 2d0 / ( 1d0 + x ) - 1d0 + x
      S2x   = S2(x)
*
      X1NSPA_ASA = 8d0 * CF * ( ( ( - 3d0 * lnx ) / 2d0 - 2d0 * ln1mx 
     1     * lnx ) * pqq - 5d0 * ( 1d0 - x ) - ( lnx**2 * ( 1d0 
     2     + x ) ) / 2d0 - lnx * ( 1.5d0 + ( 7d0 * x ) / 2d0 ) )
     3     + 8d0 * CF * ( 2d0 * pqqmx * S2x 
     4     + 4d0 * ( 1d0 - x ) + 2d0 * lnx * ( 1d0 + x ) )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the regular 2-loop piece for P_NS^-. 
*
      FUNCTION X1NSMA_ASA (X)
*
      IMPLICIT REAL*8 (A - Z)
*
      include "../commons/ColorFactors.h"
*
* ...some abbreviations
*
      lnx   = dlog(x)
      ln1mx = dlog(1d0 - x)
      pqq   = 2d0 / ( 1d0 - x ) - 1d0 - x
      pqqmx = 2d0 / ( 1d0 + x ) - 1d0 + x
      S2x   = S2(x)
*
      X1NSMA_ASA = 8d0 * CF * ( ( ( - 3d0 * lnx ) / 2d0 - 2d0 * ln1mx 
     1     * lnx ) * pqq - 5d0 * ( 1d0 - x ) - ( lnx**2 * ( 1d0 
     2     + x ) ) / 2d0 - lnx * ( 1.5d0 + ( 7d0 * x ) / 2d0 ) )
     3     - 8d0 * CF * ( 2d0 * pqqmx * S2x 
     4     + 4d0 * ( 1d0 - x ) + 2d0 * lnx * ( 1d0 + x ) )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' piece.
*
       FUNCTION X1NSC_ASA ()
*
       IMPLICIT REAL*8 (A - Z)
*
       include "../commons/consts.h"
       include "../commons/ColorFactors.h"
*
* ...The coefficient of delta(1-x)
*
       X1NSC_ASA = 2d0 * cf * (
     &     + 3D0/2D0
     &     + 24D0*zeta3
     &     - 12D0*zeta2 )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The 2-loop gluon->quark splitting functions P_qgamma^(1,1)
*
      FUNCTION X1QGAMA_ASA (X)
*
      IMPLICIT REAL*8 (A - Z)
*
      include "../commons/consts.h"
      include "../commons/ColorFactors.h"
*
      lnx   = dlog(x)
      ln1mx = dlog(1d0 - x)
      pqg   = x**2 + ( 1d0 - x )**2
*
      X1QGAMA_ASA = 2d0 * CF * ( 4d0 + 4d0 * ln1mx + ( 10d0 - 4d0 
     1      * ( ln1mx - lnx ) + 2d0 * ( - ln1mx + lnx )**2 
     2      - 2d0 * pi**2 / 3d0 ) * pqg - lnx * ( 1d0 - 4d0 * x )
     3      - lnx**2 * ( 1d0  - 2d0 * x ) - 9d0 * x )
*
      RETURN
      END
*
* ---------------------------------------------------------------------
*
*
* ..The 2-loop quark->gluon splitting functions P_gq^(1,1)
*
      FUNCTION X1GAMQA_ASA (X)
*
      IMPLICIT REAL*8 (A - Z)
*
      include "../commons/ColorFactors.h"
*
* ...some abbreviations
*
      lnx = dlog(x)
      ln1mx = dlog(1d0 - x)
      pgq   = ( 1d0 + ( 1d0 - x )**2 ) / x
*
      X1GAMQA_ASA = 4d0 * CF * ( - 2.5d0 - ( 3d0 * ln1mx + ln1mx**2 )
     1      * pgq - lnx**2 * ( 1d0 - x / 2d0 ) - ( 7d0 * x ) / 2d0 
     2      - 2d0 * ln1mx * x + lnx * ( 2d0 + ( 7d0 * x ) / 2d0 ) )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The regular piece of the 2-loop gg splitting function P_ggamma^(1,1) 
*
      FUNCTION X1GGAMA_ASA (X)
*
      IMPLICIT REAL*8 (A - Z)
*
      include "../commons/ColorFactors.h"
*
* ...some abbreviations
*
      lnx = dlog(x)
      ln1mx = dlog(1d0 - x)
*
      X1GGAMA_ASA  = 4d0 * CF *( - 16d0 + 4d0 / ( 3d0 * x ) + 8d0 * x 
     1      + ( 20d0 * x**2 ) / 3d0 - lnx**2 * ( 2d0 + 2d0 * x ) 
     2      - lnx * ( 6d0 + 10d0 * x ) )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' piece.
*
       FUNCTION X1GAMGAMC_ASA ()
*
       IMPLICIT REAL*8 (A - Z)
*
       include "../commons/ColorFactors.h"
*
* ...The coefficient of delta(1-x)
*
       X1GAMGAMC_ASA = - 4D0 * CF
*
       RETURN
       END
*
***********************************************************************
*
*     O(alpha^2) splitting functions for the NLO QED evolution
*
***********************************************************************
*
*
* ..The 2-loop pure-singlet splitting functions P_ps^(0,2)
*
      FUNCTION X1PSA_AA (X)
*
      IMPLICIT REAL*8 (A - Z)
*
* ...some abbreviations
*
      DX = 1D0/X
      LNX = DLOG(X)
      HR200 = LNX * LNX / 2D0
*
* ...The splitting function in terms of the harmonic polylogs
*
      X1PSA_AA =
     &    ( - 8D0 + 24D0*x - 224D0/9D0*x**2 + 80D0/9D0*
     &    dx + 4D0*LNX + 20D0*LNX*x + 32D0/3D0*LNX*x**2 -
     &    8D0*HR200 - 8D0*HR200*x )
*
      RETURN
      END
*
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' piece of  P_gammagamma^(0,2).
*
       FUNCTION X1GAMGAMC_AA ()
*
       IMPLICIT REAL*8 (A - Z)
*
* ...The coefficient of delta(1-x)
*
       X1GAMGAMC_AA = - 4D0
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The 2-loop remainder term to the quark->gluon splitting functions P_gammaq^(0,2)
*
      FUNCTION REM_X1GAMQA_AA (X)
*
      IMPLICIT REAL*8 (A - Z)
*
* ...some abbreviations
*
      ln1mx = dlog(1d0 - x)
      pgq   = ( 1d0 + ( 1d0 - x )**2 ) / x
*
      REM_X1GAMQA_AA = 4d0 * ( - 4d0 * x / 3d0
     1               - pgq * ( 20d0 / 9d0 + 4d0 * ln1mx / 3d0 ) )
*
       RETURN
       END
*
* =====================================================================
*
*
* ..This is the  remainder term to the regular 2-loop piece for P_NS^\pm(0,2). 
*
      FUNCTION REM_X1NSA_AA (X)
*
      IMPLICIT REAL*8 (A - Z)
*
* ..The soft coefficient for use in X2NSB and X2NSC
*
      COMMON / REM_P1SOFT_AA / A2
*
* ...some abbreviations
*
      lnx = dlog(x)
      pqq = 2d0 / ( 1d0 - x ) - 1d0 - x
      DM  = 1d0 / ( 1d0 - x )
*
      gqq1 = 4d0 * ( ( - 1.1111111111111112d0 - ( 2d0 * lnx ) 
     1     / 3d0 ) * pqq - ( 4d0 * ( 1d0 - x ) ) / 3d0 ) 
*
* ...The soft (`+'-distribution) part of the splitting function
*
      A2 = - 80D0/9D0
*
      GQQ1L = DM * A2
*
* ...The regular piece of the coefficient function
*
      REM_X1NSA_AA = GQQ1 - GQQ1L
*
      RETURN
      END
*
* ---------------------------------------------------------------------
*
*
* ..This is the singular (soft) piece.
*
      FUNCTION REM_X1NSB_AA (Y)
      IMPLICIT REAL*8 (A - Z)
*
      COMMON / REM_P1SOFT_AA / A2
*
      REM_X1NSB_AA  = A2 / ( 1d0 - Y )
*
      RETURN
      END
*
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' piece.
*
       FUNCTION REM_X1NSC_AA (Y)
*
       IMPLICIT REAL*8 (A - Z)
*
       include "../commons/consts.h"
*
* ...The coefficient of delta(1-x)
*
       P1DELT = 
     &     - 2D0/3D0
     &     - 16D0/3D0*zeta2
*
* ...The soft (`+'-distribution) part of the splitting function
*
       A2 = - 80D0/9D0
*
       REM_X1NSC_AA = LOG (1D0-Y) * A2 + P1DELT
*
       RETURN
       END
