*
* ..File: xcttns2e.f    F_T,NS  for  e+e- annihilation
*
*
* ..The exact 2-loop MS(bar) non-singlet coefficient functions for the 
*    fragmentation function  F_T  in e+e- -> h + X at mu_r = mu_f = Q. 
*    Expansion parameter: alpha_s/(4 pi).
* 
* ..The distributions (in the mathematical sense) are given as in eq.
*    (B.26) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.
*    The name-endings A, B, and C of the functions below correspond to 
*    the kernel superscripts [2], [3], and [1] in that equation.
*    The regular piece  XTTNP2A  has to be called first.
*
* ..The code uses the package of Gehrmann and Remiddi for the harmonic
*    polylogarithms published in hep-ph/0107173 = CPC 141 (2001) 296.
*
* ..This function was first calculated in 
*   P. Rijken and W. van Neerven,  hep-ph/9609377 = NP B488 (1997) 233
*
* ..Reference: A. Mitov, S. Moch and A. Vogt,  
*              hep-ph/0604053 = Phys. Lett. B638 (2006) 61
*
*
* =====================================================================
*
*
* ..The regular piece. 
*
       FUNCTION XTTNP2A (X, NF)
*
       IMPLICIT REAL*8 (A - Z)
       COMPLEX*16 HC1, HC2, HC3, HC4 
       INTEGER NF, NF2, N1, N2, NW, I1, I2, I3, N
       PARAMETER ( N1 = -1, N2 = 1, NW = 3 ) 
       DIMENSION HC1(N1:N2),HC2(N1:N2,N1:N2),HC3(N1:N2,N1:N2,N1:N2), 
     ,           HC4(N1:N2,N1:N2,N1:N2,N1:N2) 
       DIMENSION HR1(N1:N2),HR2(N1:N2,N1:N2),HR3(N1:N2,N1:N2,N1:N2), 
     ,           HR4(N1:N2,N1:N2,N1:N2,N1:N2) 
       DIMENSION HI1(N1:N2),HI2(N1:N2,N1:N2),HI3(N1:N2,N1:N2,N1:N2), 
     ,           HI4(N1:N2,N1:N2,N1:N2,N1:N2) 
       PARAMETER ( Z2 = 1.6449 34066 84822 64365 D0,
     ,             Z3 = 1.2020 56903 15959 42854 D0 )
*
* ..The soft coefficients for use in X2NP2B and X2NP2C
*
       COMMON / C2TSOFT / A0, A1, A2, A3
*
* ...Colour factors
*
       CF  = 4./3.D0
       CA  = 3.D0
*
* ...Some abbreviations
*
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
       XTTNP2A = CTeq2 - CTeq2L
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The singular (soft) piece.
*
       FUNCTION XTTNS2B (Y, NF)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
*
       COMMON / C2TSOFT / A0, A1, A2, A3
*
       DL1 = LOG (1.D0-Y)
       DM  = 1.D0/(1.D0-Y)
*
       XTTNS2B = DM * ( DL1**3 * A3 + DL1**2 * A2 + DL1 * A1 + A0)
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The 'local' piece.
*
       FUNCTION XTTNS2C (Y, NF)
*
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF, NF2
       PARAMETER ( Z2 = 1.6449 34066 84822 64365 D0,
     ,             Z3 = 1.2020 56903 15959 42854 D0 )
*
       COMMON / C2TSOFT / A0, A1, A2, A3
*
* ...Colour factors
*
       CF  = 4./3.D0
       CA  = 3.D0
*
* ...The coefficient of delta(1-x)
*
       C2DELT = 
     ,    + ca*cf * ( - 5465.D0/72.D0 + 140.D0/3.D0*z3 + 215.D0/3.D0*z2
     ,                - 49.D0/5.D0*z2**2 )
     ,    + cf**2 * ( 331.D0/8.D0 - 78.D0*z3 - 39.D0*z2 + 30.D0*z2**2 )
     ,    + cf*nf * ( 457.D0/36.D0 + 4.D0/3.D0*z3 - 38.D0/3.D0*z2 )
*
       DL1 = LOG (1.D0-Y)
*
       XTTNS2C =   DL1**4 * A3/4.D0 + DL1**3 * A2/3.D0 
     ,           + DL1**2 * A1/2.D0 + DL1 * A0 + C2DELT
*
       RETURN
       END
*
* =================================================================av==
