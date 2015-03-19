*
* ..File: xcltns2e.f    FL_NS  in  e+e- -> h + X
*
*
* ..The exact 2-loop MS(bar) non-singlet coefficient functions for the
*    fragmentation function  F_L  in e+e- -> h + X at mu_r = mu_f = Q.
*    Expansion parameter: alpha_s/(4 pi).
*
* ..The distributions (in the mathematical sense) are given as in eq.
*    (B.26) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.
*    The name-endings A, B, and C of the functions below correspond to
*    the kernel superscripts [2], [3], and [1] in that equation.
*
* ..The code uses the package of Gehrmann and Remiddi for the harmonic
*    polylogarithms published in hep-ph/0107173 = CPC 141 (2001) 296.
*
* ..This function was first calculated in
*   P. Rijken and W. van Neerven,  hep-ph/9604436 = PL B386 (1996) 422
*
* ..Reference: A. Mitov, S. Moch and A. Vogt,  
*              hep-ph/0604053 = Phys. Lett. B638 (2006) 61
*
*
* =====================================================================
*
*
* ..There is only a regular piece.
*
       FUNCTION XLTNP2A (X, NF)
*
       IMPLICIT REAL*8 (A - Z)
       COMPLEX*16 HC1, HC2, HC3, HC4, HC5
       INTEGER NF, NF2, N1, N2, NW
       PARAMETER ( N1 = -1, N2 = 1, NW = 3 )
       DIMENSION HC1(N1:N2),HC2(N1:N2,N1:N2),HC3(N1:N2,N1:N2,N1:N2),
     ,           HC4(N1:N2,N1:N2,N1:N2,N1:N2) 
       DIMENSION HR1(N1:N2),HR2(N1:N2,N1:N2),HR3(N1:N2,N1:N2,N1:N2),
     ,           HR4(N1:N2,N1:N2,N1:N2,N1:N2) 
       DIMENSION HI1(N1:N2),HI2(N1:N2,N1:N2),HI3(N1:N2,N1:N2,N1:N2),
     ,           HI4(N1:N2,N1:N2,N1:N2,N1:N2) 
       PARAMETER ( Z2 = 1.6449 34066 84822 64365 D0 )
*
* ...Colour factors and abbreviations
*
       CF  = 4./3.D0
       CA  = 3.D0
*
       DX = 1.D0/X
*
* ...Harmonic polylogs (HPLs) up to weight NW=3 by Gehrmann and Remiddi 
*
       CALL HPLOG (X, NW, HC1,HC2,HC3,HC4, HR1,HR2,HR3,HR4,
     ,             HI1,HI2,HI3,HI4, N1, N2)
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
       XLTNP2A = CLeq2
*
       RETURN
       END
*
* =================================================================av==
