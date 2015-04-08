*
* ..File: xcltsg2e.f    FL_{ps,g}  in  e+e- -> h + X
*
*
* ..The 2-loop MS(bar) pure-singlet and gluon coefficient functions for
*    the fragmentation function F_L in e+e- -> h + X at mu_r = mu_f = Q.
*    Expansion parameter: alpha_s/(4 pi).
*
* ..These functions are regular, i.e., there is only an `A' piece in 
*    the notation of Appendix B of Floratos, Kounnas, Lacaze (1981).
*
* ..The code uses the package of Gehrmann and Remiddi for the harmonic
*    polylogarithms published in hep-ph/0107173 = CPC 141 (2001) 296.
*
* ..First calculation:
*    P. Rijken and W. van Neerven, hep-ph/9604436 = PL B386 (1996) 422
*
*
* =====================================================================
*
*
* ..Pure singlet
*
       FUNCTION CLTPS2A (X, NF)
*
       IMPLICIT REAL*8 (A - Z)
       COMPLEX*16 HC1, HC2, HC3, HC4, HC5
       INTEGER NF, NF2, N1, N2, NW
       PARAMETER ( N1 = -1, N2 = 1, NW = 2 )
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
*
       DX = 1.D0/X
*
* ...Harmonic polylogs (HPLs) up to weight NW=2 by Gehrmann and Remiddi 
*
       CALL HPLOG (X, NW, HC1,HC2,HC3,HC4, HR1,HR2,HR3,HR4,
     ,             HI1,HI2,HI3,HI4, N1, N2)
*
* ...The coefficient function in terms of the HPLs
*
       cLeqps2 =
     &  + nf*cf * (  - 56.D0/3.D0 + 104.D0/3.D0*x - 8*x**2 - 8*dx + 8*
     &    z2 - 16*Hr1(0) - 16*Hr1(0)*x + 8.D0/3.D0*Hr1(0)*x**2 + 32.D0/
     &    3.D0*Hr1(0)*dx + 8*Hr1(1)*x - 8.D0/3.D0*Hr1(1)*x**2 - 16.D0/3.
     &    D0*Hr1(1)*dx + 24*Hr2(0,0) - 8*Hr2(0,1) )
*
       CLTPS2A = cLeqps2
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..Gluon
*
       FUNCTION CLTGL2A (X, NF)
*
       IMPLICIT REAL*8 (A - Z)
       COMPLEX*16 HC1, HC2, HC3, HC4, HC5
       INTEGER NF, NF2, N1, N2, NW
       PARAMETER ( N1 = -1, N2 = 1, NW = 2 )
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
* ...Harmonic polylogs (HPLs) up to weight NW=2 by Gehrmann and Remiddi 
*
       CALL HPLOG (X, NW, HC1,HC2,HC3,HC4, HR1,HR2,HR3,HR4,
     ,             HI1,HI2,HI3,HI4, N1, N2)
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
       CLTGL2A = cLeg2
*
       RETURN
       END
*
* =================================================================av==
