*
* ..File: xcttsg2e.f    FT_{ps,g}  in  e+e- -> h + X
*
*
* ..The 2-loop MS(bar) pure-singlet and gluon coefficient functions for
*    the fragmentation function F_T in e+e- -> h + X at mu_r = mu_f = Q.
*    Expansion parameter: alpha_s/(4 pi).
*
* ..These functions are regular, i.e., there is only an `A' piece in 
*    the notation of Appendix B of Floratos, Kounnas, Lacaze (1981).
*
* ..The code uses the package of Gehrmann and Remiddi for the harmonic
*    polylogarithms published in hep-ph/0107173 = CPC 141 (2001) 296.
*
* ..First calculation:
*    P. Rijken and W. van Neerven, hep-ph/9609377 = NP B488 (1997) 233
*
*
* =====================================================================
*
*
* ..Pure singlet
*
       FUNCTION CTTPS2A (X, NF)
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
       PARAMETER ( Z2 = 1.6449 34066 84822 64365 D0,
     ,             Z3 = 1.2020 56903 15959 42854 D0 )
*
* ...Colour factors and abbreviations
*
       CF  = 4./3.D0
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
       CTTPS2A = cTeqps2
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..Gluon
*
       FUNCTION CTTGL2A (X, NF)
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
       PARAMETER ( Z2 = 1.6449 34066 84822 64365 D0,
     ,             Z3 = 1.2020 56903 15959 42854 D0 )
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
       CTTGL2A = cTeg2
*
       RETURN
       END
*
* =================================================================av==
